import numpy as np
from utils.pint_registry import ureg
import equations

def enforce_quantity(val, unit):
    if hasattr(val, "units"):
        return val
    return ureg.Quantity(val, unit)

class GrowthProductionModule:
    """
    Handles microbial growth and acetate production for both biofilm and planktonic cells.
    Incorporates substrate-driven growth, yield-based product formation, maintenance costs, and biomass decay.
    """
    def __init__(self, registry):
        self.registry = registry
        self.vars = self._load_vars(registry)
        v = self.vars
        self.X_B = enforce_quantity(v['X_B'], 'mol')
        self.X_P = enforce_quantity(v['X_P'], 'mol/meter**3')  # Planktonic: concentration!
        self.acetate = 0.0 * ureg('mol')

    def _load_vars(self, registry):
        required_vars = {
            "V_catholyte": ("V_catholyte", "meter**3"),
            "X_B": ("X_B", "mol"),
            "X_P": ("X_P", "mol/meter**3"),
            "A_cat": ("A_cat", "meter**2"),
            "YAc_X_B": ("YAc_X_B", "mol/mol"),
            "YAc_X_P": ("YAc_X_P", "mol/mol"),
            "m_e": ("m_e", "1/second"),
            "m_H2": ("m_H2", "mol/mol/second"),
            "k_decay_B": ("k_decay_B", "1/second"),
            "k_decay_P": ("k_decay_P", "1/second"),
        }
        vars_out = {}
        missing = []
        for varname, (regkey, req_unit) in required_vars.items():
            try:
                value = registry.get(regkey)
                if value is None:
                    missing.append(regkey)
                    continue
                if req_unit and not hasattr(value, "to"):
                    if isinstance(value, str):
                        value = float(value)
                    value = ureg.Quantity(value, req_unit)
                elif req_unit:
                    value = value.to(req_unit)
                vars_out[varname] = value
            except Exception as e:
                missing.append(f"{regkey} ({str(e)})")
        if missing:
            raise ValueError("Missing variables in GrowthProductionModule: " + ", ".join(missing))
        return vars_out

    def step(
        self,
        dt,
        ne_DET,        # [mol e⁻ m⁻² s⁻¹] Electron flux to biofilm (from uptake module)
        nCO2_DET,      # [mol/s] CO2 uptake biofilm (from uptake)
        r_H2,          # [mol/s] H2 uptake planktonic (from uptake)
        nCO2_P,        # [mol/s] CO2 uptake planktonic (from uptake)
    ):
        v = self.vars

        dt = enforce_quantity(dt, 'second')
        ne_DET = enforce_quantity(ne_DET, 'mol / meter**2 / second')
        nCO2_DET = enforce_quantity(nCO2_DET, 'mol / second')
        r_H2 = enforce_quantity(r_H2, 'mol / second')
        nCO2_P = enforce_quantity(nCO2_P, 'mol / second')

        # Maintenance
        maintenance_e = v['m_e'] * self.X_B * dt  # [mol e-]
        maintenance_H2 = v['m_H2'] * self.X_P * v['V_catholyte'] * dt  # [mol]

        # Electron and H2 available for growth after maintenance
        avail_e = ne_DET * v['A_cat'] * dt - maintenance_e
        avail_e = max(0 * ureg('mol'), avail_e)
        e_in = avail_e / dt  # [mol e-/s]

        avail_H2 = r_H2 * dt - maintenance_H2
        avail_H2 = max(0 * ureg('mol'), avail_H2)
        H2_in = avail_H2 / dt  # [mol H2/s]

        # Growth rates
        r_growth_B = equations.r_growth_biofilm(
            ne_DET=e_in / v['A_cat'],
            nCO2_DET_dot=nCO2_DET,
            A_cat=v['A_cat']
        )  # [mol X / s]
        r_growth_P = equations.r_growth_plank(
            r_H2=H2_in,
            nCO2_P_dot=nCO2_P
        )  # [mol X / s]

        # Acetate production
        r_acetate_B = equations.r_acetate_biofilm(v['YAc_X_B'], r_growth_B)
        r_acetate_P = equations.r_acetate_plank(v['YAc_X_P'], r_growth_P)

        # Decay
        decay_B = v['k_decay_B'] * self.X_B * dt               # [mol]
        decay_P = v['k_decay_P'] * self.X_P * v['V_catholyte'] * dt   # [mol]

        # Biomass update (ensure non-negativity)
        dX_B = r_growth_B * dt - decay_B                       # [mol]
        dX_P = (r_growth_P * dt - decay_P) / v['V_catholyte']  # [mol/m³]
        self.X_B = max(0 * ureg('mol'), self.X_B + dX_B)
        self.X_P = max(0 * ureg('mol/meter**3'), self.X_P + dX_P)

        # Acetate accumulation
        d_acetate = (r_acetate_B + r_acetate_P) * dt
        self.acetate += d_acetate

        # Output all rates and new pools
        return {
            "r_growth_B": r_growth_B,
            "r_growth_P": r_growth_P,
            "r_acetate_B": r_acetate_B,
            "r_acetate_P": r_acetate_P,
            "decay_B": decay_B,
            "decay_P": decay_P,
            "dX_B": dX_B,
            "dX_P": dX_P,
            "X_B": self.X_B,
            "X_P": self.X_P,
            "d_acetate": d_acetate,
            "acetate": self.acetate,
            "maintenance_e": maintenance_e,
            "maintenance_H2": maintenance_H2,
            "avail_electrons": e_in,
            "avail_H2": H2_in
        }

    def set_biomass(self, X_B, X_P):
        self.X_B = enforce_quantity(X_B, 'mol')
        self.X_P = enforce_quantity(X_P, 'mol/meter**3')

    def get_state(self):
        return {
            "X_B": self.X_B,
            "X_P": self.X_P,
            "acetate": self.acetate
        }
