import numpy as np
from utils.pint_registry import ureg
import equations

def enforce_quantity(val, unit):
    """Ensure value is a Pint Quantity with the given unit."""
    if hasattr(val, "units"):
        return val
    return ureg.Quantity(val, unit)

class MicrobialUptakeModule:
    """
    Calculates microbial substrate uptake rates for biofilm (DET) and planktonic (PT) cells.
    Computes electron, CO2, and H2 fluxes using core equations from equations module.
    Supports ODE-style stepping to update uptake rates based on current system state and parameters.
    Uses stoichiometry: 
      - 1 mol CO2 : 4.2 mol e- (biofilm/DET)
      - 1 mol CO2 : 2.1 mol H2 (planktonic/H2)
    """

    def __init__(self, registry):
        self.registry = registry
        self.vars = self._load_vars(registry)
        v = self.vars
        # State variables, can be updated by simulation loop:
        self.X_B = enforce_quantity(v['X_B'], "mol")    # Biofilm biomass (mol, or gDW/L)
        self.X_P = enforce_quantity(v['X_P'], "mol / meter**3")  # Planktonic biomass (mol/m³, or gDW/L)
        self.A_biofilm = enforce_quantity(v['A_biofilm'], "meter**2")
        self.F = enforce_quantity(v['F'], "coulomb / mole")
        self.A_cat = enforce_quantity(v['A_cat'], "meter**2")
        self.f_DET = v['f_DET']  # fraction, unitless
        self.kL_CO2 = enforce_quantity(v['kL_CO2'], "meter / second")
        self.V_catholyte = enforce_quantity(v['V_catholyte'], "meter**3")
        self.kA_H2 = enforce_quantity(v['kA_H2'], "1 / second")

    def _load_vars(self, registry):
        required_vars = {
            # All as Pint quantities!
            "F": ("F", "coulomb / mole"),
            "A_cat": ("A_cat", "meter**2"),
            "f_DET": ("f_DET", ""),  # fraction, unitless
            "X_B": ("X_B", "mol"),
            "A_biofilm": ("A_biofilm", "meter**2"),
            "kL_CO2": ("kL_CO2", "meter / second"),
            "X_P": ("X_P", "mol / meter**3"),
            "V_catholyte": ("V_catholyte", "meter**3"),
            "kA_H2": ("kA_H2", "1 / second"),
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
            raise ValueError("Missing variables in MicrobialUptakeModule: " + ", ".join(missing))
        return vars_out

    def step(
        self,
        dt,
        I,  # applied current (ampere)
        C_CO2, C_CO2_sat,
        C_H2, C_H2_sat
    ):
        v = self.vars

        # Enforce units on all inputs
        dt = enforce_quantity(dt, "second")
        I = enforce_quantity(I, "ampere")
        C_CO2 = enforce_quantity(C_CO2, "mol / meter**3")
        C_CO2_sat = enforce_quantity(C_CO2_sat, "mol / meter**3")
        C_H2 = enforce_quantity(C_H2, "mol / meter**3")
        C_H2_sat = enforce_quantity(C_H2_sat, "mol / meter**3")

        # --- Biofilm (DET) limiting uptake: e- (I) ---
        n_electrons = max(0, self.f_DET * I / self.F)    # [mol/s]
        nCO2_det = max(0, n_electrons / 4.2)             # [mol/s]

        # --- Planktonic (H2-mediated) limiting uptake: H2 ---
        r_H2 = equations.r_H2(self.kA_H2, C_H2, C_H2_sat)  # [mol/m³/s]
        r_H2 = enforce_quantity(r_H2, "mol / meter**3 / second")
        r_H2 = max(0 * ureg('mol / meter**3 / second'), r_H2)

        n_H2 = max(0 * ureg('mol/second'), r_H2 * self.V_catholyte)
        nCO2_P = max(0 * ureg('mol/second'), n_H2 / 2.1)
        rCO2_P = max(0 * ureg('mol / meter**3 / second'), r_H2 / 2.1)

        # Specific rates (per planktonic biomass)
        X_P_safe = enforce_quantity(self.X_P, "mol / meter**3")

        if X_P_safe.magnitude > 0:
            q_H2_P = (r_H2 / X_P_safe).to('1/second')
            q_CO2_P = (rCO2_P / X_P_safe).to('1/second')
        else:
            q_H2_P = 0 * ureg('1/second')
            q_CO2_P = 0 * ureg('1/second')

        # Specific rates (useful for later growth calculations)
        X_B_safe = enforce_quantity(self.X_B, "mol")
        if X_B_safe.magnitude > 0:
            q_e_biofilm = (n_electrons / X_B_safe).to('1/second')
            q_CO2_B = (nCO2_det / X_B_safe).to('1/second')
        else:
            q_e_biofilm = 0 * ureg('1/second')
            q_CO2_B = 0 * ureg('1/second')

        # Return only Pint Quantities (never floats)
        uptake = {
            "n_electrons": enforce_quantity(n_electrons, "mol/second"),
            "nCO2_DET": enforce_quantity(nCO2_det, "mol/second"),
            "q_e_biofilm": enforce_quantity(q_e_biofilm, "1/second"),
            "q_CO2_B": enforce_quantity(q_CO2_B, "1/second"),
            "n_H2": enforce_quantity(n_H2, "mol/second"),
            "n_CO2_P": enforce_quantity(nCO2_P, "mol/second"),
            "q_H2_P": enforce_quantity(q_H2_P, "1/second"),
            "q_CO2_P": enforce_quantity(q_CO2_P, "1/second"),
        }
        return uptake

    def update_biomass(self, dX_B, dX_P):
        self.X_B += enforce_quantity(dX_B, "mol")
        self.X_P += enforce_quantity(dX_P, "mol / meter**3")

    def set_biomass(self, X_B, X_P):
        self.X_B = enforce_quantity(X_B, "mol")
        self.X_P = enforce_quantity(X_P, "mol / meter**3")

    def get_state(self):
        return {
            "X_B": self.X_B,
            "X_P": self.X_P
        }
