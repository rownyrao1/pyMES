import numpy as np
from utils.pint_registry import ureg
import equations

def enforce_quantity(val, unit):
    if hasattr(val, "units"):
        return val
    return ureg.Quantity(val, unit)

def clamp_nonnegative(val):
    # Clamps negative values to zero, keeps units
    return max(val, 0 * val.units) if hasattr(val, "units") else max(val, 0.0)

class GasDynamicsModule:
    """
    Models the gas phase and dissolved dynamics for H2 and CO2 in an MES reactor.
    Handles gas production, CO2 sparging, gas-liquid mass transfer, microbial uptake,
    bubble rise/venting, and maintains mass balance. Uses rates and functions from equations module.
    Enforces a hard cap on total headspace pressure.
    """

    def __init__(self, registry, p_max_headspace=None):        
        self.registry = registry
        self.vars = self._load_vars(registry)
        v = self.vars
        # State variables: always initialize with correct units!
        self.C_H2 = 0.0 * ureg('mol/meter**3')
        self.C_CO2 = 0.0 * ureg('mol/meter**3')
        self.n_H2_bubble = 0.0 * ureg('mol')
        self.n_CO2_bubble = 0.0 * ureg('mol')
        self.n_H2_head = 0.0 * ureg('mol')
        self.n_CO2_head = 0.0 * ureg('mol')
        vb_H2 = equations.vb_H2(
            v["rhoL"].magnitude, v["rhoG_H2"].magnitude, v["g"].magnitude,
            v["d_b_H2"].magnitude, v["mu"].magnitude)
        vb_CO2 = equations.vb_CO2(
            v["rhoL"].magnitude, v["rhoG_CO2"].magnitude, v["g"].magnitude,
            v["d_b_CO2"].magnitude, v["mu"].magnitude)
        self.kvent_H2 = equations.kvent_H2(vb_H2, v["Hreactor"].magnitude) * ureg('1/second')
        self.kvent_CO2 = equations.kvent_CO2(vb_CO2, v["Hreactor"].magnitude) * ureg('1/second')
        # Max headspace pressure (default 1.2 atm)
        if p_max_headspace is not None:
            self.p_max_headspace = (p_max_headspace * ureg.atm).to('pascal')
        else:
            self.p_max_headspace = (1.0 * ureg.atm).to('pascal')
    def _load_vars(self, registry):
        required_vars = {
            "V_catholyte": ("V_catholyte", "meter**3"),
            "Vhead": ("Vhead", "meter**3"),
            "kLa_H2": ("kLa_H2", "1/second"),
            "kLa_CO2": ("kLa_CO2", "1/second"),
            "kLa": ("kLa", "1/second"),
            "H_H2": ("H_H2", "mol/meter**3/pascal"),
            "H_CO2": ("H_CO2", "mol/meter**3/pascal"),
            "T": ("T", "kelvin"),
            "R": ("R", "joule/(mol*kelvin)"),
            "rhoL": ("rhoL", "kg/meter**3"),
            "rhoG_H2": ("rhoG_H2", "kg/meter**3"),
            "g": ("g", "meter/second**2"),
            "d_b_H2": ("d_b_H2", "meter"),
            "mu": ("mu", "pascal*second"),
            "rhoG_CO2": ("rhoG_CO2", "kg/meter**3"),
            "d_b_CO2": ("d_b_CO2", "meter"),
            "Hreactor": ("Hreactor", "meter"),
            "V_an": ("V_an", "meter**3"),
            "A_cat": ("A_cat", "meter**2"),
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
            raise ValueError("Missing variables in GasDynamicsModule: " + ", ".join(missing))
        return vars_out

    def step(self, dt, H2_production=0.0, CO2_flow_rate=0.0, CO2_conc=0.0,
             H2_uptake=0.0, CO2_uptake=0.0):

        # === Enforce units on ALL inputs ===
        dt = enforce_quantity(dt, "second")
        H2_production = enforce_quantity(H2_production, "mol/second")
        CO2_flow_rate = enforce_quantity(CO2_flow_rate, "meter**3/second")
        CO2_conc = enforce_quantity(CO2_conc, "mol/meter**3")
        H2_uptake = enforce_quantity(H2_uptake, "mol/second")
        CO2_uptake = enforce_quantity(CO2_uptake, "mol/second")

        v = self.vars

        # === 1. Gas inputs to bubble pool ===
        n_H2_gen = H2_production * dt           # mol
        n_CO2_sparge = equations.nCO2_sparge(CO2_flow_rate, CO2_conc) * dt  # mol
        self.n_H2_bubble += n_H2_gen
        self.n_CO2_bubble += n_CO2_sparge

        # Clamp negatives (shouldn't happen, but just in case)
        self.n_H2_bubble = clamp_nonnegative(self.n_H2_bubble)
        self.n_CO2_bubble = clamp_nonnegative(self.n_CO2_bubble)

        # === 2. Gas-liquid transfer (bubble <-> dissolved) ===
        p_H2_bubble = equations.p_H2(self.n_H2_bubble, v["R"], v["T"], v["Vhead"])
        p_CO2_bubble = equations.pCO2(self.n_CO2_bubble, v["R"], v["T"], v["Vhead"])
        C_H2_sat = equations.C_H2_sat(v["H_H2"], p_H2_bubble)
        C_CO2_sat = equations.C_CO2_sat(v["H_CO2"], p_CO2_bubble)
        J_H2_diss = equations.J_H2_diss(v["kLa_H2"], C_H2_sat, self.C_H2)
        J_CO2_diss = equations.J_CO2_diss(v["kLa_CO2"], C_CO2_sat, self.C_CO2)
        n_H2_diss = equations.n_H2_diss(J_H2_diss, v["V_catholyte"])
        n_CO2_diss = equations.nCO2_diss(J_CO2_diss, v["V_catholyte"])
        n_H2_diss = np.minimum(n_H2_diss.magnitude, self.n_H2_bubble.magnitude) * ureg("mol")
        n_CO2_diss = np.minimum(n_CO2_diss.magnitude, self.n_CO2_bubble.magnitude) * ureg("mol")
        self.n_H2_bubble -= n_H2_diss
        self.n_CO2_bubble -= n_CO2_diss
        self.C_H2 += n_H2_diss / v["V_catholyte"]
        self.C_CO2 += n_CO2_diss / v["V_catholyte"]

        # Clamp again after subtractions
        self.n_H2_bubble = clamp_nonnegative(self.n_H2_bubble)
        self.n_CO2_bubble = clamp_nonnegative(self.n_CO2_bubble)
        self.C_H2 = clamp_nonnegative(self.C_H2)
        self.C_CO2 = clamp_nonnegative(self.C_CO2)

        # === 3. Microbial uptake from dissolved pool ===
        n_H2_uptake = H2_uptake * dt 
        n_CO2_uptake = CO2_uptake * dt
        available_H2 = self.C_H2 * v["V_catholyte"]
        available_CO2 = self.C_CO2 * v["V_catholyte"]
        n_H2_uptake = min(n_H2_uptake, available_H2)
        n_CO2_uptake = min(n_CO2_uptake, available_CO2)
        self.C_H2 -= n_H2_uptake / v["V_catholyte"]
        self.C_CO2 -= n_CO2_uptake / v["V_catholyte"]

        # Clamp negatives again
        self.C_H2 = clamp_nonnegative(self.C_H2)
        self.C_CO2 = clamp_nonnegative(self.C_CO2)

        # === 4. Bubble venting/loss (to headspace or out) ===
        n_H2_vent = equations.n_H2_vent(self.kvent_H2, self.n_H2_bubble) * dt
        n_CO2_vent = equations.n_CO2_vent(self.kvent_CO2, self.n_CO2_bubble) * dt
        n_H2_vent = np.minimum(n_H2_vent.magnitude, self.n_H2_bubble.magnitude) * ureg("mol")
        n_CO2_vent = np.minimum(n_CO2_vent.magnitude, self.n_CO2_bubble.magnitude) * ureg("mol")
        self.n_H2_bubble -= n_H2_vent
        self.n_CO2_bubble -= n_CO2_vent

        self.n_H2_bubble = clamp_nonnegative(self.n_H2_bubble)
        self.n_CO2_bubble = clamp_nonnegative(self.n_CO2_bubble)

        # Add vented to headspace
        self.n_H2_head += n_H2_vent
        self.n_CO2_head += n_CO2_vent

        self.n_H2_head = clamp_nonnegative(self.n_H2_head)
        self.n_CO2_head = clamp_nonnegative(self.n_CO2_head)

        # === 5. Hard cap on headspace pressure ===
        n_total_head_actual = self.n_H2_head + self.n_CO2_head
        n_total_max = (self.p_max_headspace * v["Vhead"]) / (v["R"] * v["T"])
        vented_H2 = 0.0 * ureg('mol')
        vented_CO2 = 0.0 * ureg('mol')

        if n_total_head_actual > n_total_max:
            x_H2 = (self.n_H2_head / n_total_head_actual).magnitude if n_total_head_actual.magnitude > 0 else 0.0
            x_CO2 = (self.n_CO2_head / n_total_head_actual).magnitude if n_total_head_actual.magnitude > 0 else 0.0
            n_H2_head_actual = self.n_H2_head
            n_CO2_head_actual = self.n_CO2_head
            # Set new capped values
            self.n_H2_head = x_H2 * n_total_max
            self.n_CO2_head = x_CO2 * n_total_max
            # Record how much vented away
            vented_H2 = n_H2_head_actual - self.n_H2_head
            vented_CO2 = n_CO2_head_actual - self.n_CO2_head

        # Clamp again
        self.n_H2_head = clamp_nonnegative(self.n_H2_head)
        self.n_CO2_head = clamp_nonnegative(self.n_CO2_head)

        # Recalculate partial pressures
        p_H2_head = equations.p_H2(self.n_H2_head, v["R"], v["T"], v["Vhead"])
        p_CO2_head = equations.pCO2(self.n_CO2_head, v["R"], v["T"], v["Vhead"])

        return {
            "C_H2": self.C_H2,
            "n_H2_bubble": self.n_H2_bubble,
            "n_H2_head": self.n_H2_head,
            "C_CO2": self.C_CO2,
            "n_CO2_bubble": self.n_CO2_bubble,
            "n_CO2_head": self.n_CO2_head,
            "p_H2_bubble": p_H2_bubble,
            "p_CO2_bubble": p_CO2_bubble,
            "p_H2_head": p_H2_head,
            "p_CO2_head": p_CO2_head,
            "n_H2_diss": n_H2_diss,
            "n_CO2_diss": n_CO2_diss,
            "n_H2_vent": n_H2_vent,
            "n_CO2_vent": n_CO2_vent,
            "dissolved_H2_uptake": n_H2_uptake,
            "dissolved_CO2_uptake": n_CO2_uptake,
            "vented_H2": vented_H2,
            "vented_CO2": vented_CO2,
        }

    def get_state(self):
        return {
            "C_H2": self.C_H2,
            "n_H2_bubble": self.n_H2_bubble,
            "n_H2_head": self.n_H2_head,
            "C_CO2": self.C_CO2,
            "n_CO2_bubble": self.n_CO2_bubble,
            "n_CO2_head": self.n_CO2_head,
        }

    def mass_balance_check(self):
        v = self.vars
        total_H2 = self.C_H2 * v["V_catholyte"] + self.n_H2_bubble + self.n_H2_head
        total_CO2 = self.C_CO2 * v["V_catholyte"] + self.n_CO2_bubble + self.n_CO2_head
        return {"total_H2": total_H2, "total_CO2": total_CO2}
