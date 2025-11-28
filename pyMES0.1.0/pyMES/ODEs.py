import numpy as np
from scipy.integrate import solve_ivp
from utils.pint_registry import ureg

class AnodeODESolver:
    """
    Class for solving ODEs describing accumulation of electrons, protons, and O2 at the anode.
    """

    def __init__(self, F):
        """
        F: Faraday constant (Pint Quantity, e.g., 96485 * ureg.coulomb / ureg.mole)
        """
        self.F = F.to("coulomb/mole").magnitude

    def odes(self, t, y, Ival):
        """
        ODE system: y = [Ne, NH, NO2]
        Ival: Current in Ampere (float, abs value)
        """
        dNe_dt = Ival / self.F          # Electron production rate [mol/s]
        dNH_dt = Ival / self.F          # Proton production rate [mol/s]
        dNO2_dt = Ival / (4 * self.F)   # Oxygen production rate [mol/s]
        return [dNe_dt, dNH_dt, dNO2_dt]

    def integrate(self, I, t_span=(0, 3600), N0=None, num_points=100):
        """
        Integrate ODEs over the time span t_span.

        Parameters:
            I: Current (Pint Quantity, direction ignored; abs used)
            t_span: (start, stop) in seconds, default 0 to 3600 s
            N0: Initial [Ne, NH, NO2] in mol (default zeros)
            num_points: Time resolution
        Returns:
            dict: time, Ne, NH, NO2 (all Pint Quantities)
        """
        Ival = abs(I.to("A").magnitude)
        if N0 is None:
            N0 = [0.0, 0.0, 0.0]

        t_eval = np.linspace(t_span[0], t_span[1], num_points)
        sol = solve_ivp(lambda t, y: self.odes(t, y, Ival), t_span, N0, t_eval=t_eval, method="RK45")
        Ne = sol.y[0] * ureg("mol")
        NH = sol.y[1] * ureg("mol")
        NO2 = sol.y[2] * ureg("mol")
        time = sol.t * ureg("second")
        return {"time": time, "Ne": Ne, "NH": NH, "NO2": NO2}

class GasPoolODESolver:
    """
    Handles the ODEs for H2 and CO2 pools (dissolved and bubble).
    """
    def __init__(self, params):
        self.p = params  # All params (floats or numpy scalars, not pint objects)

    def rhs(self, t, y):
        # Unpack variables and parameters
        C_H2, n_H2_bubble, C_CO2, n_CO2_bubble = y
        p = self.p

        # All params must be in SI units, as floats (no Pint inside here)
        p_H2 = n_H2_bubble * p["R"] * p["T"] / p["Vhead"]
        p_CO2 = n_CO2_bubble * p["R"] * p["T"] / p["Vhead"]
        C_H2_sat = p["H_H2"] * p_H2
        C_CO2_sat = p["H_CO2"] * p_CO2

        J_H2_diss = p["kLa"] * (C_H2_sat - C_H2)
        J_CO2_diss = p["kLa"] * (C_CO2_sat - C_CO2)
        n_H2_diss = J_H2_diss * p["V_liq"]
        n_CO2_diss = J_CO2_diss * p["V_liq"]

        n_H2_vent = p["kvent_H2"] * n_H2_bubble
        n_CO2_vent = p["kvent_CO2"] * n_CO2_bubble

        # Microbial uptake and source rates should be float-returning callables or floats
        n_H2_prod = p.get("n_H2_prod_func", lambda t: p.get("n_H2_prod", 0.0))(t)
        n_CO2_prod = p.get("n_CO2_prod_func", lambda t: p.get("n_CO2_prod", 0.0))(t)
        n_H2_sparge = p.get("n_H2_sparge_func", lambda t: p.get("n_H2_sparge", 0.0))(t)
        n_CO2_sparge = p.get("n_CO2_sparge_func", lambda t: p.get("n_CO2_sparge", 0.0))(t)
        n_H2_uptake = p.get("n_H2_uptake_func", lambda C, t: p.get("n_H2_uptake", 0.0))(C_H2, t)
        n_CO2_uptake = p.get("n_CO2_uptake_func", lambda C, t: p.get("n_CO2_uptake", 0.0))(C_CO2, t)

        dC_H2_dt = (n_H2_diss - n_H2_uptake) / p["V_liq"]
        dn_H2_bubble_dt = n_H2_sparge + n_H2_prod - n_H2_diss - n_H2_vent

        dC_CO2_dt = (n_CO2_diss - n_CO2_uptake) / p["V_liq"]
        dn_CO2_bubble_dt = n_CO2_sparge + n_CO2_prod - n_CO2_diss - n_CO2_vent

        return [dC_H2_dt, dn_H2_bubble_dt, dC_CO2_dt, dn_CO2_bubble_dt]

    def integrate(self, y0, t_span=(0, 3600), num_points=1000):
        t_eval = np.linspace(t_span[0], t_span[1], num_points)
        sol = solve_ivp(self.rhs, t_span, y0, t_eval=t_eval, method="RK45")
        return sol
