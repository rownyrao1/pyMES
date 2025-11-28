import numpy as np
from scipy.optimize import root_scalar    
from utils.pint_registry import ureg
from ODEs import AnodeODESolver
import equations

class ElectrochemistryModule:
    def __init__(self, registry, mode, applied_value, logger=None):
        self.registry = registry
        self.mode = mode.strip().lower()
        self.applied_value = applied_value
        self.logger = logger
        self.vars = self._load_vars()

    def _load_vars(self):
        """
        Loads all required variables from the registry,
        converts to correct units, and returns a dict.
        Errors if any required variable is missing or cannot be converted.
        """
        required_vars = {
            "pH_cat": ("pH_cat", ""),                 # dimensionless
            "pH_an": ("pH_an", ""),
            "d_an_memb": ("d_an_memb", "m"),
            "A_an_cs": ("A_an_cs", "m^2"),
            "A_an": ("A_an", "m^2"),
            "d_cat_memb": ("d_cat_memb", "m"),
            "A_cat_cs": ("A_cat_cs", "m^2"),
            "A_cat": ("A_cat", "m^2"),
            "z_list_cat": ("z_list_cat", ""),
            "z_list_an": ("z_list_an", ""),
            "D_list_cat": ("D_list_cat", "m^2 / s"),
            "D_list_an": ("D_list_an", "m^2 / s"),
            "c_list_cat": ("c_list_cat", "mol / m^3"),
            "c_list_an": ("c_list_an", "mol / m^3"),
            "An_resistivity": ("An_resistivity", "ohm * m"),
            "l_an": ("l_an", "m"),
            "Cat_resistivity": ("Cat_resistivity", "ohm * m"),
            "l_cat": ("l_cat", "m"),
            "membrane_thickness": ("membrane_thickness", "m"),
            "k_memb": ("k_memb", "siemens / meter"),
            "A_memb": ("A_memb", "m^2"),
            "F": ("F", "coulomb / mole"),             # Faraday's constant
            "R": ("R", "joule / (mole * kelvin)"),    # Universal gas constant
            "T": ("T", "kelvin"),
            "j0_cat": ("j0_cat", "A/m^2"),
            "alpha_cat": ("alpha_cat", ""),           # dimensionless
            "n_cat": ("n_cat", ""),                   # dimensionless
            "j0_an": ("j0_an", "A/m^2"),
            "alpha_an": ("alpha_an", ""),             # dimensionless
            "n_an": ("n_an", ""),                     # dimensionless
        }

        vars_out = {}
        missing_vars = []
        errors = []

        for varname, (regkey, req_unit) in required_vars.items():
            try:
                value = self.registry.get(regkey)
            except Exception:
                missing_vars.append(regkey)
                continue

            # Attach or convert units if specified
            if req_unit and req_unit.strip():
                if not hasattr(value, "to"):
                    try:
                        value = ureg.Quantity(value, req_unit)
                    except Exception as e:
                        errors.append(f"{regkey} (could not attach unit {req_unit}: {e})")
                        continue
                else:
                    try:
                        value = value.to(req_unit)
                    except Exception as e:
                        errors.append(f"{regkey} (cannot convert to {req_unit}: {e})")
                        continue

            # Check for NaN
            if hasattr(value, "magnitude") and isinstance(value.magnitude, (float, int)):
                if np.isnan(value.magnitude):
                    errors.append(f"{regkey} is NaN!")
            elif hasattr(value, "magnitude") and hasattr(value.magnitude, "__iter__"):
                if np.any(np.isnan(value.magnitude)):
                    errors.append(f"{regkey} contains NaN values!")

            vars_out[varname] = value

        if missing_vars:
            raise ValueError(
                f"Missing required variables from registry: {', '.join(missing_vars)}"
            )
        if errors:
            raise ValueError(
                "The following errors were found loading variables:\n" + "\n".join(errors)
            )

        return vars_out

    def run(self):
        mode_map = {
            "potentiostatic": "potentiostatic",
            "pot": "potentiostatic",
            "ca": "potentiostatic",
            "galvanostatic": "galvanostatic",
            "gal": "galvanostatic",
            "cp": "galvanostatic"
        }
        mode_key = self.mode.strip().lower()
        if mode_key not in mode_map:
            raise ValueError(f"Invalid mode '{self.mode}'. Allowed: {list(mode_map.keys())}")
        selected_mode = mode_map[mode_key]

        # Enforce correct units for the applied_value
        if selected_mode == "potentiostatic":
            # should be volts
            if not hasattr(self.applied_value, "units") or not self.applied_value.check('[volt]'):
                self.applied_value = self.applied_value * ureg.volt
            return self._run_potentiostatic()
        elif selected_mode == "galvanostatic":
            # should be amps
            if not hasattr(self.applied_value, "units") or not self.applied_value.check('[current]'):
                self.applied_value = self.applied_value * ureg.ampere
            return self._run_galvanostatic()
        else:
            raise ValueError("Unexpected internal error: mode switch failed.")
      
    def _precompute(self, v):
        """Precompute all resistances, conductivities, and equilibrium potentials."""
        # 1. Equilibrium potentials
        E_eq_cat = equations.E_eq_cat(v["pH_cat"]) * ureg("V")
        E_eq_an = equations.E_eq_an(v["pH_an"]) * ureg("V")

        # 2. Ionic conductivity (use lists)
        k_an = equations.k_anolyte(
            v["z_list_an"], v["F"], v["D_list_an"], v["R"], v["T"], v["c_list_an"]
        )
        k_cat = equations.k_catholyte(
            v["z_list_cat"], v["F"], v["D_list_cat"], v["R"], v["T"], v["c_list_cat"]
        )

        # 3. Resistances
        R_anolyte = equations.R_anolyte(v["d_an_memb"], k_an, v["A_an"])
        R_catholyte = equations.R_catholyte(v["d_cat_memb"], k_cat, v["A_cat"])
        R_anode = equations.R_anode(v["An_resistivity"], v["l_an"], v["A_an_cs"])
        R_cathode = equations.R_cathode(v["Cat_resistivity"], v["l_cat"], v["A_cat_cs"])
        R_membrane = equations.R_membrane(v["membrane_thickness"], v["k_memb"], v["A_memb"])
        R_local_cat = equations.R_local_cat(R_cathode, R_catholyte, R_membrane)
        R_local_an = equations.R_local_an(R_anode, R_anolyte)
        R_total = equations.R_total(R_local_cat, R_local_an)

        resistances = {
            "R_anolyte": R_anolyte,
            "R_catholyte": R_catholyte,
            "R_anode": R_anode,
            "R_cathode": R_cathode,
            "R_membrane": R_membrane,
            "R_local_cat": R_local_cat,
            "R_local_an": R_local_an,
            "R_total": R_total,
        }

        conductivities = {
            "k_anolyte": k_an,
            "k_catholyte": k_cat,
        }

        eq_potentials = {
            "E_eq_cat": E_eq_cat,
            "E_eq_an": E_eq_an,
        }

        return eq_potentials, conductivities, resistances
        
    def _run_potentiostatic(self):
        v = self.vars
        ureg = v["A_cat"].units._REGISTRY
        eq_potentials, conductivities, resistances = self._precompute(v)

        E_app = self.applied_value.to("volt")
        E_eq_cat = eq_potentials["E_eq_cat"]
        E_eq_an = eq_potentials["E_eq_an"]
        R_local_cat_ = resistances["R_local_cat"]
        R_local_an_ = resistances["R_local_an"]
        

        A_cat = v["A_cat"]
        A_an = v["A_an"]
        j0_cat = v["j0_cat"]
        j0_an = v["j0_an"]
        alpha_cat = v["alpha_cat"]
        alpha_an = v["alpha_an"]
        n_cat = v["n_cat"]
        n_an = v["n_an"]
        F = v["F"]
        R = v["R"]
        T = v["T"]

        def current_eq(I_val):
            I = I_val
            deltaV = E_app.magnitude - E_eq_cat.magnitude - I * R_local_cat_.magnitude
            exponent = (-alpha_cat * n_cat * F * deltaV / (R * T)).to_base_units().magnitude
            j = (I / A_cat).magnitude
            j0 = j0_cat.magnitude
            return (j + j0 * np.exp(exponent))
        sol = root_scalar(current_eq, bracket=[-1, 1], method='brentq')
        if not sol.converged:
            raise RuntimeError("Root finding failed to converge.")
        I_solution = sol.root * ureg.ampere


        anode_production = self._compute_anode_production(I=I_solution)

        F = self.vars["F"]
        anode_ode_solver = AnodeODESolver(F)
        ode_profile = anode_ode_solver.integrate(I=-I_solution, t_span=(0, 3600), num_points=100)


        j_cat, j_an = self._compute_current_densities(-I_solution, A_cat, A_an)
        eta_an, E_an, V_cell = self._compute_anode_cell_params(
            I=-I_solution, E_cat=E_app, E_eq_an=E_eq_an, j_an=j_an, j0_an=j0_an,
            alpha_an=alpha_an, n_an=n_an, F=F, R=R, T=T, R_local_an=R_local_an_
        )

        return {
            "mode": "potentiostatic",
            "current": I_solution,
            "E_app": E_app,
            "E_eq_cat": E_eq_cat,
            "R_local_cat": R_local_cat_,
            "j_cat": j_cat,
            "j_an": j_an,
            "eta_an": eta_an,
            "E_an": E_an,
            "V_cell": V_cell,
            "conductivities": conductivities,
            "eq_potentials": eq_potentials,
            "resistances": resistances,
            "anode_production": anode_production,
            "anode_ode_profile": ode_profile,
        }
    

    def _run_galvanostatic(self):
        v = self.vars
        ureg = v["A_cat"].units._REGISTRY
        eq_potentials, conductivities, resistances = self._precompute(v)

        E_eq_cat = eq_potentials["E_eq_cat"]
        E_eq_an = eq_potentials["E_eq_an"]
        R_local_cat_ = resistances["R_local_cat"]
        R_local_an_ = resistances["R_local_an"]

        I_app = self.applied_value * ureg.ampere if not hasattr(self.applied_value, "units") else self.applied_value.to("A")
        A_cat = v["A_cat"]
        A_an = v["A_an"]
        j0_cat = v["j0_cat"]
        j0_an = v["j0_an"]
        alpha_cat = v["alpha_cat"]
        alpha_an = v["alpha_an"]
        n_cat = v["n_cat"]
        n_an = v["n_an"]
        F = v["F"]
        R = v["R"]
        T = v["T"]

        def voltage_eq(E_cat_val):
            E_cat = E_cat_val * ureg.volt
            deltaV = E_cat.magnitude - E_eq_cat.magnitude - I_app.magnitude * R_local_cat_.magnitude
            exponent = (-alpha_cat * n_cat * F * deltaV / (R * T)).to_base_units().magnitude
            j = (I_app / A_cat).to("A/m^2").magnitude
            j0 = j0_cat.to("A/m^2").magnitude
            return j + j0 * np.exp(exponent)
        sol = root_scalar(voltage_eq, bracket=[-2, 2], method='brentq')
        if not sol.converged:
            raise RuntimeError("Root finding for E_app failed to converge.")
        E_cat_solution = sol.root * ureg.volt

        anode_production = self._compute_anode_production(I=-I_app)

        F = self.vars["F"]
        anode_ode_solver = AnodeODESolver(F)
        ode_profile = anode_ode_solver.integrate(I=-I_app, t_span=(0, 3600), num_points=100)


        j_cat, j_an = self._compute_current_densities(-I_app, A_cat, A_an)
        eta_an, E_an, V_cell = self._compute_anode_cell_params(
            I=-I_app, E_cat=E_cat_solution, E_eq_an=E_eq_an, j_an=j_an, j0_an=j0_an,
            alpha_an=alpha_an, n_an=n_an, F=F, R=R, T=T, R_local_an=R_local_an_
        )

        return {
            "mode": "galvanostatic",
            "I_app": I_app,
            "E_cat": E_cat_solution,
            "E_eq_cat": E_eq_cat,
            "R_local_cat": R_local_cat_,
            "j_cat": j_cat,
            "j_an": j_an,
            "eta_an": eta_an,
            "E_an": E_an,
            "V_cell": V_cell,
            "conductivities": conductivities,
            "eq_potentials": eq_potentials,
            "resistances": resistances,
            "anode_production": anode_production,
            "anode_ode_profile": ode_profile,
        }

    
    def _compute_current_densities(self, I, A_cat, A_an):
        j_cat = (I / A_cat).to("A/m^2")
        j_an = (I / A_an).to("A/m^2")
        return j_cat, j_an

    def _compute_anode_cell_params(self, I, E_cat, E_eq_an, j_an, j0_an, alpha_an, n_an, F, R, T, R_local_an):

        E_eq_an = E_eq_an

        # Tafel equation
        ratio = (j_an / j0_an).to_base_units().magnitude
        if ratio <= 0 or np.isnan(ratio):
            raise ValueError(f"Invalid argument for log: j_an/j0_an = {ratio}")
        eta_an = (R * T / (alpha_an * n_an * F)) * np.log(ratio)
        eta_an = eta_an.to('volt')
        
        #E_an = E_eq_an + eta_an + (I * R_local_an)
        E_an = equations.E_an(E_eq_an, eta_an, I, R_local_an)
        E_an = E_an.to('volt')

        V_cell = equations.V_cell (E_an, E_cat)
        V_cell = V_cell.to('volt')

        return eta_an, E_an, V_cell

    def _compute_anode_production(self, I):
        """
        Calculates rates of electrons, protons, and oxygen produced at the anode (mol/s).
        """
        F = self.vars["F"]
        I = abs(I.to("A"))

        n_e_rate = (I / F).to("mol/s")
        n_H_rate = n_e_rate
        n_O2_rate = (n_e_rate / 4).to("mol/s")

        return {
            "anode_electron_prod": n_e_rate,
            "anode_proton_prod": n_H_rate,
            "anode_oxygen_prod": n_O2_rate,
        }
