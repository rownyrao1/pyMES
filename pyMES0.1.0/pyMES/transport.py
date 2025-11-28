import numpy as np
from utils.pint_registry import ureg
from scipy.integrate import solve_ivp
import equations

class TransportModule:
    """
    Models proton balances and buffer action in an MES reactor with anode and cathode compartments.
    """

    def __init__(self, registry, ec_results=None):
        self.registry = registry
        self.ec_results = ec_results or {}
        self.vars = self._load_vars()

        # Initial anolyte/catholyte H+ from pH and volume (in mol, Pint Quantities)
        self.V_an = self.vars["V_an"].to("meter**3")
        self.V_cat = self.vars["V_catholyte"].to("meter**3")

        self.N_H_an = self.initial_proton_moles(self.vars["pH_an"], self.V_an)
        self.N_H_cat = self.initial_proton_moles(self.vars["pH_cat"], self.V_cat)

        # Buffer: convert acid/base from g/L to mol/L, then total moles = conc * volume (L)
        C_acid = self.vars["C_buffer_acid"].to("g/liter")
        MW_acid = self.vars["MW_buffer_acid"].to("g/mol")
        C_base = self.vars["C_buffer_base"].to("g/liter")
        MW_base = self.vars["MW_buffer_base"].to("g/mol")
        V_cat_L = self.V_cat.to("liter")

        self.buffer_pKa = float(self.vars["buffer_pKa"])

        acid_mol_L = self.buffer_mass_to_molarity(C_acid, MW_acid)
        base_mol_L = self.buffer_mass_to_molarity(C_base, MW_base)

        self.N_buffer_acid = (acid_mol_L * V_cat_L).to("mol")
        self.N_buffer_base = (base_mol_L * V_cat_L).to("mol")

    def _load_vars(self):
        required_vars = {
            "phi_mig": ("phi_mig", ""),                  # dimensionless
            "F": ("F", "coulomb / mole"),
            "A_memb": ("A_memb", "meter**2"),
            "membrane_thickness": ("membrane_thickness", "meter"),
            "D_Hplus": ("D_Hplus", "meter**2/second"),
            "pH_an": ("pH_an", ""),
            "pH_cat": ("pH_cat", ""),
            "flow_rate": ("flow_rate", "meter**3/second"),
            "CO2_conc": ("CO2_conc", "mol/meter**3"),
            "R": ("R", "joule/(mol*kelvin)"),
            "T": ("T", "kelvin"),
            "Vhead": ("Vhead", "meter**3"),
            "H_H2": ("H_H2", "mol/(meter**3*pascal)"),
            "H_CO2": ("H_CO2", "mol/(meter**3*pascal)"),
            "kLa": ("kLa", "1/second"),
            "V_catholyte": ("V_catholyte", "meter**3"),
            "rhoL": ("rhoL", "kg/meter**3"),
            "rhoG_H2": ("rhoG_H2", "kg/meter**3"),
            "g": ("g", "meter/second**2"),
            "d_b_H2": ("d_b_H2", "meter"),
            "mu": ("mu", "pascal*second"),
            "rhoG_CO2": ("rhoG_CO2", "kg/meter**3"),
            "d_b_CO2": ("d_b_CO2", "meter"),
            "Hreactor": ("Hreactor", "meter"),
            "V_an": ("V_an", "meter**3"),
            "C_buffer_acid": ("C_buffer_acid", "g/liter"),
            "C_buffer_base": ("C_buffer_base", "g/liter"),
            "MW_buffer_acid": ("MW_buffer_acid", "g/mol"),
            "MW_buffer_base": ("MW_buffer_base", "g/mol"),
            "buffer_pKa": ("buffer_pKa", ""),             # dimensionless
        }
        vars_out = {}
        missing = []
        for varname, (regkey, req_unit) in required_vars.items():
            try:
                value = self.registry.get(regkey)
                if req_unit and not hasattr(value, "to"):
                    value = ureg.Quantity(value, req_unit)
                elif req_unit:
                    value = value.to(req_unit)
                vars_out[varname] = value
            except Exception as e:
                missing.append(f"{regkey} ({str(e)})")
        if missing:
            raise ValueError("Missing variables in TransportModule: " + ", ".join(missing))
        if self.ec_results:
            for k, v in self.ec_results.items():
                vars_out[k] = v
        vars_out.update(self.ec_results)
        return vars_out

    @staticmethod
    def buffer_mass_to_molarity(g_per_L, MW_g_per_mol):
        """Convert buffer concentration from g/L to mol/L (all Pint Quantities)."""
        return (g_per_L / MW_g_per_mol).to("mol/liter")

    @staticmethod
    def initial_proton_moles(pH, volume_m3):
        """
        Calculate total moles of free H+ from pH and volume (Pint Quantity, m^3).
        Returns: Pint Quantity (mol)
        """
        Hplus_conc = 10 ** (-pH) * ureg("mol/liter")
        volume_L = volume_m3.to("liter")
        return (Hplus_conc * volume_L).to("mol")

    @staticmethod
    def buffer_pH(base_mol_L, acid_mol_L, pKa):
        """Calculate buffer pH using Henderson–Hasselbalch equation."""
        if acid_mol_L <= 0 or base_mol_L <= 0:
            raise ValueError("Both acid and base concentrations must be > 0.")
        ratio = base_mol_L / acid_mol_L
        return pKa + np.log10(ratio)

    @staticmethod
    def buffer_update_after_proton_addition(
        N_H_cat, N_base, N_acid, V_cat_L, pKa, protons_added
    ):
        N_H_cat_new = N_H_cat + protons_added
        n_buffered = min(N_H_cat_new.magnitude, N_base.magnitude) * ureg.mol if N_H_cat_new > 0 and N_base > 0 else 0 * ureg.mol
        N_base_new = N_base - n_buffered
        N_acid_new = N_acid + n_buffered
        N_H_cat_new = N_H_cat_new - n_buffered
        # pH calculation
        if N_base_new > 0 * ureg.mol and N_acid_new > 0 * ureg.mol:
            base_conc = (N_base_new / V_cat_L).to("mol/liter")
            acid_conc = (N_acid_new / V_cat_L).to("mol/liter")
            pH = pKa + np.log10(base_conc.magnitude / acid_conc.magnitude)
        else:
            Hplus_conc = (N_H_cat_new / V_cat_L).to("mol/liter")
            pH = -np.log10(Hplus_conc.magnitude) if Hplus_conc.magnitude > 0 else 14
        return pH, N_base_new, N_acid_new, N_H_cat_new

    def anode_proton_production(self, F=None):
        if F is None:
            F = self.vars["F"].to("coulomb/mole")
        I = -self.ec_results.get("current")
        if I is None:
            raise ValueError("Current not found in ec_results")
        if not hasattr(I, "to"):
            I = ureg.Quantity(I, "A")
        return (I / F).to("mol/second")

    def proton_transfer_flux(self):
        v = self.vars
        I = -self.ec_results.get("current")
        if not hasattr(I, "to"):
            I = ureg.Quantity(I, "A")
        F = v["F"].to("coulomb/mole")
        phi_mig = v["phi_mig"]
        A_memb = v["A_memb"].to("meter**2")
        D_Hplus = v["D_Hplus"].to("meter**2/second")
        membrane_thickness = v["membrane_thickness"].to("meter")
        # All Pint Quantities below:
        H_an = (self.N_H_an / self.V_an).to("mol/meter**3")
        H_cat = (self.N_H_cat / self.V_cat).to("mol/meter**3")
        r_H_gen = self.anode_proton_production(F)
        J_mig = (phi_mig * I / (F * A_memb)).to("mol/meter**2/second")
        J_diff = (D_Hplus / membrane_thickness * (H_an - H_cat)).to("mol/meter**2/second")
        J_total = J_mig + J_diff
        r_H_trans = (J_total * A_memb).to("mol/second")
        return r_H_gen, J_mig, J_diff, J_total, r_H_trans

    def update_anolyte_state(self, r_H_gen, r_H_trans, dt):
        if not hasattr(dt, "units"):
            dt = ureg.Quantity(dt, "second")
        # Protons produced at anode (water oxidation)
        protons_produced = r_H_gen * dt

        # Water consumption: 2 H2O → 4 H+ + 4 e- + O2, so 1 mol H+ produced = 0.5 mol H2O consumed
        water_consumed = 0.5 * protons_produced

        # Water molar volume = 18 mL/mol = 18e-6 m³/mol
        dV_an = -water_consumed * 18e-6 * ureg("meter**3 / mole")

        # Update anolyte volume (subtract water lost)
        self.V_an = self.V_an + dV_an

        # Prevent negative/zero volume
        if self.V_an.magnitude < 1e-8:
            self.V_an = 1e-8 * ureg("meter**3")

        N_H_an_new = self.N_H_an + (r_H_gen - r_H_trans) * dt
        Hplus_an_conc = (N_H_an_new / self.V_an).to("mol/meter**3")
        pH_an = -np.log10(Hplus_an_conc.magnitude / 1000) if Hplus_an_conc.magnitude > 0 else 14
        self.N_H_an = N_H_an_new
        return N_H_an_new, Hplus_an_conc, pH_an

    def update_catholyte_state(self, r_H_trans, dt, H2_proton_used=None):
        if not hasattr(dt, "units"):
            dt = ureg.Quantity(dt, "second")
        protons_transferred = r_H_trans * dt
        if H2_proton_used is None:
            H2_proton_used = protons_transferred  # Default: all protons are used in HER
        elif not hasattr(H2_proton_used, "units"):
            H2_proton_used = ureg.Quantity(H2_proton_used, "mol")
        protons_excess = protons_transferred - H2_proton_used
        V_cat_L = self.V_cat.to("liter")
        if protons_excess.magnitude > 0:
            N_H_cat_new = self.N_H_cat + protons_excess
            pH_cat, N_buffer_base_new, N_buffer_acid_new, N_H_cat_final = self.buffer_update_after_proton_addition(
                N_H_cat_new, self.N_buffer_base, self.N_buffer_acid, V_cat_L, self.buffer_pKa, protons_excess
            )
            self.N_H_cat = N_H_cat_final
            self.N_buffer_base = N_buffer_base_new
            self.N_buffer_acid = N_buffer_acid_new
        else:
            N_H_cat_new = self.N_H_cat  # unchanged
            Hplus_conc = (N_H_cat_new / self.V_cat).to("mol/liter")
            pH_cat = -np.log10(Hplus_conc.magnitude) if Hplus_conc.magnitude > 0 else 14
        return pH_cat, self.N_H_cat, self.N_buffer_base, self.N_buffer_acid

    def step(self, dt, H2_proton_used=None):
        if not hasattr(dt, "units"):
            dt = ureg.Quantity(dt, "second")
        r_H_gen, J_mig, J_diff, J_total, r_H_trans = self.proton_transfer_flux()
        N_H_an_new, Hplus_an_conc, pH_an = self.update_anolyte_state(r_H_gen, r_H_trans, dt)

        # Proton-limited hydrogen production
        n_H2_proton_limited = (r_H_trans / 2).to("mol/second")

        # Final hydrogen production rate (mol/s), limited by whichever is less
        H2_production_rate = n_H2_proton_limited

        # Default: all protons used for HER is the amount actually consumed (for buffer accounting)
        if H2_proton_used is None:
            H2_proton_used = 2 * H2_production_rate * dt  # moles of protons consumed at the cathode

        pH_cat, N_H_cat, N_buffer_base, N_buffer_acid = self.update_catholyte_state(
            r_H_trans, dt, H2_proton_used=H2_proton_used
        )

        return {
            "N_H_an": N_H_an_new,
            "Hplus_an": Hplus_an_conc,
            "pH_an": pH_an,
            "N_H_cat": N_H_cat,
            "Hplus_cat": (N_H_cat / self.V_cat).to("mol/meter**3"),
            "pH_cat": pH_cat,
            "N_buffer_base": N_buffer_base,
            "N_buffer_acid": N_buffer_acid,
            "J_mig": J_mig,
            "J_diff": J_diff,
            "J_total": J_total,
            "r_H_gen": r_H_gen,
            "r_H_trans": r_H_trans,
            "H2_production": H2_production_rate,
            "V_an": self.V_an,
            "n_H2_proton_limited": n_H2_proton_limited,
        }
