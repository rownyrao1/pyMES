import numpy as np
from scipy.optimize import root_scalar
from pint import Quantity
from utils.pint_registry import ureg


# 1. Equilibrium cathode potential (HER)
def E_eq_cat(pH_cat):
    return -0.059 * pH_cat  # [V vs SHE]

def E_eq_an(pH_an):
    return 1.23 + 0.059 * pH_an  # [V vs SHE]

# 3. Anolyte resistance
def R_anolyte(d_an_memb, k_anolyte, A_an):
    return d_an_memb / (k_anolyte * A_an)  # [Ω]

# 4. Catholyte resistance
def R_catholyte(d_cat_memb, k_catholyte, A_cat):
    return d_cat_memb / (k_catholyte * A_cat)  # [Ω]

# 5. Ionic conductivity (general)
def kappa_total(z_list, F, D_list, R, T, c_list):
    """
    Computes ionic conductivity (S/m) given all lists as Pint quantities with correct units.

    Args:
        z_list (list): List of ion charges (dimensionless)
        F (pint.Quantity): Faraday constant (coulomb/mole)
        D_list (list): List of diffusion coefficients (meter^2/second)
        R (pint.Quantity): Universal gas constant (joule/kelvin/mole)
        T (pint.Quantity): Temperature (kelvin)
        c_list (list): List of concentrations (mole/meter^3)

    Returns:
        pint.Quantity: Ionic conductivity (siemens/meter)
    """

    assert all(hasattr(Di, 'units') for Di in D_list), "D_list elements must have Pint units"
    assert all(hasattr(ci, 'units') for ci in c_list), "c_list elements must have Pint units"

    kappa = 0 * F**2 * D_list[0] * c_list[0] / (R * T)  # Initialize with correct units

    for zi, Di, ci in zip(z_list, D_list, c_list):
        term = (zi ** 2) * F ** 2 * Di * ci / (R * T)
        kappa += term
    kappa = kappa.to_base_units()

    try:
        kappa_S_per_m = kappa.to("S/m")
    except Exception as e:
        print("Conversion error:", e)
        kappa_S_per_m = kappa  # fallback to base units

    return kappa_S_per_m
    

def k_anolyte(z_list_an, F, D_list_an, R, T, c_list_an):
    return kappa_total(z_list_an, F, D_list_an, R, T, c_list_an)

def k_catholyte(z_list_cat, F, D_list_cat, R, T, c_list_cat):
    return kappa_total(z_list_cat, F, D_list_cat, R, T, c_list_cat)

# 6. Anode resistance
def R_anode(An_resistivity, l_an, A_an):
    return An_resistivity * l_an / A_an  # [Ω]

# 7. Cathode resistance
def R_cathode(Cat_resistivity, l_cat, A_cat):
    return Cat_resistivity * l_cat / A_cat  # [Ω]

# 8. Membrane resistance
def R_membrane(membrane_thickness, k_memb, A_memb):
    return membrane_thickness / (k_memb * A_memb)  # [Ω]

# 9. Local resistance at cathodic chamber
def R_local_cat(R_cathode, R_catholyte, R_membrane):
    return R_cathode + R_catholyte + R_membrane  # [Ω]

# 10. Local resistance at anodic chamber
def R_local_an(R_anode, R_anolyte):
    return R_anode + R_anolyte  # [Ω]

# 11. Total resistance of the system
def R_total(R_local_cat, R_local_an):
    return R_local_cat + R_local_an  # [Ω]

# 12. Cathodic current density (galvanostatic)
def j_cat(I, A_cat):
    return I / A_cat  # [A/cm²]

# 13. Cathodic overpotential (potentiostatic)
def eta_cat_pot(E_app, E_eq_cat, I, R_local_cat):
    return E_app - E_eq_cat - I * R_local_cat  # [V vs SHE]

# 14. Cathodic overpotential (galvanostatic)
def eta_cat_gal(I_app, j0_cat, alpha_cat, n_cat, F, R, T, A_cat):
    j = I_app / A_cat
    return -R * T / (alpha_cat * n_cat * F) * np.log(-j / j0_cat)  # [V vs SHE]

def I_pot(j0_cat, alpha_cat, n_cat, F, E_app, E_eq_cat, R_local_cat, R, T, A_cat):
    j0 = j0_cat
    alpha = float(alpha_cat)
    n = float(n_cat)
    F_ = F
    Eapp = E_app
    Eeq = E_eq_cat
    Rloc = R_local_cat
    R_ = R
    T_ = T
    A = A_cat

    def equation(I_pot):
        exp_term = -alpha * n * F_ * (Eapp - Eeq - I_pot * Rloc) / (R_ * T_)
        return (I_pot / A) + j0 * np.exp(exp_term)

    # Solve for I_pot numerically (bracket may need adjustment)
    sol = root_scalar(equation, bracket=[-100, 100], method='brentq')
    if not sol.converged:
        raise RuntimeError("Failed to solve for I_pot")
    
# 16. Calculated cathode potential (galvanostatic)
def E_cat(E_app, eta_cat_gal):
    return E_app - eta_cat_gal  # [V vs SHE]

# 17. Anodic current density
def j_an(I, A_an):
    return I / A_an  # [A/cm²]

# 18. Anodic overpotential
def eta_an(j_an, j0_an, alpha_an, n_an, F, R, T):
    return R * T / (alpha_an * n_an * F) * np.log(j_an / j0_an)  # [V vs SHE]

# 19. Calculated anode potential
def E_an(E_eq_an, eta_an, I, R_local_an):
    return E_eq_an + eta_an + I * R_local_an  # [V vs SHE]

# 20. Cell potential
def V_cell(E_an, E_cat):
    return E_an - E_cat  # [V vs SHE]

# 21. Proton flux (membrane)
def J_Hplus(phi_mig, I, F, A_memb, membrane_thickness, D_Hplus, Hplus_an, Hplus_cat):
    migration = phi_mig * I / (F * A_memb)
    diffusion = D_Hplus * (Hplus_an - Hplus_cat)/membrane_thickness
    return migration + diffusion  # [mol m⁻² s⁻¹]

# 22. Proton supply rate at the cathode
def n_Hplus_in(J_Hplus, A_memb):
    return J_Hplus * A_memb  # [mol s⁻¹]

# 23. Max possible H₂ formation (from electrons)
def n_H2_max(n_electrons, n_Hplus_in):
    return 0.5 * min (n_electrons, n_Hplus_in) # [mol s⁻¹]

# 24. CO₂ input (sparging)
def nCO2_sparge(flow_rate, CO2_conc):
    return flow_rate * CO2_conc  # [mol s⁻¹]

# 25. Partial pressure of H₂ (bubble pool, gas law)
def p_H2(n_H2_bubble, R, T, Vhead):
    return n_H2_bubble * R * T / Vhead  # [Pa]

# 26. Partial pressure of CO₂ (bubble pool)
def pCO2(nCO2_bubble, R, T, Vhead):
    return nCO2_bubble * R * T / Vhead  # [Pa]

# 27. H₂ saturation concentration
def C_H2_sat(H_H2, p_H2):
    return H_H2 * p_H2  # [mol m⁻³]

# 28. CO₂ saturation concentration
def C_CO2_sat(H_CO2, pCO2):
    return H_CO2 * pCO2  # [mol m⁻³]

# 29. H₂ gas dissolution flux
def J_H2_diss(kLa, C_H2_sat, C_H2):
    return kLa * (C_H2_sat - C_H2)  # [mol m⁻³ s⁻¹]

# 30. CO₂ gas dissolution flux
def J_CO2_diss(kLa, C_CO2_sat, C_CO2):
    return kLa * (C_CO2_sat - C_CO2)  # [mol m⁻³ s⁻¹]

# 31. H₂ dissolution rate
def n_H2_diss(J_H2_diss, V_catholyte):
    return J_H2_diss * V_catholyte  # [mol s⁻¹]

# 32. CO₂ dissolution rate
def nCO2_diss(J_CO2_diss, V_catholyte):
    return J_CO2_diss * V_catholyte  # [mol s⁻¹]

# 33. H₂ venting rate
def n_H2_vent(kvent_H2, n_H2_bubble):
    return kvent_H2 * n_H2_bubble  # [mol s⁻¹]

# 34. CO₂ venting rate
def n_CO2_vent(kvent_CO2, nCO2_bubble):
    return kvent_CO2 * nCO2_bubble  # [mol s⁻¹]

# 35. Gas venting constant
def kvent_H2(vb_H2, Hreactor):
    return vb_H2 / Hreactor
def kvent_CO2(vb_CO2, Hreactor):
    return vb_CO2 / Hreactor

# 36. Bubble rise velocity (Stokes' law)
def vb_H2(rhoL, rhoG_H2, g, d_b_H2, mu):
    return (rhoL - rhoG_H2) * g * d_b_H2**2 / (18 * mu)
def vb_CO2(rhoL, rhoG_CO2, g, d_b_CO2, mu):
    return (rhoL - rhoG_CO2) * g * d_b_CO2**2 / (18 * mu)

# 37. Electron flux from cathode to biofilm - (DET)
def ne_DET(I, F, A_cat, f_DET):
    return f_DET * I / (F * A_cat)   # [mol e⁻ m⁻² s⁻¹]

# 38. CO₂ uptake rate (diffusion-limited, biofilm)
def nCO2_DET(kL_CO2, Abiofilm, C_CO2_sat, C_CO2):
    # Attach units if C_CO2 is just a float
    if not hasattr(C_CO2, 'to'):
        C_CO2 = ureg.Quantity(C_CO2, 'mol/meter**3')
    if not hasattr(C_CO2_sat, 'to'):
        C_CO2_sat = ureg.Quantity(C_CO2_sat, 'mol/meter**3')
    return kL_CO2 * Abiofilm * (C_CO2_sat - C_CO2)  # [mol s⁻¹]

# 39a. Specific electron uptake rate (biofilm)
def q_e(ne_DET, X_B):
    return ne_DET / X_B  # [mol e⁻ m⁻² (gDW)⁻¹ s⁻¹]

# 39b. Specific CO₂ uptake rate (biofilm)
def q_CO2_B(nCO2_DET_dot, X_B):
    return nCO2_DET_dot / X_B  # [mol CO₂ (C-mol X)⁻¹ s⁻¹]

# 40. CO₂ uptake by planktonic cells (Monod)from pint import Quantity

def n_CO2_P(X_P, V_catholyte, Vmax_CO2, KCO2, C_CO2):
    # Attach units if C_CO2 is just a float
    if not hasattr(C_CO2, 'to'):
        C_CO2 = ureg.Quantity(C_CO2, 'mol/meter**3')
    if not hasattr(KCO2, 'to'):
        KCO2 = ureg.Quantity(KCO2, 'mol/meter**3')
    if not hasattr(Vmax_CO2, 'to'):
        Vmax_CO2 = ureg.Quantity(Vmax_CO2, '1/second')
    # Do the calculation
    return V_catholyte * X_P * ((Vmax_CO2 * C_CO2) / (C_CO2 + KCO2))


# 41. H₂ uptake by planktonic cells
def r_H2(kA, C_H2, C_H2_sat):
    return kA * (C_H2_sat - C_H2)  # [mol H₂ m⁻³ s⁻¹]

# 42a. Specific H₂ uptake rate (planktonic)
def q_H2(r_H2, X_P):
    return r_H2 / X_P  # [mol H₂ (gDW)⁻¹ s⁻¹]

# 42b. Specific CO₂ uptake rate (planktonic)
def q_CO2_P(n_CO2_P, X_P):
    return n_CO2_P / X_P  # [mol CO₂ (C-mol X)⁻¹ s⁻¹]

# 43. Growth rate of biofilm cells
def r_growth_biofilm(ne_DET, nCO2_DET_dot, A_cat):
    return min((ne_DET * A_cat) / 4.2, nCO2_DET_dot / 1)  # [mol X s⁻¹]

# 44. Acetate production rate by biofilm
def r_acetate_biofilm(YAc_X_B, r_growth_biofilm):
    return YAc_X_B * r_growth_biofilm  # [mol s⁻¹]

# 45. Growth rate of planktonic cells
def r_growth_plank(r_H2, nCO2_P_dot):
    return min(r_H2 / 2.1, nCO2_P_dot / 1)  # [mol X s⁻¹]

# 46. Acetate production rate by planktonic cells
def r_acetate_plank(YAc_X_P, r_growth_plank):
    return YAc_X_P * r_growth_plank  # [mol s⁻¹]
