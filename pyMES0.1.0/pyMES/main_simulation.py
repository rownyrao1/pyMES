import numpy as np
from utils.pint_registry import PintRegistry, ureg
from electrochemistry import ElectrochemistryModule
from transport import TransportModule
from gas_dynamics import GasDynamicsModule
from microbial_uptake import MicrobialUptakeModule
from growth_production import GrowthProductionModule
import equations

def to_quantity(val, unit):
    if hasattr(val, "to"):
        return val.to(unit)
    return ureg.Quantity(val, unit)

# === 1. Load Registry from Excel ===
reg = PintRegistry.from_excel("utils/variables.xlsx")

# === 2. Run Electrochemistry Module ===
ec = ElectrochemistryModule(reg, mode="ca", applied_value=-1.0)
ec_results = ec.run()  # Must provide a 'current' key

# === 3. Initialize Transport Module (if needed) ===
tm = TransportModule(reg, ec_results=ec_results)
tm_results = tm.step(3600)  # example step (1 hour)

# === 4. Initialize Gas Dynamics Module ===
gd = GasDynamicsModule(reg, p_max_headspace=1.0)

# === 5. Initialize Microbial Uptake Module ===
uptake_module = MicrobialUptakeModule(reg)

# === 6. Initialize Growth and Production Module ===
growth_module = GrowthProductionModule(reg)

# === 7. Simulation Settings ===
dt = 60 * 60  # 1 hour in seconds
n_steps = 24  # 1 day, hourly

# Generation and feed rates
H2_gen = tm_results['H2_production']             # mol/s from transport
CO2_feed = reg['flow_rate']                      # m³/s (gas flow)
CO2_conc = reg['CO2_conc']                       # mol/m³ (CO₂ in feed)

print(
    f"{'Hour':>4} | {'C_H2(diss)':>10} | {'n_H2_head':>10} | {'C_CO2(diss)':>10} | {'n_CO2_head':>10} | "
    f"{'p_H2(atm)':>10} | {'p_CO2(atm)':>10} | {'pTot(atm)':>10} | {'Biofilm_X':>10} | {'Plank_X':>10} | {'Acetate':>10}"
)

for step in range(n_steps):
    # === A. Get current concentrations and bubble pools for uptake calculation ===
    C_H2 = gd.C_H2
    C_CO2 = gd.C_CO2
    v = gd.vars

    # Calculate saturation values using current bubble pools
    p_H2 = equations.p_H2(gd.n_H2_head, v["R"], v["T"], v["Vhead"])
    p_CO2 = equations.pCO2(gd.n_CO2_head, v["R"], v["T"], v["Vhead"])
    C_H2_sat = equations.C_H2_sat(v["H_H2"], p_H2)
    C_CO2_sat = equations.C_CO2_sat(v["H_CO2"], p_CO2)

    # === B. Calculate dynamic uptake rates using microbial uptake module ===
    uptake = uptake_module.step(
        dt=dt,
        I=-ec_results["current"],
        C_CO2=C_CO2,
        C_CO2_sat=C_CO2_sat,
        C_H2=C_H2,
        C_H2_sat=C_H2_sat
    )
    nCO2_DET = to_quantity(uptake.get('nCO2_DET', 0), 'mol/second')
    n_CO2_P  = to_quantity(uptake.get('n_CO2_P',  0), 'mol/second')
    n_H2     = to_quantity(uptake.get('n_H2', 0), 'mol/second')
    n_electrons = to_quantity(uptake.get('n_electrons', 0), 'mol/s')
    q_e_biofilm = to_quantity(uptake.get('q_e_biofilm', 0), '1/s')
    q_CO2_B = to_quantity(uptake.get('q_CO2_B', 0), '1/s')
    q_H2_P = to_quantity(uptake.get('q_H2_P', 0), '1/s')
    q_CO2_P = to_quantity(uptake.get('q_CO2_P', 0), '1/s')

    # Non-negativity enforcement
    if nCO2_DET.magnitude < 0:
        nCO2_DET = 0 * ureg('mol/second')
    if n_CO2_P.magnitude < 0:
        n_CO2_P = 0 * ureg('mol/second')
    if n_H2.magnitude < 0:
        n_H2 = 0 * ureg('mol/second')

    CO2_uptake = (nCO2_DET + n_CO2_P).to('mol/second')
    H2_uptake = n_H2

    A_cat = v["A_cat"]   # already loaded as a Pint Quantity from gas dynamics vars
    ne_DET = n_electrons / A_cat  # [mol/s] / [m²] = [mol/(m²·s)]
    # Now call growth module:
    growth_result = growth_module.step(
        dt=ureg.Quantity(dt, 'second'),
        ne_DET=ne_DET,
        nCO2_DET=nCO2_DET,
        r_H2=n_H2,
        nCO2_P=n_CO2_P
    )

    # Update uptake module with new biomass for next step
    uptake_module.set_biomass(growth_result["X_B"], growth_result["X_P"])

    # === D. Gas Dynamics Step ===
    res = gd.step(
        dt,
        H2_production=H2_gen,
        CO2_flow_rate=CO2_feed,
        CO2_conc=CO2_conc,
        H2_uptake=H2_uptake,
        CO2_uptake=CO2_uptake
    )

    # === E. Output and Diagnostics ===
    p_H2_atm = res['p_H2_head'].to('atm').magnitude
    p_CO2_atm = res['p_CO2_head'].to('atm').magnitude
    pTot = (res['p_H2_head'] + res['p_CO2_head']).to('atm').magnitude

    print(
        f"{step+1:4d} | {res['C_H2'].magnitude:10.3e} | {res['n_H2_head'].magnitude:10.3e} | "
        f"{res['C_CO2'].magnitude:10.3e} | {res['n_CO2_head'].magnitude:10.3e} | "
        f"{p_H2_atm:10.3e} | {p_CO2_atm:10.3e} | {pTot:10.3e} | "
        f"{growth_result['X_B'].magnitude:10.3e} | {growth_result['X_P'].magnitude:10.3e} | {growth_result['acetate'].magnitude:10.3e}"
    )
