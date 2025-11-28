import numpy as np
from utils.pint_registry import PintRegistry, ureg
from electrochemistry import ElectrochemistryModule
from transport import TransportModule
from gas_dynamics import GasDynamicsModule
import equations

# === 1. Load Registry from Excel ===
reg = PintRegistry.from_excel("utils/variables.xlsx")

# === 2. Run Electrochemistry/Transport Modules (optional, for input rates) ===
ec = ElectrochemistryModule(reg, mode="ca", applied_value=-1.0)
ec_results = ec.run()
tm = TransportModule(reg, ec_results=ec_results)
tm_results = tm.step(3600)  # for 1 hour, just for example

# === 3. Initialize Gas Dynamics Module ===
# You can change p_max_headspace as needed; e.g., (1.2 * ureg.atm)
gd = GasDynamicsModule(reg, p_max_headspace=1.2)

# === 4. Simulation Settings ===
dt = 60 * 60  # 1 hour in seconds
n_steps = 48  # 24 hours

# Generation and uptake rates
H2_gen = tm_results['H2_production']             # mol/s from transport
CO2_feed = reg['flow_rate']                      # m³/s (gas flow)
CO2_conc = reg['CO2_conc']                       # mol/m³ (CO₂ in feed)
H2_uptake = 0.00000014 * ureg("mol/second")
CO2_uptake = 0.000000025 * ureg("mol/second")

print(
    f"{'Hour':>4} | {'C_H2(diss)':>10} | {'n_H2_head':>10} | {'C_CO2(diss)':>10} | {'n_CO2_head':>10} | "
    f"{'p_H2(atm)':>10} | {'p_CO2(atm)':>10} | {'V_H2(mmol)':>10} | {'V_CO2(mmol)':>11} | {'pTot(atm)':>10}"
)

# === 5. Main Simulation Loop ===
for step in range(n_steps):
    res = gd.step(
        dt,
        H2_production=H2_gen,
        CO2_flow_rate=CO2_feed,
        CO2_conc=CO2_conc,
        H2_uptake=H2_uptake,
        CO2_uptake=CO2_uptake
    )
    p_H2_atm = res['p_H2_head'].to('atm').magnitude
    p_CO2_atm = res['p_CO2_head'].to('atm').magnitude
    pTot = (res['p_H2_head'] + res['p_CO2_head']).to('atm').magnitude

    v_H2 = res['vented_H2'].to('mmol').magnitude
    v_CO2 = res['vented_CO2'].to('mmol').magnitude

    print(
        f"{step+1:4d} | {res['C_H2'].magnitude:10.3e} | {res['n_H2_head'].magnitude:10.3e} | "
        f"{res['C_CO2'].magnitude:10.3e} | {res['n_CO2_head'].magnitude:10.3e} | "
        f"{p_H2_atm:10.3e} | {p_CO2_atm:10.3e} | {pTot:10.3e}"
    )
    # Optionally, print mass balance, etc.
    # mass_bal = gd.mass_balance_check()
    # print(f"Total H2: {mass_bal['total_H2']}, Total CO2: {mass_bal['total_CO2']}")

