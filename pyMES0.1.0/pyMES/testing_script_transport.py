import numpy as np
from utils.pint_registry import PintRegistry
from electrochemistry import ElectrochemistryModule
from transport import TransportModule  # Rename to transport if needed
import equations

# === USER PARAMETERS ===
dt = 3600.0  # time step in seconds (1 hour)
n_steps = 5 * 24  # 5 days, 24 hours per day
gas_t_per_step = dt  # If integrating gas pool ODEs
gas_num_points = 30  # ODE integration resolution

# === 1. Load variable registry from Excel/CSV ===
reg = PintRegistry.from_excel("utils/variables.xlsx")

# === 2. Run the Electrochemistry module for chosen mode and applied value ===
ec = ElectrochemistryModule(reg, mode="ca", applied_value=-1.0)
ec_results = ec.run()

# === 3. Initialize the Transport module with EC results ===
tm = TransportModule(reg, ec_results=ec_results)

# === 4. (OPTIONAL, if your TransportModule has gas pool logic) ===
# gas_state = np.array([0.0, 0.0, 0.0, 0.0])  # [C_H2, n_H2_bubble, C_CO2, n_CO2_bubble]
# If not using gas pools, skip gas_state throughout.

# === 5. Arrays to store simulation results ===
results = []
time_hours = []

# === 6. Main Simulation Loop ===
for step in range(n_steps):
    t_hr = step * dt / 3600
    # Remove gas_state and gas_* if not implemented in your transport module
    result = tm.step(dt)
    results.append(result)
    time_hours.append(t_hr)

    # If using gas pools, update gas_state for next step
    # if "gas_state" in result:
    #     gas_state = result["gas_state"]

    # Print daily summaries
    if step % 24 == 0 or step == n_steps - 1:
        print(f"\n=== Day {step // 24}, Hour {step % 24} ===")
        print(f"pH_an: {result['pH_an']:.2f}, pH_cat: {result['pH_cat']:.2f}")
        print(f"Protons (anolyte): {result['N_H_an']:.3e}, Protons (catholyte): {result['N_H_cat']:.3e}")
        print(f"Buffer base: {result['N_buffer_base']:.3e}, Buffer acid: {result['N_buffer_acid']:.3e}")
        print(f"H2 production: {result['H2_production']:.3e}, Proton transport rate: {result['r_H_trans']:.3e}")
        # If gas tracking is implemented:
        # print(f"H2 (bubble): {result['gas_state'][1]:.3e} mol, CO2 (bubble): {result['gas_state'][3]:.3e} mol")
        # print(f"H2 Pressure: {result.get('p_H2', 0)/101325:.3f} atm, CO2 Pressure: {result.get('p_CO2', 0)/101325:.3f} atm")

# === 7. (Optional) Convert results to arrays/dataframes for plotting ===
# Example:
# import pandas as pd
# df = pd.DataFrame(results)
# df["time_hr"] = time_hours
# df.plot(x="time_hr", y=["pH_an", "pH_cat"])
