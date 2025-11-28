# pyMES

**Version:** 0.1.0  
**Authors:** Ravineet Yadav, Mohammed Qasim, Dr. Siddharth Gadkari and Dr. Sunil A. Patil

---

## Overview

**pyMES** is an integrated Python framework for multiscale modeling and optimization of Microbial Electrosynthesis (MES) systems.  
It provides modular, customizable, and extensible tools to simulate and analyze the coupled electrochemical, gas, microbial, and transport phenomena in MES from CO2.

---

## Features

- **Electrochemistry module:** Simulate potentiostatic/galvanostatic operation, compute current/voltage profiles, Tafel analysis, and system resistances.
- **Transport module:** Proton transport, buffering, and dynamic pH balance in anode/cathode chambers.
- **Gas dynamics module:** Gas phase (H2/CO2) dissolution, bubble dynamics, mass transfer, and headspace constraints.
- **Microbial uptake module:** Calculate electron, H2, and CO2 uptake by biofilm and planktonic populations.
- **Growth/production module:** Model microbial growth, decay, acetate production, and biomass partitioning.
- **Registry-based variable management:** Flexible, Excel/CSV-driven parameter registry.
- **Comprehensive utilities:** Logging, exporting, plotting, and registry handling.
- **Extensible architecture:** Easily add new modules or modify existing equations.

---

## Installation

> **Requirements:**  
> - Python 3.7+  
> - numpy, scipy, pandas, pint, matplotlib

## Usage Example

```bash
import pyMES

# Load your variable registry (see "Variable Registry" below)

from pyMES.utils.pint_registry import PintRegistry

reg = PintRegistry.from_excel("input_registry.xlsx")

# Initialize modules
ec = pyMES.ElectrochemistryModule(reg, mode="potentiostatic", applied_value=-1.0)  # -1.0 V
gd = pyMES.GasDynamicsModule(reg)
mu = pyMES.MicrobialUptakeModule(reg)
gp = pyMES.GrowthProductionModule(reg)

# Run electrochemistry simulation
ec_results = ec.run()
print(ec_results)

# List available modules
print(pyMES.list_available_modules())

# Access a module dynamically
TransportCls = pyMES.get_module_class("transport")
tm = TransportCls(reg, ec_results=ec_results)
```

## Variable Registry
pyMES relies on a parameter registry (Excel or CSV) with columns:

**Variable:** Name of variable (must match code)

**Value:** Initial value (number, string, or dict)

**Unit:** Units (compatible with pint)

**Description:** Short description

**Min/Max:** (Optional) Valid range

## Example (input_registry.xlsx):
| Variable    | Value   | Unit | Description                    | Min        | Max        |
|-------------|---------|------|-------------------------------|------------|------------|
| pH_cat      | 7       | 1    | pH at cathode                  | 5.00E+00   | 9.00E+00   |
| pH_an       | 2.5     | 1    | pH at anode                    | 2.00E+00   | 9.00E+00   |
| d_an_memb   | 0.03    | m    | Anode-membrane distance        | 1.00E-03   | 1.00E-01   |
| A_an_cs     | 0.00005 | m²   | Anode area (cross sectional)   |            |            |
| A_an        | 0.00125 | m²   | Anode area                     | 1.00E-04   | 1.00E+03   |
| d_cat_memb  | 0.03    | m    | Cathode-membrane distance      | 1.00E-03   | 1.00E-01   |

## Extending pyMES
> **Add a New Module**
> - Create a new .py file in pyMES/ (e.g. my_module.py).
>
> - Define your class or functions.
>
> - Import it in pyMES/__init__.py, and add to MODULE_REGISTRY and __all__.
> 
> **Add a New Function/Equation**
> - Add your function to the appropriate file (e.g. equations.py), document it, and use as needed.
>
> **Add Variables**
> - Add new variables to your registry Excel/CSV with all required metadata.
>
> **Testing**
> - Add or update test cases in the tests/ directory.
>
> - Run tests as appropriate.
>
> **Documentation**
> - All core modules and utilities are documented with Python docstrings.
>
> - See inline documentation for usage, units, and parameters.
>
>For advanced usage or contributions, see CONTRIBUTING.md (if provided).
>

```bash
python -m venv venv
source venv/bin/activate    # or venv\Scripts\activate on Windows
pip install -r requirements.txt
```
> - Make changes and run your tests in tests/.
>
> - Submit issues or PRs via GitHub.

## License
MIT License

## Acknowledgments
> **Developed by Ravineet Yadav and contributors.**
>
> For academic, research, and educational use.
>
> Inspired by open research in electrochemical and microbial process modeling.

## Contact
For queries, feature requests, or collaborations, reach out via GitHub Issues or email the author.