from utils.pint_registry import PintRegistry
from electrochemistry import ElectrochemistryModule

reg = PintRegistry.from_excel("utils/variables.xlsx")
ec = ElectrochemistryModule(reg, mode="ca", applied_value=-1.0)
results = ec.run()
print(results)