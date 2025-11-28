"""
pyMES Initialization Module

Exposes core modules and registry for dynamic access.
"""

__version__ = "0.1.0"
__author__ = "Ravineet Yadav"
__description__ = (
    "pyMES: An Integrated Python Framework for Multiscale Modeling and Optimization of Microbial "
    "Electrosynthesis from CO2, integrating chemical and biological systems."
)

# --- Imports: Expose core module classes directly ---
from .electrochemistry import ElectrochemistryModule
from .transport import TransportModule
from .gas_dynamics import GasDynamicsModule
from .microbial_uptake import MicrobialUptakeModule
from .growth_production import GrowthProductionModule

# Key: module name as string; Value: class reference
MODULE_REGISTRY = {
    "electrochemistry": ElectrochemistryModule,
    "transport": TransportModule,
    "gas_dynamics": GasDynamicsModule,
    "microbial_uptake": MicrobialUptakeModule,
    "growth_production": GrowthProductionModule,
}

def list_available_modules():
    """
    Return a list of available core modules in the framework.
    """
    return list(MODULE_REGISTRY.keys())

def get_module_class(module_name):
    """
    Get the module class by name.
    Usage: get_module_class("electrochemistry")
    """
    return MODULE_REGISTRY.get(module_name, None)

__all__ = [
    "ElectrochemistryModule",
    "TransportModule",
    "GasDynamicsModule",
    "MicrobialUptakeModule",
    "GrowthProductionModule",
    "list_available_modules",
    "get_module_class",
    "__version__",
    "__author__",
    "__description__",
]
