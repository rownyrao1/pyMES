class Registry:
    """
    Object-oriented wrapper for variable registry.
    Supports scalar and entity-specific variable access.
    """
    def __init__(self, registry_dict):
        self.registry = registry_dict

    @classmethod
    def from_excel(cls, file_path):
        """
        Create a Registry instance directly from an Excel file.
        """
        from registry_utils import create_registry_from_excel
        registry_dict = create_registry_from_excel(file_path)
        return cls(registry_dict)

    def get(self, var, entity=None):
        val = self.registry[var]['value']
        if entity is not None:
            try:
                return val[entity]
            except (KeyError, TypeError):
                raise KeyError(f"Entity '{entity}' not found in '{var}'")
        return val

    def extract(self, var_list):
        """
        Extract multiple scalar variables as a dict.
        """
        return {var: self.get(var) for var in var_list}

    def available_variables(self):
        """
        Return a list of all variable names in the registry.
        """
        return list(self.registry.keys())

    def __getitem__(self, var):
        """
        Dictionary-like access: reg['anode_area']
        """
        return self.get(var)

    def __contains__(self, var):
        return var in self.registry

    def __repr__(self):
        return f"<Registry with {len(self.registry)} variables>"

    def to_dataframe(self):
        """
        Return registry as a pandas DataFrame.
        """
        from registry_utils import registry_to_dataframe
        return registry_to_dataframe(self.registry)
