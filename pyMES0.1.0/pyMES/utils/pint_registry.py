import pandas as pd
import ast
import pint

ureg = pint.UnitRegistry()

class PintRegistry:
    ureg = pint.UnitRegistry()
    """
    Registry wrapper that returns pint.Quantity for all variables,
    supporting both scalar and multi-entity (dict) values.
    """
    def __init__(self, registry_dict):
        self.registry = registry_dict

    @classmethod
    def from_excel(cls, file_path):
        """
        Create a PintRegistry instance directly from an Excel file.
        """
        registry = {}
        df = pd.read_excel(file_path)
        for _, row in df.iterrows():
            value = row['Value']
            # Parse string dictionaries (for entity-specific variables)
            if isinstance(value, str) and value.strip().startswith("["):
                try:
                    value = ast.literal_eval(value)
                except Exception as e:
                    print(f"Warning: Could not parse dict for {row['Variable']}: {e}")
            registry[row['Variable']] = {
                'value': value,
                'unit': row['Unit'],
                'description': row['Description'],
                'min': row.get('Min'),
                'max': row.get('Max')
            }
        return cls(registry)
    
    def get(self, var, entity=None):
        """
        Returns pint.Quantity with units for the variable.
        If the variable is entity-specific (dict), supply 'entity'.
        If no unit, returns a pint dimensionless quantity.
        """
        try:
            entry = self.registry[var]
        except KeyError:
            raise KeyError(f"Variable '{var}' not found in registry.")
        value = entry['value']
        unit = entry.get('unit', '')
        if entity is not None:
            try:
                value = value[entity]
            except Exception as e:
                raise KeyError(f"Entity '{entity}' not found for '{var}': {e}")

        # Make sure unit is a string; treat None, numbers, or blank as dimensionless
        if not isinstance(unit, str) or unit.strip() == "":
            unit = ""  # This tells pint to treat as dimensionless
        return value * ureg(unit)

    def extract(self, var_list):
        """
        Returns a dict of {var: pint.Quantity} for the given list.
        """
        return {var: self.get(var) for var in var_list}

    def get_dim(self, var, entity=None):
        """Returns the dimensionality string (e.g., '[length]', '[mass]/[time]')"""
        return self.get(var, entity).dimensionality

    def assert_dim(self, var, dim, entity=None):
        """
        Assert a variable's dimension matches 'dim' (string or pint Dimensionality).
        E.g., dim='[length]', '[current]'
        """
        q = self.get(var, entity)
        if not q.check(dim):
            raise ValueError(f"Variable '{var}' with unit '{q.units}' is not dimensionally '{dim}'.")
        return True

    def available_variables(self):
        return list(self.registry.keys())

    def __getitem__(self, var):
        return self.get(var)

    def __contains__(self, var):
        return var in self.registry

    def __repr__(self):
        return f"<PintRegistry with {len(self.registry)} variables>"
    
    def convert_units(value, from_unit, to_unit):
        q = value * ureg(from_unit)
        return q.to(to_unit).magnitude
    