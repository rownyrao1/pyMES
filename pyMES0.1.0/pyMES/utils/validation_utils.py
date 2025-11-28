class Validator:
    """
    Validation utility class for registry and variable checks.
    All methods are static for convenient use.
    """

    @staticmethod
    def validate_presence(registry, required_vars):
        """
        Ensures all required variables are present in the registry.
        Raises ValueError if any are missing.
        """
        missing = [var for var in required_vars if var not in registry]
        if missing:
            raise ValueError(f"Missing required variables: {missing}")

    @staticmethod
    def validate_type(value, expected_type, varname="variable"):
        """
        Checks if value is of the expected type.
        Raises TypeError if not.
        """
        if not isinstance(value, expected_type):
            raise TypeError(f"{varname} must be of type {expected_type}, got {type(value)}.")

    @staticmethod
    def validate_range(value, min_val=None, max_val=None, varname="variable"):
        """
        Validates that value is within [min_val, max_val] (if provided).
        Raises ValueError if not.
        """
        if min_val is not None and value < min_val:
            raise ValueError(f"{varname} ({value}) is below minimum {min_val}.")
        if max_val is not None and value > max_val:
            raise ValueError(f"{varname} ({value}) is above maximum {max_val}.")

    @staticmethod
    def validate_unit(quantity, expected_unit, varname="variable"):
        """
        Validates a pint.Quantity has the expected unit (exact match).
        Raises ValueError if not.
        """
        if not hasattr(quantity, "units"):
            raise ValueError(f"{varname} is not a pint.Quantity.")
        if str(quantity.units) != expected_unit:
            raise ValueError(f"{varname} has unit '{quantity.units}', expected '{expected_unit}'.")

    @staticmethod
    def validate_dimension(quantity, expected_dim, varname="variable"):
        """
        Validates a pint.Quantity has the expected dimensionality (e.g., '[length]').
        Raises ValueError if not.
        """
        if not hasattr(quantity, "check"):
            raise ValueError(f"{varname} is not a pint.Quantity.")
        if not quantity.check(expected_dim):
            raise ValueError(f"{varname} with units '{quantity.units}' does not match dimension '{expected_dim}'.")

# --- Example usage ---
if __name__ == "__main__":
    import pint
    ureg = pint.UnitRegistry()

    registry = {'A': 5, 'B': 2}
    Validator.validate_presence(registry, ['A', 'B'])
    Validator.validate_type(registry['A'], int, "A")
    Validator.validate_range(registry['A'], min_val=0, max_val=10, varname="A")

    # Unit/dimensionality checks (using pint)
    q = 5 * ureg("m^2")
    Validator.validate_unit(q, "meter ** 2", "area")
    Validator.validate_dimension(q, "[length] ** 2", "area")
