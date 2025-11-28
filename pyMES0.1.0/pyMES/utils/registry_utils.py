import pandas as pd

def create_registry_from_excel(file_path):
    """
    Create a variable registry from an Excel or CSV file.

    Args:
        file_path (str): Path to the Excel/CSV file.
    Returns:
        dict: Registry dictionary keyed by variable name.
    """
    df = pd.read_excel(file_path)
    registry = {}
    for _, row in df.iterrows():
        registry[row['Variable']] = {
            'value': row['Value'],
            'unit': row['Unit'],
            'description': row['Description'],
            'min': row.get('Min'),
            'max': row.get('Max')
        }
    return registry

def update_registry_from_excel(registry, user_file_path):
    """
    Update the variable registry with values from a user-edited file.

    Args:
        registry (dict): Existing registry dictionary.
        user_file_path (str): Path to updated Excel/CSV file.
    Returns:
        dict: Updated registry.
    """
    user_df = pd.read_excel(user_file_path)
    for _, row in user_df.iterrows():
        var = row['Variable']
        user_val = row['Value']
        # Only update if variable exists and user provided a value
        if var in registry and pd.notnull(user_val):
            min_val = registry[var].get('min')
            max_val = registry[var].get('max')
            # Validate
            if min_val is not None and user_val < min_val:
                print(f"Warning: {var} below min. Setting to min ({min_val}).")
                user_val = min_val
            if max_val is not None and user_val > max_val:
                print(f"Warning: {var} above max. Setting to max ({max_val}).")
                user_val = max_val
            registry[var]['value'] = user_val
    return registry

def registry_to_dataframe(registry):
    """
    Convert the registry dict back to a DataFrame (for exporting or review).
    """
    data = []
    for var, props in registry.items():
        data.append({
            'Variable': var,
            'Value': props['value'],
            'Unit': props['unit'],
            'Description': props['description'],
            'Min': props.get('min', ''),
            'Max': props.get('max', '')
        })
    return pd.DataFrame(data)