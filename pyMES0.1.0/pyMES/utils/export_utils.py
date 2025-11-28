import pandas as pd
import json
import os

def export_results_to_csv(results, file_path, append=False):
    """
    Export results to a CSV file.
    Args:
        results (dict or list of dict): Model results.
        file_path (str): Path to output CSV file.
        append (bool): If True, appends to existing file (if any).
    """
    df = pd.DataFrame(results if isinstance(results, list) else [results])
    if append and os.path.exists(file_path):
        df.to_csv(file_path, mode='a', header=False, index=False)
    else:
        df.to_csv(file_path, index=False)

def export_results_to_excel(results, file_path):
    """
    Export results to an Excel file (.xlsx).
    Args:
        results (dict or list of dict): Model results.
        file_path (str): Path to output Excel file.
    """
    df = pd.DataFrame(results if isinstance(results, list) else [results])
    df.to_excel(file_path, index=False)

def export_results_to_json(results, file_path):
    """
    Export results to a JSON file.
    Args:
        results (dict or list of dict): Model results.
        file_path (str): Path to output JSON file.
    """
    with open(file_path, 'w') as f:
        json.dump(results, f, indent=2)

