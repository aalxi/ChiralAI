import pandas as pd
from datetime import datetime
import os

def save_data_to_csv(data: list[dict], base_filename: str) -> str:
    """
    Saves a list of dictionaries to a timestamped CSV file.
    The filename will be 'output/{base_filename}_{timestamp}.csv'.
    """
    if not data:
        return "No data to save."

    # Create an 'output' directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{base_filename}_{timestamp}.csv"
    full_path = os.path.join(output_dir, filename)

    # Use pandas to easily handle list of dicts and save to CSV
    df = pd.DataFrame(data)
    df.to_csv(full_path, index=False, encoding='utf-8')

    return full_path