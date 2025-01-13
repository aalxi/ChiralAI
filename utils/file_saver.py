import pandas as pd
from datetime import datetime

def save_suggestions_to_csv(suggestions):
    """
    Saves a list of molecule suggestions to a timestamped CSV file.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"suggestions_{timestamp}.csv"

    flattened_data = []
    for suggestion in suggestions:
        flat_dict = {}
        for key, value in suggestion.items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flat_dict[f"{key}_{sub_key}"] = str(sub_value)
            else:
                flat_dict[key] = str(value)
        flattened_data.append(flat_dict)

    df = pd.DataFrame(flattened_data)
    df.to_csv(filename, index=False)
    return filename
