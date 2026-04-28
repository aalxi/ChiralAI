import json
import os
import pandas as pd
from datetime import datetime


def _flatten_scoring(scoring: dict) -> dict:
    te = scoring.get("top_enzyme") or {}
    return {
        "scoring_composite_score": scoring.get("composite_score"),
        "scoring_confidence": scoring.get("confidence"),
        "scoring_top_enzyme_ec": te.get("ec_number"),
        "scoring_top_enzyme_ee": te.get("ee_value"),
        "scoring_top_enzyme_source": te.get("ee_source"),
        "scoring_stereo_confirmed": scoring.get("stereo_confirmed"),
        "scoring_feasibility_flux": scoring.get("feasibility_flux"),
        "scoring_notes": "; ".join(scoring.get("scoring_notes") or []),
    }


def save_suggestions_to_csv(suggestions, out_dir: str = ".") -> tuple[str, str]:
    """
    Saves suggestions to a timestamped CSV (flat columns) and a parallel JSON
    sidecar (full nested structure including enzyme_rankings).

    Returns (csv_path, json_path).
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path  = os.path.join(out_dir, f"suggestions_{timestamp}.csv")
    json_path = os.path.join(out_dir, f"suggestions_{timestamp}.json")

    flattened_data = []
    for suggestion in suggestions:
        flat_dict = {}
        for key, value in suggestion.items():
            if key == "scoring" and isinstance(value, dict):
                flat_dict.update(_flatten_scoring(value))
            elif isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flat_dict[f"{key}_{sub_key}"] = str(sub_value)
            else:
                flat_dict[key] = str(value)
        flattened_data.append(flat_dict)

    df = pd.DataFrame(flattened_data)
    df.to_csv(csv_path, index=False)

    with open(json_path, "w") as f:
        json.dump(suggestions, f, indent=2, default=str)

    return csv_path, json_path
