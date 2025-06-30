import sys, json, os
from datetime import datetime

from ChiraLLM.query_handler    import ask_gpt_chirality
from utils.file_saver          import save_suggestions_to_csv
from ChiraLLM.feasibility      import process_csv_and_analyze


def main() -> None:
    print("Welcome to ChiraLLM (Discovery Engine of ChiralAI)!")
    query = input("Natural-language query (e.g. 'biodegradable polymer precursor'): ")

    # ─── Stage A: LLM suggestions → raw CSV ──────────────────────────
    raw = ask_gpt_chirality(query)
    print("Raw GPT output:", raw)

    try:
        parsed = json.loads(raw)
        suggestions = [parsed] if isinstance(parsed, dict) else parsed
    except Exception as exc:
        print(f"JSON-parsing failed → {exc}")
        sys.exit(1)

    suggestions_csv = save_suggestions_to_csv(suggestions)
    print(f"✅  Raw suggestions saved to {suggestions_csv}")

    # ─── Stage B: Feasibility analysis from CSV ─────────────────────
    host = input("Host organism for feasibility check (e.g. 'E. coli' or 'yeast'): ").strip()

    analysis_results = process_csv_and_analyze(suggestions_csv, host)
    if isinstance(analysis_results, dict) and analysis_results.get("error"):
        # fatal error (e.g., file not found)
        print("Feasibility step failed →", analysis_results["error"])
        sys.exit(1)

    # Flatten + save annotated output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file  = f"feasibility_{host.replace(' ', '_')}_{timestamp}.csv"
    save_suggestions_to_csv(analysis_results).replace("suggestions", "feasibility")

    # Cheeky way to ensure the filename uses our prefix
    os.rename([f for f in os.listdir() if f.startswith("suggestions_") and f.endswith(".csv")][-1], out_file)

    print(f"🎉  Feasibility results saved to {out_file}")


if __name__ == "__main__":
    main()