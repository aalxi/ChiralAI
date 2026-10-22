"""ChiralAI CLI.

  python main.py --query "chiral monomer for biodegradable polyester"
  python main.py --smiles 'C[C@H](O)C(=O)O' --name 'L-lactic acid'
  python main.py --smiles '...' --host yeast --thermo

Stage A (the language model) is optional and is never trusted: every structure
it proposes is checked against KEGG before it is assessed. Stage B is
database-grounded throughout, and says so when a layer is unavailable.
"""
from __future__ import annotations
import argparse, json, sys

from ChiraLLM.feasibility import DEFAULT_MAX_MW, assess
from utils.display import display_results_in_terminal
from utils.file_saver import save_data_to_csv


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ChiralAI: stereochemistry-aware biosynthetic assessment")
    p.add_argument("--query", help="natural-language request; uses the LLM discovery stage")
    p.add_argument("--smiles", help="assess one structure directly, skipping the LLM")
    p.add_argument("--name", help="name for --smiles (used only to cross-check the structure)")
    p.add_argument("--kegg-id", help="KEGG compound id to cross-check against, if known")
    p.add_argument("--host", default="E. coli", help="host organism (default: E. coli)")
    p.add_argument("--any-pathway", action="store_true",
                   help="accept any KEGG pathway as evidence, not just biosynthesis maps. "
                        "Reproduces the old behaviour; measured 13 points less specific")
    p.add_argument("--max-mw", type=float, default=DEFAULT_MAX_MW,
                   help=f"molecular-weight flag in Da (default {DEFAULT_MAX_MW:g})")
    p.add_argument("--thermo", action="store_true", help="compute reaction ΔG'° (needs equilibrator-api)")
    p.add_argument("--brenda", action="store_true", help="consult a local BRENDA dump if present")
    p.add_argument("--no-cross-check", action="store_true",
                   help="skip structure/name formula agreement (not recommended)")
    p.add_argument("--json", action="store_true", help="print the full records as JSON")
    return p


def main() -> int:
    args = build_parser().parse_args()
    opts = dict(require_biosynthesis=not args.any_pathway, max_mw=args.max_mw,
                with_thermo=args.thermo, with_brenda=args.brenda,
                cross_check=not args.no_cross_check)

    if args.smiles:
        records = [assess(args.name or args.smiles, args.smiles, args.host,
                          cited_kegg_id=args.kegg_id, **opts)]
    else:
        query = args.query or input("Natural-language query: ")
        from ChiraLLM.query_handler import ask_gpt_chirality, parse_suggestions
        try:
            suggestions = parse_suggestions(ask_gpt_chirality(query))
        except Exception as exc:
            print(f"Stage A failed → {type(exc).__name__}: {exc}", file=sys.stderr)
            return 1
        print(f"Model proposed {len(suggestions)} molecule(s); each will be "
              f"structure-checked before assessment.\n")
        saved = save_data_to_csv(suggestions, "suggestions")
        print(f"Raw suggestions saved to {saved}")
        records = [
            assess(s.get("name") or s.get("Name") or "?",
                   s.get("smiles") or s.get("SMILES") or "",
                   args.host, cited_kegg_id=s.get("kegg_id") or s.get("KEGG_ID"), **opts)
            for s in suggestions
        ]

    if args.json:
        print(json.dumps(records, indent=2, default=str))
    else:
        display_results_in_terminal(records)
    save_data_to_csv([{k: v for k, v in r.items() if not isinstance(v, (dict, list))}
                      for r in records], f"feasibility_{args.host.replace(' ', '_')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
