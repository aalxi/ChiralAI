#!/usr/bin/env python
"""Which enzyme makes THIS enantiomer?

  python -m stereo.query --smiles 'C[C@H](O)C(=O)O'
  python -m stereo.query --name limonene
  python -m stereo.query --smiles '...' --sequences

A thin command line over ``stereo.api.lookup``. The matching logic lives in
the library so the pipeline and this CLI cannot drift apart.
"""
import argparse, csv, sys, urllib.parse, urllib.request

from rdkit import Chem

from stereo.api import lookup, _index

LABEL = {
    "EXACT":      "EXACT      — a curated reaction acts on this exact stereoisomer.",
    "WRONG-FORM": "WRONG-FORM — this scaffold is curated, but not your stereoisomer.",
    "SKELETON":   "SKELETON   — scaffold curated, stereochemistry unspecified in your query.",
    "NONE":       "NONE       — no curated reaction matches this connectivity.",
}


def uniprot_for(rhea_ids, limit=5):
    out = []
    for rid in rhea_ids[:3]:
        q = urllib.parse.quote(f'(rhea:{rid.replace("RHEA:", "")})')
        url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
               f"&fields=accession,id,protein_name,organism_name,ec&format=tsv&size={limit}")
        try:
            with urllib.request.urlopen(url, timeout=40) as fh:
                lines = fh.read().decode().strip().split("\n")
            for ln in lines[1:]:
                p = ln.split("\t")
                if len(p) >= 4:
                    out.append((rid, p[0], p[2][:52], p[3][:34]))
        except Exception as exc:
            out.append((rid, "ERROR", str(exc)[:40], ""))
    return out


def report(smiles: str, want_seq: bool) -> int:
    r = lookup(smiles)
    if r["outcome"] == "ERROR":
        print("could not parse that SMILES")
        return 2

    print(f"target            : {Chem.MolToSmiles(Chem.MolFromSmiles(smiles))}")
    print(f"InChIKey          : {r['inchikey']}")
    print(f"stereocentres     : {r['n_stereocentres']} total, {r['n_defined']} defined")
    if r["n_stereocentres"] and not r["n_defined"]:
        print("\nWARNING: your structure has stereocentres but none are specified.")
        print("         A biosynthetic route cannot be selected without that.")
    print(f"\nRESULT: {LABEL[r['outcome']]}\n")

    if r["outcome"] == "EXACT":
        print(f"  compound      : {r['compound']}")
        print(f"  makes it      : {', '.join(r['enzymes_for_requested_form'])}")
        print(f"  Rhea          : {';'.join(r['rhea_ids'])[:120]}")
        risk = r["selectivity_risk"]
        if risk:
            print("\n  SELECTIVITY RISK — the enantiomer is also enzymatically accessible:")
            print(f"    {risk['enantiomer']}")
            print(f"    made by     : {', '.join(risk['made_by_ec']) or 'see pair table'}")
            print("    Choosing between these two EC families IS the stereochemical decision.")
        else:
            print("\n  No enantiocomplementary enzyme is recorded in Rhea for this scaffold.")
        if want_seq:
            print("\n  UniProt entries:")
            for rid, acc, nm, org in uniprot_for(r["rhea_ids"]):
                print(f"    {rid:12s} {acc:10s} {nm:52s} {org}")
        return 0

    for e in r["enzymes_for_other_forms"]:
        print(f"    {e['compound'][:58]:58s} EC {', '.join(e['ec'])[:34]}")
    print(f"\n  {r['detail']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles")
    ap.add_argument("--name")
    ap.add_argument("--sequences", action="store_true", help="fetch UniProt accessions")
    a = ap.parse_args()
    if a.name and not a.smiles:
        rows, _, _ = _index()
        hits = [x for x in rows if a.name.lower() in x["name"].lower()]
        if not hits:
            print(f"no indexed compound matching '{a.name}'")
            return 1
        print(f"{len(hits)} indexed compound(s) matching '{a.name}':\n")
        for x in hits[:12]:
            flag = " [enantiocomplementary pair known]" if x["enantiocomplementary_pair_known"] == "True" else ""
            print(f"  {x['name'][:44]:44s} {x['chebi']:14s} EC {x['ec_acting_on_this_isomer'][:30]}{flag}")
            print(f"      {x['smiles'][:84]}")
        return 0
    if not a.smiles:
        ap.print_help()
        return 1
    return report(a.smiles, a.sequences)


if __name__ == "__main__":
    sys.exit(main())
