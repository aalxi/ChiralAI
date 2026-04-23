import requests

def query_kegg(compound_id):
    """
    Queries KEGG for compound data and parses the flat-file response.
    """
    url = f"http://rest.kegg.jp/get/{compound_id}"
    response = requests.get(url)

    if response.status_code != 200:
        return {"status": "error", "message": f"KEGG returned {response.status_code} for {compound_id}"}

    return _parse_kegg_flat_file(response.text)


def _parse_kegg_flat_file(text):
    """
    Parses a KEGG flat-file response into a structured dict.

    KEGG flat-file format: field labels start at column 0, continuation
    lines are indented. Multiple values within a field are whitespace-separated
    or semicolon-separated depending on field type.
    """
    fields = {"ENTRY": [], "NAME": [], "FORMULA": [], "PATHWAY": [], "ENZYME": [], "REACTION": []}
    current_field = None

    for line in text.splitlines():
        if line.startswith("///"):
            break
        if not line:
            continue

        if not line[0].isspace():
            # New field — label occupies the first 12 characters
            label = line.split()[0]
            value = line[12:].strip()
            if label in fields:
                current_field = label
                if value:
                    fields[current_field].append(value)
            else:
                current_field = None
        elif current_field is not None:
            # Continuation line
            fields[current_field].append(line.strip())

    # Extract ENTRY id (first token of first value)
    entry = fields["ENTRY"][0].split()[0] if fields["ENTRY"] else compound_id

    # NAME: semicolon-separated synonyms may appear on one line
    names = []
    for val in fields["NAME"]:
        names.extend([n.strip().rstrip(";") for n in val.split(";") if n.strip()])

    formula = fields["FORMULA"][0] if fields["FORMULA"] else None

    # PATHWAY lines: "mapXXXXX  Pathway name"
    pathways = []
    for val in fields["PATHWAY"]:
        parts = val.split(None, 1)
        pathways.append({"id": parts[0], "name": parts[1] if len(parts) > 1 else ""})

    # ENZYME lines: space-separated EC numbers
    enzymes = []
    for val in fields["ENZYME"]:
        enzymes.extend(val.split())

    # REACTION lines: space-separated reaction IDs
    reactions = []
    for val in fields["REACTION"]:
        reactions.extend(val.split())

    return {
        "status": "success",
        "entry": entry,
        "names": names,
        "formula": formula,
        "pathways": pathways,
        "enzymes": enzymes,
        "reactions": reactions,
    }
