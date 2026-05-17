"""Curated benchmark set: chiral targets with literature-documented biosynthesis routes.

Each entry encodes ground truth derived from textbook biosynthesis or peer-reviewed
biocatalysis literature. The benchmark scores whether predict_route recovers
recognizable routes against this ground truth.

Selection criteria:
  - Compound is chiral (defined R/S configuration)
  - Compound is NOT itself a central metabolite — i.e., not in
    ChiraLLM.route_predictor.CENTRAL_METABOLITES (which would cause the search
    to terminate at 0 steps because the target is already at the biosynthetic anchor)
  - Compound exists in KEGG with at least one annotated reaction
  - Literature route is documented (textbook, primary paper, or industrial process)
  - Mix of pathway lengths and target types (industrial chiral building blocks,
    pharma intermediates, plant secondary metabolites)

Fields:
  kegg_id              — KEGG compound ID
  name                 — Common chemical name
  category             — Coarse classification (sanity / industrial / pharma /
                         plant_secondary / alkaloid)
  expected_terminals   — Acceptable terminal precursors (KEGG IDs); any match counts as recall
  expected_ec_numbers  — EC numbers known to participate in the literature route
  expected_steps       — Approximate literature step count (KEGG-granularity)
  literature_summary   — One-sentence pathway description
  citation             — Primary source: textbook chapter, paper, or industrial process
"""

BENCHMARK_TARGETS = [
    {
        "kegg_id": "C00186",
        "name": "L-Lactate",
        "category": "sanity",
        "expected_terminals": ["C00022"],  # pyruvate
        "expected_ec_numbers": ["1.1.1.27", "1.1.1.28"],  # L-LDH, D-LDH (some entries)
        "expected_steps": 1,
        "literature_summary": "Pyruvate → L-lactate by L-lactate dehydrogenase (NADH-dependent KRED).",
        "citation": "Berg, Tymoczko, Stryer Biochemistry, 9e — fermentation chapter",
    },
    {
        "kegg_id": "C00599",
        "name": "(R)-Pantolactone",
        "category": "industrial",
        "expected_terminals": ["C00141", "C00966", "C00522"],  # KIV, 2-dehydropantoate, pantoate
        "expected_ec_numbers": ["1.1.1.169", "2.1.2.11"],
        "expected_steps": 3,
        "literature_summary": "KIV → 2-dehydropantoate → (R)-pantoate → (R)-pantolactone (industrial KRED route to pantothenate).",
        "citation": "Hata et al. 1989, Agric Biol Chem 53:1249-1254; Liao et al. 2017, Bioresour Bioprocess",
    },
    {
        "kegg_id": "C00355",
        "name": "L-DOPA",
        "category": "pharma",
        "expected_terminals": ["C00082"],  # L-tyrosine
        "expected_ec_numbers": ["1.14.16.2", "1.14.18.1"],
        "expected_steps": 1,
        "literature_summary": "L-tyrosine → L-DOPA by tyrosine hydroxylase (P450) or tyrosinase; Parkinson's drug precursor.",
        "citation": "Min et al. 2015, Microb Cell Fact 14:71 (E. coli L-DOPA production review)",
    },
    {
        "kegg_id": "C01089",
        "name": "(R)-3-Hydroxybutyrate",
        "category": "industrial",
        "expected_terminals": ["C00024"],  # acetyl-CoA
        "expected_ec_numbers": ["2.3.1.9", "1.1.1.36", "3.1.2.20"],
        "expected_steps": 3,
        "literature_summary": "Acetyl-CoA → acetoacetyl-CoA → (R)-3-hydroxybutyryl-CoA → (R)-3-hydroxybutyrate (PHB monomer pathway).",
        "citation": "Steinbüchel & Lütke-Eversloh 2003, Biochem Eng J 16:81-96",
    },
    {
        "kegg_id": "C01984",
        "name": "(R)-Mandelate",
        "category": "industrial",
        "expected_terminals": ["C02137", "C00601"],  # phenylglyoxylate, benzaldehyde
        "expected_ec_numbers": ["1.1.1.272", "5.1.2.2"],
        "expected_steps": 2,
        "literature_summary": "Phenylglyoxylate → (R)-mandelate by D-mandelate dehydrogenase (Pseudomonas pathway).",
        "citation": "Schmidt & Sieber 1990, Eur J Biochem 191:425-432",
    },
    {
        "kegg_id": "C00418",
        "name": "(R)-Mevalonate",
        "category": "industrial",
        "expected_terminals": ["C00024"],  # acetyl-CoA
        "expected_ec_numbers": ["2.3.1.9", "2.3.3.10", "1.1.1.34"],  # thiolase, HMG-CoA synthase, HMGR
        "expected_steps": 3,
        "literature_summary": "Acetyl-CoA → acetoacetyl-CoA → HMG-CoA → (R)-mevalonate (statin and isoprenoid precursor pathway).",
        "citation": "Goldstein & Brown 1990, Nature 343:425-430 (HMG-CoA reductase classic review)",
    },
    {
        "kegg_id": "C00547",
        "name": "(R)-Noradrenaline",
        "category": "pharma",
        "expected_terminals": ["C00082"],  # L-tyrosine
        "expected_ec_numbers": ["1.14.17.1", "4.1.1.28", "1.14.16.2"],  # DBM, AADC, TH
        "expected_steps": 3,
        "literature_summary": "L-tyrosine → L-DOPA → dopamine → (R)-noradrenaline (catecholamine biosynthesis; sympathetic neurotransmitter).",
        "citation": "Eisenhofer et al. 2004, Pharmacol Rev 56:331-349 (catecholamine metabolism review)",
    },
    {
        "kegg_id": "C00509",
        "name": "(S)-Naringenin",
        "category": "plant_secondary",
        "expected_terminals": ["C00079", "C00423"],  # L-Phe, trans-cinnamate
        "expected_ec_numbers": ["4.3.1.24", "1.14.14.91", "6.2.1.12", "2.3.1.74", "5.5.1.6"],
        "expected_steps": 5,
        "literature_summary": "L-Phe → cinnamate → 4-coumarate → 4-coumaroyl-CoA → naringenin chalcone → (S)-naringenin (general flavonoid pathway).",
        "citation": "Winkel-Shirley 2001, Plant Physiol 126:485-493 (flavonoid biosynthesis review)",
    },
    {
        "kegg_id": "C09136",
        "name": "(S)-Norcoclaurine",
        "category": "alkaloid",
        "expected_terminals": ["C00082"],  # L-tyrosine
        "expected_ec_numbers": ["4.2.1.78", "4.1.1.28", "1.14.16.2"],  # NCS, AADC, TH
        "expected_steps": 4,
        "literature_summary": "L-tyrosine → L-DOPA → dopamine + 4-HPAA → (S)-norcoclaurine (gateway to benzylisoquinoline alkaloids; morphine, codeine).",
        "citation": "Hagel & Facchini 2013, Plant Cell Physiol 54:647-672 (BIA biosynthesis review)",
    },
    {
        "kegg_id": "C00590",
        "name": "Coniferyl alcohol",
        "category": "plant_secondary",
        "expected_terminals": ["C00079"],  # L-Phe
        "expected_ec_numbers": ["4.3.1.24", "1.14.14.91", "1.2.1.44", "1.1.1.195"],  # PAL, C4H, CCR, CAD
        "expected_steps": 6,
        "literature_summary": "L-Phe → cinnamate → ... → feruloyl-CoA → coniferaldehyde → coniferyl alcohol (lignin monomer biosynthesis).",
        "citation": "Boerjan, Ralph, Baucher 2003, Annu Rev Plant Biol 54:519-546 (lignin biosynthesis review)",
    },
]
