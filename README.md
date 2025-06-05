# 🧬 ChiralAI

**From “How” to “What” in Biomanufacturing.**  
_Reimagine bio-based production: Ask for the molecule you want, and ChiralAI maps out how to make it using the power of synthetic biology and AI._

---

## 🚀 What is ChiralAI?

ChiralAI is an open-source platform at the intersection of synthetic biology and artificial intelligence, designed to unlock biology’s untapped potential for producing chiral molecules—compounds where the arrangement of atoms matters, often with huge impact in **pharmaceuticals, agrochemicals, and advanced materials**.

> **“Type ‘chiral nitrogen fertilizer for hydroponics’ and ChiralAI spits out a feasibility heat-map, the CRISPR plasmid design, and—if a step’s missing—an auto-generated enzyme to bridge the gap.”**

**Key Insight:**  
Biology is inherently chiral and excels at making enantiopure compounds. ChiralAI leverages this, using AI to answer not just “how do we make more X?” but “what *new* X’s does biology make best?”

---

## 🧠 Core Modules

### 1. **Discovery Engine** (LLM-powered)
- **Input:** Plain-English queries (e.g., “A biodegradable chiral molecule for crop protection with a cyclopropane ring”)
- **Output:**  
  - Suggested target compounds (known or hypothetical)
  - Biologically plausible precursors/scaffolds
  - Pathway candidates (if known)
  - Relevant literature/patents

### 2. **Feasibility Filter** (Chem/Bio ML models)
- **Predictions:**  
  - Can the molecule be biosynthesized in common hosts (E. coli, yeast)?
  - Likely biosynthetic pathway (e.g., shikimate, mevalonate)
  - Required cofactors, toxic intermediates
  - Known enzymes or functional analogs

### 3. **Optimization Module**
- **Assesses:**  
  - Pathway efficiency and robustness
  - Thermodynamics (e.g., via eQuilibrator-style tools)
  - Enzyme engineering opportunities (using models like ProGen, ESMFold)
  - Codon optimization, regulatory element suggestions

### 4. **Retrospective Benchmarking**
- **Validation:**  
  - Backtests with known biosynthesized chiral molecules
  - Compares synthetic routes (ChiralAI vs. literature)
  - Gathers expert feedback from synthetic biologists

---

## ✨ How to Use ChiralAI

1. **Clone the Repo**
   ```bash
   git clone https://github.com/aalxi/ChiralAI.git
   cd ChiralAI
   ```

2. **Set Up Your Environment**
   - Requires Python 3.8+
   - (Strongly recommended) Set up a virtual environment:
     ```bash
     python -m venv venv
     source venv/bin/activate  # On Windows: venv\Scripts\activate
     ```

3. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```
   _Core dependencies include:_
   - `transformers`, `torch` (LLMs & ML models)
   - `rdkit` (cheminformatics)
   - `networkx` (pathway graphing)
   - `requests`, `tqdm`, etc.

4. **Run the Main Interface**
   ```bash
   python main.py
   ```
   - **OR** launch the web UI (if available):
     ```bash
     streamlit run app.py
     ```

5. **Try Example Queries**
   - _“Chiral building block for β-lactam antibiotics”_
   - _“Enzyme pathway for enantiopure arylpropionates in yeast”_

   The system will return:
   - **Feasibility heatmap**
   - **Biosynthetic pathway sketch**
   - **Enzyme/plasmid suggestions**
   - **Gaps and engineering opportunities**

---

## 📦 Dependencies

- Python 3.8+
- [RDKit](https://www.rdkit.org/) (cheminformatics)
- [Transformers](https://huggingface.co/transformers/) (LLMs)
- [PyTorch](https://pytorch.org/) (ML backend)
- [NetworkX](https://networkx.org/) (graph analysis)
- [Streamlit](https://streamlit.io/) (optional UI)
- See `requirements.txt` for full list.

---

## 🏗️ Project Structure

```
ChiralAI/
├── main.py               # Command-line interface
├── app.py                # Web UI (Streamlit)
├── chiralai/
│   ├── discovery.py      # LLM query parsing & molecule suggestions
│   ├── feasibility.py    # ML models for pathway prediction
│   ├── optimization.py   # Pathway/enzyme optimization routines
│   └── benchmarking.py   # Validation/benchmarking scripts
├── data/                 # Example datasets, pre-built models
├── requirements.txt
└── README.md
```

---

## 🛠️ Long-Term Goals

- **Automated pathway inference** for any chiral target, including novel bio-retrosynthesis.
- **LLM fine-tuning** on chemical patents, enzymatic reactions, and metabolic pathways for domain alignment.
- **Integration with pathway databases** (KEGG, MetaCyc, BRENDA) for richer predictions.
- **Interactive visualization:** Metabolic maps, chiral similarity networks, and pathway heatmaps.
- **Industry relevance:** Focus on high-value classes like:
  - Arylpropionates
  - Cycloalkyl ketones
  - β-lactams
  - Chiral amines and alcohols

---

## 🧪 What Makes ChiralAI Unique?

- **Natural Language to Pathway:** Go from a plain-English chemical wish-list to a lab-ready playbook.
- **Bridges AI & Synthetic Biology:** Fuses LLMs, ML, and metabolic databases for end-to-end design.
- **Designed for extensibility:** Modular codebase for adding new models, organisms, and compound classes.

---

## 💡 Contributing & Collaborating

Interested in helping push the boundaries of biomanufacturing?  
Check out our [CONTRIBUTING.md](CONTRIBUTING.md) or open an issue!  
Want a more detailed onboarding or LLM fine-tuning guide? Let us know.

---

## 🔮 Roadmap & Vision

ChiralAI is just getting started. Our vision is an “oracle” for biosynthetic innovation:  
- Ask for any chiral molecule—get a real, actionable plan for bio-based production.
- Democratize metabolic engineering by making biosynthetic design as intuitive as searching Google.

---

## 👩‍🔬 Developed by BettermindLabs

- **Useful molecular descriptors:** Chiral centers, redox state, MW, solubility, logP, rotatable bonds, H-bond donors/acceptors.
- **Pathway inference:** Integrate with KEGG, MetaCyc, BRENDA; leverage Rhea, UniProt for enzyme mapping.
- **Modeling approaches:** GNNs for compound feasibility, LLMs for text-to-pathway, enzyme design via transformer models.
- **Evaluation:** Expert validation, experimental backtesting, simulated pathway scoring.
- **Visualization:** Interactive metabolic maps, chiral similarity networks.
- **Industry targets:** Prioritize pharmaceuticals, agrochemicals, advanced materials.

---

## 📫 Get in Touch

Questions? Ideas? Join the discussion or reach out via Issues or [email@example.com](mailto:email@example.com).

---

> **“Biology doesn’t ask how to make more of X—it asks what X’s it could make best. ChiralAI lets you ask the same.”**
