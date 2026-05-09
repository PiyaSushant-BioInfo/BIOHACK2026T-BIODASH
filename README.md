Create README.md in the project root with exactly this content:

# TBAnalytica 🧬
### Real-Time TB Variant Analysis & Treatment Decision Support System

TBAnalytica is a clinical bioinformatics tool that identifies 
Mycobacterium tuberculosis variants, analyzes drug resistance 
profiles, and generates treatment recommendations in real time 
using live data from NCBI, UniProt, WHO mutation catalogue, 
and AlphaFold.

Built for Nepal's TB crisis — where MDR-TB cases increased 
from 14.6% (2015) to 59% (2020).

---

## 🚀 Quick Start

### Prerequisites
- Python 3.11 or higher
- Git

### Installation

# 1. Clone the repository
```
git clone https://github.com/PiyaSushant-BioInfo/BIOHACK2026T-BIODASH.git
cd BIOHACK2026T-BIODASH
```

# 2. Install dependencies
```
pip install -r requirements.txt
```
# 3. Set up environment variables
```
cp .env.example .env
# Edit .env and add your NCBI API key (free at https://www.ncbi.nlm.nih.gov/account/)
```
# 4. Run the web interface
```
python app.py
```

# 5. Open browser
# Navigate to http://localhost:5000

---

## 💻 Usage

### Web Interface 
```
python app.py
```
Then open http://localhost:5000

### Command Line
```
python main.py
```
---

## 🔬 Features

- **Known Variant Analysis** — Full MDR profile, risk score, 
  and treatment recommendation for all major TB lineages

- **New Variant Analysis** — Live NCBI/UniProt sequence fetch, 
  dual-layer protein + gene comparison against known variants

- **Silent Mutation Detection** — Flags precursor mutations 
  one step away from resistance

- **Data Quality Scoring** — Every data point scored 0-100 
  for source credibility, citation support, completeness, recency

- **Similarity Inference** — 
  >90% match → same treatment protocol
  60-90% match → modified protocol with caution flags  
  <60% match → novel variant alert, escalate immediately

- **Clinical Reports** — Doctor version (technical) + 
  Patient version (plain language)

- **Visualizations** — Risk gauge, mutation map, 
  similarity chart, resistance heatmap

---

## 📁 Project Structure

TBAnalytica/
├── app.py                  # Flask web interface
├── main.py                 # CLI + core orchestration
├── schema.py               # Shared data schemas
├── modules/
│   ├── variant_db.py       # Known variant database
│   ├── api_calls.py        # NCBI, UniProt, PDB, AlphaFold
│   ├── protein_compare.py  # Protein similarity analysis
│   ├── gene_analysis.py    # Nucleotide mutation analysis
│   ├── similarity_score.py # Weighted scoring engine
│   ├── treatment.py        # Treatment recommendation
│   ├── report_generator.py # PDF report generation
│   └── charts.py           # Visualizations
├── templates/              # Flask HTML templates
├── static/                 # CSS and JS
├── data/
│   ├── known_variants.json
│   ├── resistance_mutations.json
│   └── treatment_protocols.json
├── output/                 # Generated reports and charts
├── .env.example
└── requirements.txt

---

## 🔑 API Keys

TBAnalytica uses the following free APIs:

| API | Required | Get Key |
|-----|----------|---------|
| NCBI Entrez | Recommended | https://www.ncbi.nlm.nih.gov/account/ |
| UniProt | No key needed | Free public API |
| PDB | No key needed | Free public API |
| AlphaFold | No key needed | Free public API |

NCBI works without a key but is rate-limited to 3 
requests/second. With a free key: 10 requests/second.

Add to your .env file:
NCBI_API_KEY=your_key_here
NCBI_EMAIL=your_email@example.com

---

## 📊 TB Variants in Database

| Lineage | Family | Nepal Prevalence |
|---------|--------|-----------------|
| Lineage 2 | Beijing | 48.4% |
| Lineage 3 | CAS/Delhi | 30.7% |
| Lineage 4 | Euro-American | 14.5% |
| Lineage 1 | Indo-Oceanic | 6.4% |
| MDR-TB strains | Various | - |
| XDR-TB strains | Various | - |

---

## 🧪 Resistance Genes Analyzed

| Gene | Drug Target | Resistance |
|------|-------------|------------|
| rpoB | Rifampicin | MDR-TB |
| katG | Isoniazid | MDR-TB |
| inhA | Isoniazid | MDR-TB |
| gyrA | Fluoroquinolones | XDR-TB |
| gyrB | Fluoroquinolones | XDR-TB |
| pncA | Pyrazinamide | MDR-TB |
| embB | Ethambutol | MDR-TB |

---

## ⚠️ Disclaimer

TBAnalytica is a decision SUPPORT tool — not a replacement 
for clinical judgment. All recommendations should be verified 
by a qualified clinician. Treatment decisions must comply with 
local guidelines (Nepal NTP / WHO).

---

## 📖 Scientific Basis

Built on:
- WHO TB Mutation Catalogue 2022
- Nepal MDR-TB surveillance data (TUTH, Kathmandu)
- NCBI RefSeq Mycobacterium tuberculosis reference genome
- UniProt/SwissProt reviewed TB protein entries

Key references:
- Shrestha et al. (2018) — MDR-TB lineage distribution Nepal
- WHO Global TB Report 2023
- CRyPTIC Consortium — rpoB resistance variants

---

## 🤝 Contributing

Pull requests welcome. For major changes please open 
an issue first.

---

## 📄 License

MIT License — see LICENSE file for details.

---

*Built for the TB crisis in Nepal 🇳🇵*
*Because when the bacteria evolves, your diagnostic 
tools need to evolve faster.*
