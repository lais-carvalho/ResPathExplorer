# 🧬 ResPathExplorer

A modular Python library for functional analysis of microbial genes using curated bioinformatics databases. This tool supports gene mapping to KEGG pathways, resistance profiles via CARD, and virulence factors through VFDB — providing a comprehensive view of microbial functionality in contexts such as immunology, food safety, and antimicrobial resistance.

## Key Features

- Functional Mapping: Annotate genes using KEGG, CARD, and VFDB.

- Enrichment Analysis: Perform pathway enrichment from gene lists.

- Pathway-Level Analysis: Visualize gene distribution across biological pathways.

- Resistance & Virulence Profiling: Contextualize resistance and virulence genes functionally.

- Custom Data Input: Supports gene lists and annotation tables in standard formats.

- Extensible Pipeline: Modular design for easy integration into broader omics workflows.

- Visualization Tools: Built-in plots for resistance profiles, virulence categories, enrichment results, and KEGG pathways.

## Examples
You can find full working examples in the `Exemples/` folder:

#### `Exemple 01/`:

Comparative Analysis of Transcriptomic Response of Escherichia coli K-12 MG1655 to Nine Representative Classes of Antibiotics

doi: 10.1128/spectrum.00317-23

- The exemple was done just for the IPM antibiotic

`pepare_data_exemplo01.ipynb`: Preparing data for analysis.

`Analysis_exemple01.ipynb`: ResPathExplorer aplication.

## Structure
```text
📁 ResPathExplorer/
├── 📁 src/ResPathExplorer/
│   ├── __init__.py
│   ├── CARDAnalysis.py
│   ├── KeggAnalysis.py
│   ├── URL_pathway.py
│   ├── VFDBAnalysis.py
│   ├── mapper_KeggFunctions.py
│   ├── rename_file.py
│   ├── save_df_as_html.py
│   └── validate_color_code.py
├── 📁 Exemples/
│   ├── 📁 Exemple 01/
│   └── 📁 Exemple 02/
├── 📁 tests/
```

#### `Exemple 02/`:

Establishment and characterization of persistent Pseudomonas aeruginosa infections in air–liquid interface cultures of human airway epithelial cells.

doi: https://doi.org/10.1128/iai.00603-24

- The exemple was done for the Calu-3 PAO1 Day 5 vs Inoculum

`pepare_data_exemplo02.ipynb`: Preparing data for analysis.

`Analysis_exemple02.ipynb`: ResPathExplorer aplication.

## 🛠 Installation

To install the library directly from GitHub:

```bash
pip install git+https://github.com/lais-carvalho/ResPathExplorer.git
```

## Acknowledgements
- European Food Safety Authority (EFSA) – support via the “Pathogens-in-Foods Database” project.

- Centro de Investigação da Montanha (CIMO), Portugal.

- University of Minho – MSc in Bioinformatics program.

## Contact
For questions or collaborations, open a GitHub Issue or contact: laiscarvalho@ipb.pt
                                                                 linkedin.com/in/laiscristinecarvalho



