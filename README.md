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
You can find full working examples in the `Examples/` folder:

#### `Foodborne bacteria/`:

Transcriptome Analysis of *Listeria monocytogenes* Exposed to Beef Fat Reveals Antimicrobial and Pathogenicity Attenuation Mechanisms

doi: https://doi.org/10.1128/AEM.03027-20

- The example was done just for C18:2n-6

`pepare_data.ipynb`: Preparing data for analysis.

`Analysis.ipynb`: ResPathExplorer aplication.

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
├── 📁 Examples/
│   ├── 📁 Foodborne bacteria/
├── 📁 tests/
```

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
                                                                 laismagalhaescarvalho@hotmail.com
                                                                 linkedin.com/in/laiscristinecarvalho




