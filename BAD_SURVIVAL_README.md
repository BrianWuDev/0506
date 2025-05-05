# Bad Survival Tree Visualization

A Python tool for creating a hierarchical survival tree visualization for poor prognosis analysis, using gene expression data from multiple tumor types.

## Features

- **Hierarchical Tree Visualization**: Tree-like structure with "Poor Prognosis" as the central node, tumor types as first-level branches, and genes as leaf nodes
- **Risk-Based Coloring**: Genes are color-coded based on risk level (high, moderate, low) derived from PCC correlation values
- **Visual Encoding**: Node size and edge width scaled by correlation strength for enhanced information density
- **Cross-Tumor Gene Highlighting**: Special diamond-shaped nodes for genes appearing in multiple tumor types
- **Interactive Controls**: Cluster, expand, reset, and position freezing functionality
- **High-Quality Export**: Options for regular and high-resolution PNG downloads
- **Detailed Tooltips**: Hover information showing gene symbols, correlation values, risk levels, and tumor types

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/bad-survival-tree.git
   cd bad-survival-tree
   ```

2. Create a virtual environment (optional but recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Data Format

The tool expects CSV files in the `data` directory with the following format:

- Each file should be named after a tumor type (e.g., `BRCA Tumor.csv`)
- Each CSV should have at least the following columns:
  - `Gene Symbol`: The name/ID of the gene
  - `PCC`: Pearson Correlation Coefficient with the GCH1 gene

Example CSV structure:
```
Gene Symbol,Gene ID,PCC
WNT4,ENSG00000162552.14,0.75
LINC01451,ENSG00000279141.3,0.74
...
```

## Usage

1. Place your tumor data CSV files in the `data` directory.

2. Run the script:
   ```bash
   python bad_survival_tree.py
   ```

3. The script will:
   - Automatically discover available tumor types from CSV files
   - Generate the interactive survival tree visualization
   - Open the visualization in your default web browser
   - Save the HTML file to the `output` directory

## Visualization Design

The visualization features:

- **Central Node**: "Poor Prognosis" node displayed prominently at the center with size 45 and bright red color
- **Tumor Nodes**: Tumor type nodes arranged in a circle around the central node, each with a unique color and size 30
- **Gene Nodes**: 
  - High Risk Genes (PCC ≥ 0.75): Dark red color
  - Moderate Risk Genes (PCC ≥ 0.65): Orange-red color
  - Low Risk Genes (PCC ≥ 0.5): Orange color
- **Cross-Tumor Genes**: Diamond-shaped nodes with 1.5x larger size
- **Node Size**: Scaled between 8-16 based on correlation strength
- **Edge Width**: Scaled between 0.5-2.5 based on correlation strength

## Interaction Controls

The interactive visualization includes:

- **Draggable Nodes**: All gene nodes can be dragged to adjust their positions
- **Cluster Nodes**: Button to tighten the tree layout
- **Expand Nodes**: Button to spread the tree layout
- **Reset & Stabilize**: Button to return to the original layout
- **Freeze/Unfreeze Positions**: Buttons to fix node positions or allow movement
- **Download PNG**: Regular and high-resolution PNG download options

## Requirements

- Python 3.6+
- pandas
- numpy
- pyvis
- A modern web browser for the interactive visualization

## Scientific Background

This visualization is designed to help researchers identify genes that contribute to poor prognosis across multiple tumor types. The hierarchical structure makes it easy to:

1. Identify high-risk genes specific to certain tumor types
2. Discover cross-tumor genes that may indicate common mechanisms of poor prognosis
3. Compare the relative importance of different genes based on their correlation values
4. Explore potential therapeutic targets for improving survival outcomes

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- This tool uses [PyVis](https://pyvis.readthedocs.io/) for interactive network visualization
- [html2canvas](https://html2canvas.hertzen.com/) for high-quality PNG exports

## Citation

If you use this tool in your research, please cite:

```
Author et al. (2024). Bad Survival Tree Visualization: A tool for visualizing 
hierarchical gene correlations related to poor prognosis across multiple tumor types.
GitHub: https://github.com/username/bad-survival-tree
``` 