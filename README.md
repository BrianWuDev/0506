# Multi-Tumor Network Visualization

A Python tool for creating an interactive gene network visualization with GCH1 as the central node, connected to multiple tumor types and related genes.

## Features

- **Interactive Visualization**: Network graph with GCH1 as the central fixed node, tumor nodes fixed in position, and gene nodes freely draggable
- **Enhanced Visual Design**: GCH1 node prominently positioned in center with larger size, each tumor type displaying a unique vibrant color
- **Gene Correlation**: Visualize correlations between genes and tumor types based on PCC values
- **Cross-Tumor Gene Highlighting**: Special diamond-shaped nodes for genes appearing in multiple tumor types
- **Visual Encoding**: Node size and edge width scaled by correlation strength with enhanced visibility
- **Interactive Controls**: Clustering, expansion, position freezing, and layout reset functionality
- **High-Quality Export**: Options for regular and high-resolution PNG downloads
- **Responsive Design**: Hover tooltips for detailed information and fully interactive interface

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/multi-tumor-network.git
   cd multi-tumor-network
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
Gene Symbol,PCC,Other_column1,Other_column2
GATA3,0.87,value1,value2
FOXA1,0.76,value1,value2
...
```

## Usage

1. Place your tumor data CSV files in the `data` directory.

2. Run the script:
   ```bash
   python multi_tumor_network.py
   ```

3. The script will:
   - Automatically discover available tumor types from CSV files
   - Generate the interactive network visualization
   - Open the visualization in your default web browser
   - Save the HTML file to the `output` directory

## Visual Design

The visualization features:

- **Central GCH1 Node**: Prominently displayed with size 45 and bright red color
- **Tumor Nodes**: Each tumor type has a unique vibrant color with size 30
- **Gene Nodes**: Sized between 8-16 based on correlation strength
- **Cross-Tumor Genes**: Displayed as diamond-shaped nodes with larger size and distinct color
- **Node Labels**: Enhanced font sizes for better readability
- **Edge Thickness**: Scaled based on correlation strength (0.5-2.5)

## Interaction Controls

The interactive visualization includes:

- **Drag Nodes**: All gene nodes can be freely moved with the mouse while GCH1 and tumor nodes remain fixed
- **Zoom**: Use the scroll wheel to zoom in/out
- **Cluster Nodes**: Button to tighten the network layout
- **Expand Nodes**: Button to spread the network layout
- **Reset Layout**: Button to return to the original layout
- **Freeze/Unfreeze Positions**: Buttons to fix node positions or allow movement
- **Download PNG**: Regular and high-resolution PNG download options

## Requirements

- Python 3.6+
- pandas
- numpy
- pyvis
- A modern web browser for the interactive visualization

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- This tool uses [PyVis](https://pyvis.readthedocs.io/) for interactive network visualization
- [html2canvas](https://html2canvas.hertzen.com/) for high-quality PNG exports

## Citation

If you use this tool in your research, please cite:

```
YourName et al. (2024). Multi-Tumor Gene Network Visualization: A tool for 
visualizing gene correlations across multiple tumor types.
GitHub: https://github.com/username/multi-tumor-network
``` 