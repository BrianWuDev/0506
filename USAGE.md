# Usage Guide - Multi-Tumor Network Visualization Tool

This document provides detailed instructions on how to use the Multi-Tumor Network Visualization Tool.

## Directory Structure

```
multi-tumor-network/
│
├── multi_tumor_network.py  # Main program
├── data/                   # Data folder
│   ├── ACC Tumor.csv       # Tumor 1 data
│   ├── BRCA Tumor.csv      # Tumor 2 data
│   └── ...                 # Other tumor data
│
├── output/                 # Output folder
│   └── multi_tumor_network.html  # Generated network graph
│
├── README.md               # Project description
├── USAGE.md                # Usage guide (this file)
└── requirements.txt        # Dependencies
```

## Running the Program

1. Ensure all necessary dependencies are installed:
   ```bash
   pip install -r requirements.txt
   ```

2. Prepare your data files:
   - Save each tumor type's data as a CSV file
   - Name format should be "Tumor Type Tumor.csv", e.g., "BRCA Tumor.csv"
   - Make sure each CSV contains at least the "Gene Symbol" and "PCC" columns

3. Run the program:
   ```bash
   python multi_tumor_network.py
   ```

4. The program will automatically generate the network graph and open it in your browser

## Interactive Operation Guide

### Basic Operations

- **Node Dragging**: All gene nodes can be freely dragged to adjust their positions (GCH1 and tumor nodes remain fixed)
- **Zooming**: Use the mouse wheel to zoom in/out
- **Panning**: Click on empty space and drag to pan the entire view
- **View Details**: Hover over nodes to see detailed information

### Control Button Functions

The bottom-left corner of the view has five control buttons:

1. **Reposition Nodes**: Reorganize nodes, clustering gene nodes more tightly
2. **Expand Layout**: Expand the layout, spreading gene nodes apart
3. **Reset & Stabilize**: Reset the layout and stabilize, returning to the initial state
4. **Fix All Positions**: Fix the positions of all nodes
5. **Enable Dragging**: Enable node dragging functionality

### Download Options

The top-right corner of the view has two download buttons:

1. **Download PNG**: Download standard resolution PNG image
2. **Download Hi-Res PNG**: Download high-resolution PNG image (suitable for publication)

## Node Type Explanation

The network graph contains three types of nodes:

1. **GCH1 Node**: Large red central node, representing the core gene GCH1
2. **Tumor Nodes**: Medium-sized colored nodes, each tumor type has a unique color
3. **Gene Nodes**:
   - **Tumor-specific Genes**: Genes related to only one tumor type, inheriting the tumor's color
   - **Cross-tumor Genes**: Diamond-shaped nodes, genes related to multiple tumor types

## Visual Encoding Explanation

- **Node Size**: Reflects the strength of correlation between genes and GCH1 (PCC value)
- **Edge Width**: Reflects connection strength, width proportional to PCC value
- **Node Color**: Each tumor type has a unique color
- **Node Shape**: Circular for regular genes, diamond for cross-tumor genes

## Custom Configuration

To modify visual parameters, edit the following sections in the `multi_tumor_network.py` file:

- `self.tumor_colors`: Modify tumor node colors
- `self.node_size_range`: Modify gene node size range
- `self.edge_width_range`: Modify edge width range
- `self.central_node_size`: Modify GCH1 node size
- `self.tumor_node_size`: Modify tumor node size
- `self.min_correlation`: Modify minimum correlation threshold

## Frequently Asked Questions

**Q: Why are some tumor nodes not displayed?**  
A: This could be because there is no corresponding data file, or the data file has no genes exceeding the threshold.

**Q: How do I add a new tumor type?**  
A: Create a new CSV file named "Tumor Type Tumor.csv" and ensure it contains the necessary columns.

**Q: How do I modify the position of the GCH1 node?**  
A: The GCH1 node is fixed at the center position. To modify this, edit the relevant sections in the code.

**Q: Why are some gene nodes particularly large?**  
A: Node size reflects the strength of correlation with GCH1; higher PCC values result in larger nodes.

**Q: How do I adjust the minimum correlation threshold?**  
A: Modify the `self.min_correlation` parameter (default is 0.5). 