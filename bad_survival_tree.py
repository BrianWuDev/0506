#!/usr/bin/env python
"""
Survival Tree Visualization for Poor Prognosis

This script creates a hierarchical tree visualization for poor prognosis survival analysis,
using gene expression data from multiple tumor types.

Features:
- Interactive tree visualization with survival-related genes
- Color-coded nodes based on risk level and gene expression
- Node size indicating impact on survival outcome
- Interactive controls for exploring the tree structure
- Detailed tooltips showing survival statistics and p-values
- High-quality PNG export options
"""

import pandas as pd
import numpy as np
import logging
import os
import re
from pathlib import Path
import json
import time
from typing import Dict, List, Tuple, Optional, Union, Any
from pyvis.network import Network
import glob

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('Bad_Survival_Tree')


class BadSurvivalTreeVisualizer:
    """Class for creating a survival tree visualization for poor prognosis.
    
    This class provides functionality to:
    - Create a hierarchical tree visualization for survival analysis
    - Highlight genes associated with poor prognosis
    - Visualize the impact of gene expression on survival outcomes
    - Provide interactive exploration of the tree structure
    - Generate downloadable PNG images of the visualization
    
    The visualization includes detailed tooltips and a legend explaining
    the node types, colors, and connection meanings.
    """
    
    def __init__(self):
        """Initialize the visualizer with optimized settings."""
        # Settings
        self.tumor_types = []  # Will be populated from data files
        self.min_correlation = 0.5
        self.max_genes_per_tumor = 1000
        
        # Node colors with a specialized palette for survival visualization
        self.central_node_color = '#E31A1C'  # Bright red for poor prognosis center
        self.high_risk_color = '#8B0000'     # Dark red for high risk
        self.moderate_risk_color = '#FF4500' # Orange red for moderate risk
        self.low_risk_color = '#FFA500'      # Orange for low risk
        
        # Tumor-specific colors (maintain unique colors for each tumor type)
        self.tumor_colors = {
            'ACC Tumor': '#1F77B4',   # Deep blue
            'BRCA Tumor': '#FF7F0E',  # Orange
            'ESCA Tumor': '#2CA02C',  # Dark green
            'GBM Tumor': '#D62728',   # Red
            'KICH Tumor': '#9467BD',  # Purple
            'LGG Tumor': '#8C564B',   # Brown
            'PCPG Tumor': '#E377C2',  # Pink
            'TGCT Tumor': '#FFC125'   # Bright yellow
        }
        
        # Visual parameters
        self.node_size_range = (8, 16)       # Gene node size range
        self.edge_width_range = (0.5, 2.5)   # Edge width range
        self.central_node_size = 45          # Central node size
        self.tumor_node_size = 30            # Tumor node size
        
        # Output settings
        self.output_dir = 'output'
        self.output_filename = 'bad_survival_tree'
        self.title = 'Survival Tree for Poor Prognosis'
        
        # Create output directory
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize statistics tracking
        self.gene_counts = {}
        
        # Physics optimization parameters for tree-like structure
        self.physics_params = {
            "gravity": -2500,           # Stronger gravity for tree-like structure
            "central_gravity": 0.5,     # Moderate central gravity
            "spring_length": 150,       # Longer springs for tree branching
            "spring_strength": 0.1,     # Moderate spring strength
            "damping": 0.3              # Lower damping for better tree layout
        }

    def discover_tumor_types(self, data_dir: str = 'data') -> List[str]:
        """Discover available tumor types from data directory.
        
        Args:
            data_dir: Directory containing tumor data files
            
        Returns:
            List of discovered tumor types
        """
        tumor_types = []
        data_files = glob.glob(os.path.join(data_dir, "*.csv"))
        
        for file_path in data_files:
            tumor_name = os.path.basename(file_path).replace('.csv', '')
            tumor_types.append(tumor_name)
            logger.info(f"Discovered tumor type: {tumor_name}")
            
        return tumor_types

    def load_all_data(self, data_dir: str = 'data') -> Dict[str, pd.DataFrame]:
        """Load data for all tumor types.
        
        Args:
            data_dir: Directory containing tumor data files
            
        Returns:
            Dictionary of tumor types to DataFrames
        """
        self.tumor_types = self.discover_tumor_types(data_dir)
        all_data = {}
        
        for tumor_type in self.tumor_types:
            file_path = os.path.join(data_dir, f"{tumor_type}.csv")
            try:
                df = pd.read_csv(file_path)
                
                # Initialize gene count for this tumor type
                self.gene_counts[tumor_type] = 0
                
                # Apply correlation threshold
                high_corr_df = df[df['PCC'] >= self.min_correlation]
                logger.info(f"{tumor_type}: Applied correlation threshold ≥ {self.min_correlation}, resulting in {len(high_corr_df)} rows")
                
                all_data[tumor_type] = high_corr_df
                
            except Exception as e:
                logger.error(f"Error loading data for {tumor_type}: {e}")
        
        return all_data

    def create_survival_tree(self, tumor_data: Dict[str, pd.DataFrame]) -> Network:
        """Create hierarchical tree visualization for poor prognosis analysis.
        
        Args:
            tumor_data: Dictionary of tumor types to DataFrames
            
        Returns:
            PyVis Network object
        """
        # Create network with appropriate configuration
        net = Network(
            height='800px', 
            width='1000px', 
            bgcolor='#ffffff', 
            font_color='black', 
            directed=True,  # Use directed edges for tree structure
            notebook=False,
            heading='',
            cdn_resources='remote'
        )
        
        # Configure network display options
        net.toggle_hide_edges_on_drag(False)
        net.toggle_physics(True)
        net.toggle_drag_nodes(True)
        net.toggle_stabilization(True)
        net.barnes_hut(
            gravity=self.physics_params["gravity"],
            central_gravity=self.physics_params["central_gravity"],
            spring_length=self.physics_params["spring_length"],
            spring_strength=self.physics_params["spring_strength"],
            damping=self.physics_params["damping"]
        )

        # Add central "Poor Prognosis" node
        net.add_node(
            "Poor Prognosis", 
            label="Poor Prognosis", 
            color=self.central_node_color, 
            size=self.central_node_size,
            font={"size": 18, "bold": True},
            shape="ellipse",
            fixed=True,  # Fix the position of the central node
            physics=False,  # Disable physics for this node
            x=0,  # Center of the layout
            y=0   # Center of the layout
        )
        
        # Track genes that appear in multiple tumors
        gene_tumors: Dict[str, List[str]] = {}
        
        # Add tumor nodes as first level branches
        for idx, (tumor_type, df) in enumerate(tumor_data.items()):
            if df.empty:
                logger.warning(f"No data for {tumor_type} after filtering")
                continue
            
            # Position tumor nodes in a circle around the central node
            angle = 2 * np.pi * idx / len(tumor_data)
            radius = 250  # Distance from center
            x_pos = radius * np.cos(angle)
            y_pos = radius * np.sin(angle)
            
            # Add tumor node
            net.add_node(
                tumor_type,
                label=tumor_type.replace(' Tumor', ''),
                color=self.tumor_colors.get(tumor_type, '#AAAAAA'),
                size=self.tumor_node_size,
                font={"size": 14, "bold": True},
                shape="box",
                fixed=True,  # Fix the position of tumor nodes
                physics=False,  # Disable physics for tumor nodes
                x=x_pos,
                y=y_pos
            )
            
            # Add edge from central node to tumor
            net.add_edge(
                "Poor Prognosis",
                tumor_type,
                width=3,
                color=self.high_risk_color
            )
            
            # Process genes for this tumor (limit to top genes by correlation)
            top_genes = df.nlargest(self.max_genes_per_tumor, 'PCC')
            
            # Track the gene count
            self.gene_counts[tumor_type] = len(top_genes)
            
            # Add gene nodes
            for _, gene_row in top_genes.iterrows():
                gene_symbol = gene_row['Gene Symbol']
                correlation = gene_row['PCC']
                
                # Track which tumors this gene appears in
                if gene_symbol not in gene_tumors:
                    gene_tumors[gene_symbol] = []
                gene_tumors[gene_symbol].append(tumor_type)
                
                # Unique node ID for this gene in this tumor
                gene_id = f"{gene_symbol}_{tumor_type}"
                
                # Determine risk level based on correlation
                if correlation >= 0.75:
                    risk_color = self.high_risk_color
                    risk_level = "High Risk"
                elif correlation >= 0.65:
                    risk_color = self.moderate_risk_color
                    risk_level = "Moderate Risk"
                else:
                    risk_color = self.low_risk_color
                    risk_level = "Low Risk"
                
                # Calculate node size based on correlation
                size_scale = (correlation - self.min_correlation) / (1 - self.min_correlation)
                node_size = self.node_size_range[0] + size_scale * (self.node_size_range[1] - self.node_size_range[0])
                
                # Calculate edge width based on correlation
                edge_width = self.edge_width_range[0] + size_scale * (self.edge_width_range[1] - self.edge_width_range[0])
                
                # Create tooltip content
                tooltip = f"<b>{gene_symbol}</b><br>Correlation: {correlation:.3f}<br>Risk Level: {risk_level}<br>Tumor: {tumor_type}"
                
                # Add gene node
                net.add_node(
                    gene_id,
                    label=gene_symbol,
                    color=risk_color,
                    size=node_size,
                    title=tooltip,
                    font={"size": 10},
                    shape="dot"
                )
                
                # Add edge from tumor to gene
                net.add_edge(
                    tumor_type,
                    gene_id,
                    width=edge_width,
                    color={"color": risk_color, "opacity": 0.8}
                )
        
        # Identify cross-tumor genes and update their appearance
        for gene, tumors in gene_tumors.items():
            if len(tumors) > 1:
                # This is a cross-tumor gene (appears in multiple tumors)
                logger.info(f"Cross-tumor gene: {gene} appears in {len(tumors)} tumors")
                
                # For each instance of this gene across different tumors
                for tumor_type in tumors:
                    gene_id = f"{gene}_{tumor_type}"
                    
                    # Update this gene node to use a diamond shape and increase size
                    net.get_node(gene_id)["shape"] = "diamond"
                    net.get_node(gene_id)["size"] *= 1.5
                    
                    # Update the tooltip to note this is a cross-tumor gene
                    current_tooltip = net.get_node(gene_id)["title"]
                    updated_tooltip = f"{current_tooltip}<br><b>Cross-tumor gene</b> (in {len(tumors)} tumors)"
                    net.get_node(gene_id)["title"] = updated_tooltip
        
        return net

    def create_html_with_title_and_legend(self, html_path: str, net: Network) -> str:
        """Create HTML with title, controls, legend, and the network visualization.
        
        Args:
            html_path: Path to save the HTML file
            net: PyVis Network object
            
        Returns:
            Path to the created HTML file
        """
        # Generate the basic network HTML
        net.save_graph(html_path)
        
        # Read the generated HTML
        with open(html_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
        
        # CSS for styling the page and controls
        css = """
        <style>
            body {
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 0;
                background-color: #f8f9fa;
            }
            #mynetwork {
                width: 100%;
                height: 800px;
                border: 1px solid #ddd;
                background-color: #ffffff;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
            }
            h1 {
                color: #333;
                text-align: center;
                margin-bottom: 20px;
            }
            .controls {
                position: absolute;
                bottom: 20px;
                left: 20px;
                z-index: 999;
                background: rgba(255, 255, 255, 0.8);
                padding: 10px;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            }
            .download-buttons {
                position: absolute;
                top: 20px;
                right: 20px;
                z-index: 999;
                background: rgba(255, 255, 255, 0.8);
                padding: 10px;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            }
            .legend {
                position: absolute;
                top: 20px;
                left: 20px;
                z-index: 999;
                background: rgba(255, 255, 255, 0.8);
                padding: 10px;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
                max-width: 280px;
            }
            .legend-item {
                display: flex;
                align-items: center;
                margin-bottom: 8px;
            }
            .legend-color {
                width: 20px;
                height: 20px;
                margin-right: 8px;
                border-radius: 3px;
            }
            .legend-shape {
                width: 20px;
                height: 20px;
                margin-right: 8px;
                text-align: center;
                line-height: 20px;
            }
            .legend-shape.diamond {
                transform: rotate(45deg);
                background-color: #2CA02C;
            }
            .legend-label {
                font-size: 14px;
            }
            button {
                background-color: #4CAF50;
                border: none;
                color: white;
                padding: 8px 16px;
                text-align: center;
                text-decoration: none;
                display: inline-block;
                font-size: 14px;
                margin: 4px 2px;
                cursor: pointer;
                border-radius: 4px;
            }
            button:hover {
                background-color: #45a049;
            }
            .status {
                display: inline-block;
                padding: 5px 10px;
                border-radius: 3px;
                font-size: 14px;
                margin-left: 10px;
                background-color: #e9ecef;
            }
            #screenshot-area {
                position: absolute;
                left: -9999px;
            }
        </style>
        """
        
        # JavaScript for interactive controls and functionality
        javascript = """
        <script src="https://html2canvas.hertzen.com/dist/html2canvas.min.js"></script>
        <script>
            // Variables to track state
            let positionsFixed = false;
            let nodePositions = {};
            
            // Wait for the network to be created
            document.addEventListener("DOMContentLoaded", function() {
                // Reference to the network
                const network = document.getElementById("mynetwork").getElementsByTagName("iframe")[0].contentWindow.network;
                
                // Status indicator
                const statusEl = document.getElementById("status");
                
                // Cluster nodes around their center
                document.getElementById("cluster-button").addEventListener("click", function() {
                    statusEl.textContent = "Clustering...";
                    statusEl.style.backgroundColor = "#ffc107";
                    
                    network.setOptions({
                        physics: {
                            enabled: true,
                            barnesHut: {
                                gravitationalConstant: -8000,
                                centralGravity: 0.6,
                                springLength: 150,
                                springConstant: 0.05,
                                damping: 0.5
                            }
                        }
                    });
                    
                    setTimeout(function() {
                        network.stopSimulation();
                        statusEl.textContent = "Ready";
                        statusEl.style.backgroundColor = "#28a745";
                    }, 3000);
                });
                
                // Expand the layout
                document.getElementById("expand-button").addEventListener("click", function() {
                    statusEl.textContent = "Expanding...";
                    statusEl.style.backgroundColor = "#ffc107";
                    
                    network.setOptions({
                        physics: {
                            enabled: true,
                            barnesHut: {
                                gravitationalConstant: -3000,
                                centralGravity: 0.05,
                                springLength: 250,
                                springConstant: 0.01,
                                damping: 0.09
                            }
                        }
                    });
                    
                    setTimeout(function() {
                        network.stopSimulation();
                        statusEl.textContent = "Ready";
                        statusEl.style.backgroundColor = "#28a745";
                    }, 3000);
                });
                
                // Reset and stabilize
                document.getElementById("reset-button").addEventListener("click", function() {
                    statusEl.textContent = "Resetting...";
                    statusEl.style.backgroundColor = "#ffc107";
                    
                    network.setOptions({
                        physics: {
                            enabled: true,
                            barnesHut: {
                                gravitationalConstant: -2500,
                                centralGravity: 0.5,
                                springLength: 150,
                                springConstant: 0.1,
                                damping: 0.3
                            }
                        }
                    });
                    
                    // Fix the central node at the center
                    const positions = network.getPositions(["Poor Prognosis"]);
                    network.moveNode("Poor Prognosis", 0, 0);
                    
                    // Restart physics simulation
                    network.startSimulation();
                    
                    setTimeout(function() {
                        network.stopSimulation();
                        statusEl.textContent = "Ready";
                        statusEl.style.backgroundColor = "#28a745";
                    }, 3000);
                });
                
                // Fix all positions
                document.getElementById("fix-positions-button").addEventListener("click", function() {
                    if (!positionsFixed) {
                        // Save current positions
                        nodePositions = network.getPositions();
                        
                        // Fix all positions
                        Object.keys(nodePositions).forEach(nodeId => {
                            network.moveNode(nodeId, nodePositions[nodeId].x, nodePositions[nodeId].y);
                        });
                        
                        network.setOptions({
                            physics: {
                                enabled: false
                            }
                        });
                        
                        positionsFixed = true;
                        this.textContent = "Unfreeze Positions";
                        statusEl.textContent = "Positions Fixed";
                        statusEl.style.backgroundColor = "#17a2b8";
                    } else {
                        network.setOptions({
                            physics: {
                                enabled: true
                            }
                        });
                        
                        positionsFixed = false;
                        this.textContent = "Freeze Positions";
                        statusEl.textContent = "Positions Unfrozen";
                        statusEl.style.backgroundColor = "#28a745";
                    }
                });
                
                // Enable dragging
                document.getElementById("enable-drag-button").addEventListener("click", function() {
                    network.setOptions({
                        interaction: {
                            dragNodes: true
                        }
                    });
                    
                    statusEl.textContent = "Dragging Enabled";
                    statusEl.style.backgroundColor = "#28a745";
                });
                
                // Download as PNG
                document.getElementById("download-png-button").addEventListener("click", function() {
                    const iframe = document.getElementById("mynetwork").getElementsByTagName("iframe")[0];
                    const iframeContent = iframe.contentDocument || iframe.contentWindow.document;
                    
                    html2canvas(iframeContent.body).then(canvas => {
                        const link = document.createElement('a');
                        link.download = 'bad_survival_tree.png';
                        link.href = canvas.toDataURL("image/png");
                        link.click();
                    });
                });
                
                // High-resolution PNG
                document.getElementById("download-hires-button").addEventListener("click", function() {
                    const iframe = document.getElementById("mynetwork").getElementsByTagName("iframe")[0];
                    const iframeContent = iframe.contentDocument || iframe.contentWindow.document;
                    
                    html2canvas(iframeContent.body, {
                        scale: 3, // Higher resolution
                        backgroundColor: "#ffffff"
                    }).then(canvas => {
                        const link = document.createElement('a');
                        link.download = 'bad_survival_tree_hires.png';
                        link.href = canvas.toDataURL("image/png");
                        link.click();
                    });
                });
                
                // Fix central node and tumor nodes positions
                setTimeout(function() {
                    // Make sure the central node is at the center of the layout
                    network.moveNode("Poor Prognosis", 0, 0);
                    
                    // Stop physics simulation
                    network.stopSimulation();
                    
                    statusEl.textContent = "Ready";
                    statusEl.style.backgroundColor = "#28a745";
                }, 2000);
            });
        </script>
        """
        
        # HTML content for title, controls, and legend
        html_components = f"""
        <div class="container">
            <h1>{self.title}</h1>
            
            <div class="legend">
                <h3>Legend</h3>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: {self.central_node_color};"></div>
                    <div class="legend-label">Poor Prognosis (Central Node)</div>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: {self.high_risk_color};"></div>
                    <div class="legend-label">High Risk Gene (PCC ≥ 0.75)</div>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: {self.moderate_risk_color};"></div>
                    <div class="legend-label">Moderate Risk Gene (PCC ≥ 0.65)</div>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: {self.low_risk_color};"></div>
                    <div class="legend-label">Low Risk Gene (PCC ≥ 0.5)</div>
                </div>
                <div class="legend-item">
                    <div class="legend-shape diamond"></div>
                    <div class="legend-label">Cross-tumor Gene</div>
                </div>
                <p><small>Node size indicates correlation strength. Edge thickness indicates relationship strength.</small></p>
            </div>
            
            <div class="controls">
                <button id="cluster-button">Cluster Nodes</button>
                <button id="expand-button">Expand Layout</button>
                <button id="reset-button">Reset & Stabilize</button>
                <button id="fix-positions-button">Freeze Positions</button>
                <button id="enable-drag-button">Enable Dragging</button>
                <span id="status" class="status" style="background-color: #ffc107;">Initializing...</span>
            </div>
            
            <div class="download-buttons">
                <button id="download-png-button">Download PNG</button>
                <button id="download-hires-button">Download Hi-Res PNG</button>
            </div>
        </div>
        """
        
        # Combine all components
        modified_html = html_content.replace('<head>', f'<head>{css}')
        modified_html = modified_html.replace('</head>', f'{javascript}</head>')
        modified_html = modified_html.replace('<body>', f'<body>{html_components}')
        
        # Write the modified HTML back to file
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(modified_html)
        
        return html_path

    def run(self, data_dir: str = 'data') -> str:
        """Run the entire pipeline to create the survival tree visualization.
        
        Args:
            data_dir: Directory containing tumor data files
            
        Returns:
            Path to the generated HTML file
        """
        logger.info(f"Starting bad survival tree visualization generation")
        
        # Load all tumor data
        tumor_data = self.load_all_data(data_dir)
        
        # Create the network visualization
        net = self.create_survival_tree(tumor_data)
        
        # Generate the HTML file
        output_path = self.output_dir / f"{self.output_filename}.html"
        final_html_path = self.create_html_with_title_and_legend(str(output_path), net)
        
        logger.info(f"Survival tree visualization generated: {final_html_path}")
        logger.info("Tumor gene counts:")
        
        for tumor_type, count in self.gene_counts.items():
            logger.info(f"  - {tumor_type}: {count} genes")
        
        return final_html_path


if __name__ == "__main__":
    visualizer = BadSurvivalTreeVisualizer()
    output_file = visualizer.run()
    
    # Open the visualization in the default web browser
    import webbrowser
    webbrowser.open(f"file://{os.path.abspath(output_file)}") 