#!/usr/bin/env python
"""
Multi-Tumor Network Visualization with GCH1 as Central Node

This script creates a network visualization with GCH1 as the central node,
connected to multiple tumor types, with genes related to each tumor type.

Features:
- Interactive network visualization with GCH1 as central node
- Multiple tumor types displayed as second-level nodes
- Tumor-specific genes connected to their respective tumor types
- Cross-tumor genes (appearing in multiple tumors) shown with special diamond shape
- Node size and edge width scaled by correlation strength
- Interactive controls for clustering, expanding, and resetting layout
- Position freezing functionality to fix node positions
- PNG download options (regular and high-resolution)
- Responsive design with hover tooltips for additional information
"""

import pandas as pd
import numpy as np
import logging
import os
import re
from pathlib import Path
import json
import time
from typing import Dict, List, Tuple, Optional, Union
from pyvis.network import Network
import glob

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('Multi_Tumor_Network')


class MultiTumorNetworkVisualizer:
    """Class for creating multi-tumor gene network with GCH1 as central node.
    
    This class provides functionality to:
    - Automatically discover available tumor types from data files
    - Load gene correlation data for multiple tumor types
    - Create an interactive network visualization with GCH1 as central node
    - Display tumor-specific genes connected to respective tumor nodes
    - Highlight cross-tumor genes (genes appearing in multiple tumor types)
    - Provide interactive controls for manipulating the network layout
    - Allow freezing and unfreezing of node positions for better exploration
    - Generate downloadable PNG images of the visualization
    - Optimize node placement using physics-based algorithms
    
    The visualization includes status indicators, control buttons, and
    a detailed legend explaining the network components.
    """
    
    def __init__(self):
        """Initialize the visualizer with optimized settings."""
        # Settings
        self.tumor_types = []  # Will be populated from data files
        self.min_correlation = 0.5
        self.max_genes_per_tumor = 1000
        
        # Node colors - professional color palette
        self.central_node_color = '#FF4136'  # Red for GCH1
        self.tumor_colors = {
            'ACC Tumor': '#3D9970',    # Green
            'BRCA Tumor': '#0074D9',   # Blue
            'ESCA Tumor': '#FF851B',   # Orange
            'GBM Tumor': '#B10DC9',    # Purple
            'KICH Tumor': '#FF4081',   # Pink
            'LGG Tumor': '#2ECC40',    # Light green
            'PCPG Tumor': '#F012BE',   # Magenta
            'TGCT Tumor': '#01FF70',   # Lime
            'SARC Tumor': '#85144b',   # Maroon
            'OV Tumor': '#FFDC00',     # Yellow
            'BLCA Tumor': '#39CCCC'    # Teal
        }
        
        # Visual parameters - 減小節點大小範圍
        self.node_size_range = (3, 10)  # 原本是 (5, 12)
        self.edge_width_range = (0.3, 1.5)  # 稍微減小邊的寬度
        
        # Output settings
        self.output_dir = 'output'
        self.output_filename = 'multi_tumor_network'
        self.title = 'Multi-Tumor Network with GCH1 as Central Node'
        
        # Create output directory
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize statistics tracking
        self.gene_counts = {}
        
        # Physics optimization parameters
        self.physics_params = {
            "gravity": -2000,       # 減少重力以避免過度聚集
            "central_gravity": 0.05, # 減少中心引力
            "spring_length": 150,    # 增加彈簧長度使節點分散
            "spring_strength": 0.05, # 減少彈簧強度使連接更鬆散
            "damping": 0.6,          # 增加阻尼以減少振動
            "avoid_overlap": 0.8     # 增加節點間間隔
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

    def create_network(self, tumor_data: Dict[str, pd.DataFrame]) -> Network:
        """Create network with nodes and edges for multiple tumor types.
        
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
            directed=False,
            notebook=False,
            heading='',
            cdn_resources='remote'
        )
        
        # Configure network display options - 關閉物理引擎但啟用節點拖曳
        net.toggle_hide_edges_on_drag(False)
        net.toggle_physics(False)  # 關閉物理引擎，防止節點亂動
        net.toggle_drag_nodes(True)  # 明確啟用節點拖曳
        net.toggle_stabilization(True)
        
        # 設置交互選項，確保節點可拖曳
        options_str = """
        {
            "interaction": {
                "dragNodes": true,
                "dragView": true,
                "zoomView": true,
                "selectable": true,
                "hover": true,
                "navigationButtons": true
            },
            "physics": {
                "enabled": false
            }
        }
        """
        net.set_options(options_str)
        
        # Add GCH1 as center node with absolute fixed position
        net.add_node(
            'GCH1', 
            label='GCH1', 
            size=30,  # 增大中心節點大小
            color=self.central_node_color, 
            title='GCH1 (Central Gene)', 
            shape='dot', 
            borderWidth=4,  # 增加邊框粗細
            font={'size': 18, 'face': 'Arial', 'color': 'white', 'strokeWidth': 3, 'strokeColor': '#000000'},
            x=0,
            y=0,
            fixed={  # 強制固定位置的詳細設置
                'x': True,
                'y': True
            },
            physics=False  # 完全禁用物理引擎
        )
        
        # Position tumor nodes radially around GCH1
        num_tumors = len(self.tumor_types)
        radius = 400  # 增加距離中心的半徑，原本是300
        
        # Collect all genes and their occurrences in tumors
        all_genes = {}  # Format: {gene_id: {tumor_type: pcc, ...}, ...}
        
        # Step 1: Collect all genes and their PCC values in different tumors
        for tumor_type, df in tumor_data.items():
            for _, row in df.iterrows():
                gene_id = row['Gene Symbol']
                pcc = row['PCC']
                
                if gene_id not in all_genes:
                    all_genes[gene_id] = {}
                
                all_genes[gene_id][tumor_type] = pcc
        
        # Identify cross-tumor genes (appearing in multiple tumor types)
        cross_tumor_genes = {gene: tumors for gene, tumors in all_genes.items() if len(tumors) > 1}
        logger.info(f"Found {len(cross_tumor_genes)} cross-tumor genes")
        
        # Add tumor nodes
        tumor_positions = {}  # Store tumor node positions
        
        for idx, tumor_type in enumerate(self.tumor_types):
            # Skip if no data for this tumor type
            if tumor_type not in tumor_data:
                continue
                
            # Calculate position on circle
            angle = 2 * np.pi * idx / num_tumors
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            tumor_positions[tumor_type] = (x, y)  # Save position
            
            # Use default color if tumor type not in color map
            color = self.tumor_colors.get(tumor_type, '#666666')
            
            # Add tumor node with absolute fixed position
            net.add_node(
                tumor_type, 
                label=tumor_type, 
                size=20,  # 增大腫瘤節點大小
                color=color,
                title=f"{tumor_type}", 
                shape='dot',
                borderWidth=3,  # 增加邊框粗細
                font={'size': 14, 'face': 'Arial', 'color': 'white', 'strokeWidth': 2, 'strokeColor': '#000000'},
                x=x,
                y=y,
                fixed={  # 強制固定位置的詳細設置
                    'x': True,
                    'y': True
                },
                physics=False  # 完全禁用物理引擎
            )
            
            # Connect tumor node to GCH1
            net.add_edge('GCH1', tumor_type, width=2, color='rgba(150,150,150,0.8)')
            
            # Initialize gene count for this tumor
            self.gene_counts[tumor_type] = 0
        
        # Add tumor-specific genes (appearing in only one tumor type)
        for tumor_type in self.tumor_types:
            if tumor_type not in tumor_data:
                continue
                
            # Get tumor position
            x, y = tumor_positions[tumor_type]
            color = self.tumor_colors.get(tumor_type, '#666666')
            
            # Get tumor-specific genes (only appearing in this tumor)
            tumor_specific_genes = tumor_data[tumor_type]
            tumor_specific_genes = tumor_specific_genes[
                tumor_specific_genes['Gene Symbol'].apply(
                    lambda g: g in all_genes and len(all_genes[g]) == 1
                )
            ]
            
            # Add tumor-specific genes
            spiral_factor = 0.35  # 增加螺旋因子，原本是0.2
            # Check if DataFrame is empty to avoid sorting error
            if not tumor_specific_genes.empty:
                tumor_specific_genes = tumor_specific_genes.sort_values('PCC', ascending=False)
            logger.info(f"{tumor_type}: Adding {len(tumor_specific_genes)} tumor-specific genes")
            
            for gene_idx, (_, gene_row) in enumerate(tumor_specific_genes.iterrows()):
                gene_id = gene_row['Gene Symbol']
                pcc = gene_row['PCC']
                
                # Calculate node size based on correlation
                size_min, size_max = self.node_size_range
                size = size_min + (pcc - self.min_correlation) * (size_max - size_min) / (1 - self.min_correlation)
                
                # Calculate edge width based on correlation
                width_min, width_max = self.edge_width_range
                width = width_min + (pcc - self.min_correlation) * (width_max - width_min)
                
                # Spiral layout - 增加基因間距
                spiral_angle = 2 * np.pi * gene_idx / 15  # 減少每圈的基因數，原本是20
                spiral_distance = 80 + spiral_factor * gene_idx * 3  # 增加起始距離和增長速度
                gene_x = x + spiral_distance * np.cos(spiral_angle)
                gene_y = y + spiral_distance * np.sin(spiral_angle)
                
                # Add tooltip information
                tooltip = f"{gene_id}<br>PCC: {pcc:.3f}<br>Tumor: {tumor_type}"
                
                net.add_node(
                    gene_id,  # Use gene ID as node ID
                    label=gene_id, 
                    title=tooltip,
                    size=size, 
                    color={'background': color, 'border': color},
                    shape='dot',
                    borderWidth=1,
                    font={'size': 6},  # 減小字體，原本是8
                    x=gene_x,
                    y=gene_y
                    # 移除固定位置設置，允許拖曳
                )
                
                # Connect gene to tumor
                net.add_edge(
                    tumor_type, 
                    gene_id, 
                    width=width, 
                    title=f"PCC: {pcc:.3f}",
                    color={'color': color, 'opacity': 0.5}
                )
                
                self.gene_counts[tumor_type] += 1
        
        # Add cross-tumor genes
        logger.info("Adding cross-tumor genes:")
        for gene_id, tumor_pccs in cross_tumor_genes.items():
            # Use gene with highest PCC value in the main tumor
            main_tumor = max(tumor_pccs.items(), key=lambda x: x[1])[0]
            main_pcc = tumor_pccs[main_tumor]
            
            # Only process genes with minimum correlation threshold
            if main_pcc < self.min_correlation:
                continue
                
            # Get main tumor position and color
            main_x, main_y = tumor_positions[main_tumor]
            
            # Calculate position - 改進跨腫瘤基因位置計算
            gene_x = main_x * 0.5  # 靠近中心點更多，原本是0.6
            gene_y = main_y * 0.5  # 靠近中心點更多，原本是0.6
            
            # 添加隨機偏移以減少重疊
            offset_x = np.random.uniform(-40, 40)
            offset_y = np.random.uniform(-40, 40)
            gene_x += offset_x
            gene_y += offset_y
            
            # Calculate node size based on highest correlation
            size_min, size_max = self.node_size_range
            size = size_min + (main_pcc - self.min_correlation) * (size_max - size_min) / (1 - self.min_correlation)
            
            # Multiple colors mixed, create unique node color
            combined_color = '#666666'  # Default color
            tumors_str = ', '.join(tumor_pccs.keys())
            
            tooltip = f"{gene_id}<br>Cross-tumor gene<br>Exists in: {tumors_str}<br>"
            tooltip += "<br>".join([f"{t}: PCC={p:.3f}" for t, p in tumor_pccs.items()])
            
            # Add cross-tumor gene node
            net.add_node(
                gene_id,
                label=gene_id,
                title=tooltip,
                size=size + 2,  # Slightly enlarge cross-tumor gene node
                color={'background': '#FF9800', 'border': '#E65100'},  # Use orange to highlight
                shape='diamond',  # Use diamond shape to distinguish cross-tumor gene
                borderWidth=2,
                font={'size': 7},  # 減小字體，原本是9
                x=gene_x,
                y=gene_y
                # 移除固定位置設置，允許拖曳
            )
            
            # Add connections for each associated tumor type
            for tumor_type, pcc in tumor_pccs.items():
                if pcc >= self.min_correlation:
                    tumor_color = self.tumor_colors.get(tumor_type, '#666666')
                    
                    # Calculate edge width based on correlation
                    width_min, width_max = self.edge_width_range
                    width = width_min + (pcc - self.min_correlation) * (width_max - width_min)
                    
                    net.add_edge(
                        tumor_type,
                        gene_id,
                        width=width,
                        title=f"{tumor_type} - {gene_id}: PCC={pcc:.3f}",
                        color={'color': tumor_color, 'opacity': 0.6}
                    )
                    
                    # Increase count
                    self.gene_counts[tumor_type] += 1
            
            logger.info(f"  {gene_id}: Exists in {len(tumor_pccs)} tumor types ({', '.join(tumor_pccs.keys())})")
        
        # Statistics and logging output
        logger.info("Gene network statistics:")
        print("\nGene network statistics:")
        
        # Tumor-specific gene statistics
        tumor_specific_count = sum(1 for gene, tumors in all_genes.items() if len(tumors) == 1)
        multi_tumor_count = len(cross_tumor_genes)
        
        print(f"Total genes: {len(all_genes)}")
        print(f"Tumor-specific genes: {tumor_specific_count}")
        print(f"Cross-tumor gene: {multi_tumor_count}")
        
        for tumor, count in self.gene_counts.items():
            logger.info(f"  {tumor}: {count} genes")
            print(f"  {tumor}: {count} genes")
        
        # Add initialization script to fix all nodes but allow dragging
        init_script = f"""
        <script type="text/javascript">
        document.addEventListener('DOMContentLoaded', function() {{
            // 在DOM載入後立即固定腫瘤節點和GCH1
            setTimeout(function() {{
                try {{
                    if (window.network) {{
                        console.log("Network object found, fixing tumor nodes and GCH1...");
                        const nodesDataset = network.body.data.nodes;
                        
                        // 只固定腫瘤節點和GCH1的位置
                        nodesDataset.forEach(function(node) {{
                            if (node.id === 'GCH1' || node.id.includes('Tumor')) {{
                                nodesDataset.update({{
                                    id: node.id,
                                    fixed: {{
                                        x: true,
                                        y: true
                                    }},
                                    physics: false
                                }});
                            }}
                        }});
                        
                        // 確保所有基因節點可拖曳
                        nodesDataset.forEach(function(node) {{
                            if (node.id !== 'GCH1' && !node.id.includes('Tumor')) {{
                                nodesDataset.update({{
                                    id: node.id,
                                    fixed: {{
                                        x: false,
                                        y: false
                                    }}
                                }});
                            }}
                        }});
                        
                        // 關閉物理引擎
                        network.setOptions({{ physics: false }});
                        
                        // 確保節點可拖曳
                        network.setOptions({{ 
                            interaction: {{ 
                                dragNodes: true 
                            }} 
                        }});
                        
                        console.log("Tumor nodes and GCH1 fixed. Gene nodes can be dragged freely.");
                        document.getElementById('physics-status').textContent = 'Tumor Nodes Fixed - Drag Gene Nodes Freely';
                        document.getElementById('physics-status').className = 'status-indicator frozen';
                    }}
                }} catch(e) {{
                    console.error("Error fixing nodes:", e);
                }}
            }}, 500);
        }});
        </script>
        """
        
        # Add the script to the HTML
        net.html += init_script
        
        # Simplify stabilization message
        stabilize_script = """
        <script type="text/javascript">
        setTimeout(function() {
            try {
                if (window.network) {
                    // 初始穩定化
                    network.stabilize(100);
                    console.log("Initial stabilization done.");
                }
            } catch(e) {
                console.error("Error in stabilization:", e);
            }
        }, 1000);
        </script>
        """
        
        # Add stabilization script
        net.html += stabilize_script
        
        return net

    def create_html_with_title_and_legend(self, html_path: str, net: Network) -> str:
        """Create HTML with title, legend, and enhanced download buttons.
        
        This method enhances the basic PyVis HTML output with:
        - Custom title and subtitle
        - Interactive legend with color coding for different node types
        - Enhanced control buttons for node clustering, expansion, and layout reset
        - Position freezing/unfreezing functionality to control node movement
        - Status indicator showing physics simulation state
        - Regular and high-resolution PNG download options
        - Responsive layout with visual enhancements
        - Browser compatibility improvements for various environments
        
        Args:
            html_path: Path to save the HTML file
            net: PyVis Network object
            
        Returns:
            Path to the modified HTML file
        """
        # Save basic network HTML
        net.save_graph(html_path)
        
        # Read the HTML content with explicit utf-8 encoding
        with open(html_path, 'r', encoding='utf-8') as file:
            html_content = file.read()
        
        # Remove any horizontal rules or dividers in the original HTML
        html_content = html_content.replace('<hr>', '')
        html_content = html_content.replace('<hr/>', '')
        html_content = html_content.replace('<hr />', '')
        
        # Add our CSS to the head
        css_styles = '''
        <style>
            /* Container styles */
            .network-container {
                max-width: 1200px;
                width: 100%;
                margin: 0 auto;
                padding: 0;
                display: flex;
                flex-direction: column;
                align-items: center;
                border: none !important;
                background-color: #f9f9f9;
            }
            
            /* Title styles */
            .header-container {
                text-align: center;
                margin-bottom: 20px;
                width: 100%;
                border-bottom: none !important;
                padding-bottom: 0 !important;
                box-shadow: none !important;
            }
            
            .main-title {
                color: #333;
                font-size: 24px;
                margin-bottom: 10px;
            }
            
            .subtitle {
                color: #666;
                font-size: 16px;
            }
            
            /* Network styles */
            #mynetwork {
                width: 1000px !important;
                height: 800px !important;
                margin: 0 auto !important;
                border: none !important;
                outline: none !important;
                box-shadow: none !important;
            }
            
            /* Legend styles */
            .legend-container {
                position: absolute;
                top: 20px;
                left: 20px;
                background-color: white;
                padding: 15px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                z-index: 1000;
                max-width: 220px;
                font-size: 14px;
            }
            
            .legend-title {
                text-align: center;
                margin-bottom: 10px;
                font-weight: bold;
            }
            
            .legend-item {
                display: flex;
                align-items: center;
                margin-bottom: 5px;
            }
            
            .legend-color {
                width: 15px;
                height: 15px;
                border-radius: 50%;
                margin-right: 10px;
            }
            
            .legend-section {
                margin-top: 10px;
                font-weight: bold;
            }
            
            .legend-info {
                margin-top: 10px;
                font-size: 12px;
            }
            
            .legend-controls {
                margin-top: 10px;
                border-top: 1px solid #eee;
                padding-top: 10px;
                font-size: 12px;
            }
            
            .controls-title {
                font-weight: bold;
                margin-bottom: 5px;
            }
            
            /* Status indicator */
            .status-indicator {
                position: absolute;
                bottom: 70px;
                left: 20px;
                background-color: rgba(255, 255, 255, 0.9);
                padding: 8px 14px;
                border-radius: 5px;
                font-size: 14px;
                font-weight: bold;
                z-index: 1000;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                border: 1px solid rgba(0,0,0,0.05);
                transition: all 0.3s ease;
            }
            
            .status-indicator.frozen {
                color: #2196F3;
                border-left: 4px solid #2196F3;
            }
            
            .status-indicator.active {
                color: #4CAF50;
                border-left: 4px solid #4CAF50;
            }
            
            .status-indicator.error {
                color: #F44336;
                border-left: 4px solid #F44336;
                background-color: rgba(244, 67, 54, 0.1);
            }
            
            /* Control buttons */
            .control-buttons {
                position: absolute;
                bottom: 20px;
                left: 20px;
                z-index: 1000;
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
            }
            
            .control-buttons button {
                padding: 8px 12px;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                font-weight: bold;
                transition: all 0.2s ease;
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            }
            
            #cluster-btn {
                background-color: #9C27B0;
            }
            
            #expand-btn {
                background-color: #FF9800;
            }
            
            #reset-btn {
                background-color: #607D8B;
            }
            
            #freeze-btn {
                background-color: #2196F3;
            }
            
            #unfreeze-btn {
                background-color: #4CAF50;
            }
            
            .control-buttons button:hover {
                filter: brightness(1.1);
                box-shadow: 0 4px 8px rgba(0,0,0,0.3);
                transform: translateY(-2px);
            }
            
            .control-buttons button:active {
                transform: translateY(1px);
                box-shadow: 0 1px 3px rgba(0,0,0,0.2);
            }
            
            /* Download buttons */
            .download-buttons {
                position: absolute;
                top: 20px;
                right: 20px;
                z-index: 1000;
                display: flex;
            }
            
            #download-btn, #download-hires-btn {
                padding: 8px 12px;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                font-weight: bold;
                transition: all 0.2s ease;
            }
            
            #download-btn {
                background-color: #4CAF50;
                margin-right: 10px;
            }
            
            #download-hires-btn {
                background-color: #2196F3;
            }
            
            #download-btn:hover, #download-hires-btn:hover,
            .control-buttons button:hover {
                filter: brightness(1.1);
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            }
        </style>
        '''
        
        # Create title HTML
        title_html = f'''
        <div class="header-container">
            <h1 class="main-title">{self.title}</h1>
            <p class="subtitle">Visualization of GCH1 gene correlations with multiple tumor types and associated genes</p>
        </div>
        '''
        
        # Create legend HTML
        legend_html = f'''
        <div class="legend-container">
            <div class="legend-title">Legend</div>
            <div class="legend-item">
                <div class="legend-color" style="background-color: {self.central_node_color};"></div>
                <span>GCH1 (Central)</span>
            </div>
            
            <div class="legend-section">Cancer Types:</div>
        '''
        
        # Add tumor types to the legend
        for tumor_type in self.tumor_types:
            color = self.tumor_colors.get(tumor_type, '#666666')
            legend_html += f'''
            <div class="legend-item">
                <div class="legend-color" style="background-color: {color};"></div>
                <span>{tumor_type}</span>
            </div>
            '''
        
        # Add correlation info to legend
        legend_html += f'''
            <div class="legend-info">
                <div>Node size: PCC correlation strength</div>
                <div>Edge width: Connection strength</div>
                <div>Minimum correlation: {self.min_correlation}</div>
                <div>Diamond shape: Cross-tumor genes</div>
            </div>
            
            <div class="legend-controls">
                <div class="controls-title">Interactive Controls:</div>
                <div>- Drag nodes to rearrange</div>
                <div>- Scroll to zoom in/out</div>
                <div>- Hover over nodes for details</div>
                <div>- Use buttons below for layout control</div>
            </div>
        </div>
        '''
        
        # Status indicator
        status_indicator = '''
        <div id="physics-status" class="status-indicator active">
            Physics: Active
        </div>
        '''
        
        # Add control buttons
        control_buttons = '''
        <div class="control-buttons">
            <button id="cluster-btn">Reposition Nodes</button>
            <button id="expand-btn">Expand Layout</button>
            <button id="reset-btn">Reset & Stabilize</button>
            <button id="freeze-btn">Fix All Positions</button>
            <button id="unfreeze-btn">Enable Dragging</button>
        </div>
        '''
        
        # Add JavaScript for control buttons
        control_script = '''
        <script>
        // Get the network instance
        let network = null;
        let nodesDataset = null;
        let edgesDataset = null;
        let isPhysicsEnabled = false;  // 物理引擎預設關閉
        let savedPositions = {};
        
        // Wait for the network to be created
        document.addEventListener('DOMContentLoaded', function() {
            // The network object is stored in the global window object by vis.js
            // Wait a bit for the network to initialize
            setTimeout(function() {
                try {
                    // Try multiple ways to access the network in different browsers
                    if (window.network) {
                        network = window.network;
                        console.log("Network object found via window.network");
                    } else if (document.querySelector('.vis-network') && document.querySelector('.vis-network').__vis_network__) {
                        network = document.querySelector('.vis-network').__vis_network__;
                        console.log("Network object found via DOM element");
                    }
                    
                    if (network) {
                        // Get the nodes and edges datasets
                        try {
                            nodesDataset = network.body.data.nodes;
                            edgesDataset = network.body.data.edges;
                            
                            // Save initial positions immediately for reset functionality
                            saveCurrentNodePositions();
                            
                            // Update status indicator to match initial state
                            updatePhysicsStatus(false);
                            
                            // Initialize control buttons
                            initControlButtons();
                        } catch (e) {
                            console.error("Error accessing network data:", e);
                        }
                    } else {
                        console.warn("Network object not found on initial attempt, will retry");
                        // Retry after a longer delay
                        setTimeout(initializeNetworkAccess, 2000);
                    }
                } catch (e) {
                    console.error("Error during network initialization:", e);
                    // Retry after a longer delay
                    setTimeout(initializeNetworkAccess, 2000);
                }
            }, 1000);
        });
        
        function initializeNetworkAccess() {
            try {
                // Final attempt to get network object
                if (window.network) {
                    network = window.network;
                } else {
                    // Try to find the network through the DOM
                    const visNetworkElements = document.querySelectorAll('.vis-network');
                    for (let i = 0; i < visNetworkElements.length; i++) {
                        if (visNetworkElements[i].__vis_network__) {
                            network = visNetworkElements[i].__vis_network__;
                            break;
                        }
                    }
                }
                
                if (network) {
                    console.log("Network object found on retry");
                    nodesDataset = network.body.data.nodes;
                    edgesDataset = network.body.data.edges;
                    
                    // Save initial positions immediately for reset functionality
                    saveCurrentNodePositions();
                    
                    // Update status indicator to match initial state
                    updatePhysicsStatus(false);
                    
                    // Initialize control buttons
                    initControlButtons();
                } else {
                    console.error("Could not find network object after multiple attempts");
                    // Show error message to user
                    showNetworkError();
                }
            } catch (e) {
                console.error("Error during network retry initialization:", e);
                showNetworkError();
            }
        }
        
        function showNetworkError() {
            const statusEl = document.getElementById('physics-status');
            if (statusEl) {
                statusEl.textContent = 'Error: Network not initialized';
                statusEl.className = 'status-indicator error';
                statusEl.style.color = 'red';
            }
        }
        
        function initControlButtons() {
            // Add event listeners to buttons only after network is available
            document.getElementById('cluster-btn').addEventListener('click', clusterNodes);
            document.getElementById('expand-btn').addEventListener('click', expandNodes);
            document.getElementById('reset-btn').addEventListener('click', resetLayout);
            document.getElementById('freeze-btn').addEventListener('click', freezePositions);
            document.getElementById('unfreeze-btn').addEventListener('click', unfreezePositions);
            
            console.log("Control buttons initialized");
        }
        
        function updatePhysicsStatus(isEnabled) {
            const statusEl = document.getElementById('physics-status');
            if (statusEl) {
                if (isEnabled) {
                    statusEl.textContent = 'Physics: Active, Stabilizing...';
                    statusEl.className = 'status-indicator active';
                } else {
                    statusEl.textContent = 'Physics: Off, Drag Nodes Freely';
                    statusEl.className = 'status-indicator frozen';
                }
            }
        }
        
        function saveCurrentNodePositions() {
            if (!network || !nodesDataset) {
                console.warn("Cannot save positions: network or nodes not available");
                return;
            }
            
            try {
                const positions = network.getPositions();
                savedPositions = positions;
                console.log("Saved current positions", positions);
            } catch (e) {
                console.error("Error saving node positions:", e);
            }
        }
        
        function applyFixedPositions() {
            if (!network || !nodesDataset) {
                console.warn("Cannot apply fixed positions: network or nodes not available");
                return;
            }
            
            try {
                const nodeIds = Object.keys(savedPositions);
                
                nodeIds.forEach(nodeId => {
                    // 跳過 GCH1 和腫瘤節點，它們已經是固定的
                    if (nodeId === 'GCH1' || nodeId.includes("Tumor")) return;
                    
                    const pos = savedPositions[nodeId];
                    if (pos && typeof pos.x === 'number' && typeof pos.y === 'number') {
                        nodesDataset.update({
                            id: nodeId,
                            x: pos.x,
                            y: pos.y,
                            fixed: {
                                x: true,
                                y: true
                            }
                        });
                    }
                });
                
                console.log("Applied fixed positions to gene nodes only");
            } catch (e) {
                console.error("Error applying fixed positions:", e);
            }
        }
        
        function releaseFixedPositions() {
            if (!network || !nodesDataset) {
                console.warn("Cannot release fixed positions: network or nodes not available");
                return;
            }
            
            try {
                const ids = nodesDataset.getIds();
                
                ids.forEach(id => {
                    // 跳過 GCH1 和腫瘤節點，保持它們固定
                    if (id === 'GCH1' || id.includes("Tumor")) return;
                    
                    nodesDataset.update({
                        id: id,
                        fixed: {
                            x: false,
                            y: false
                        }
                    });
                });
                
                console.log("Released fixed positions on gene nodes only");
            } catch (e) {
                console.error("Error releasing fixed positions:", e);
            }
        }
        
        // Cluster nodes
        function clusterNodes() {
            if (!network) {
                console.warn("Cannot cluster nodes: network not available");
                return;
            }
            
            try {
                // Increase gravity to pull nodes together
                network.physics.options.forceAtlas2Based.gravitationalConstant = -500;
                network.physics.options.forceAtlas2Based.centralGravity = 0.4;
                network.physics.options.forceAtlas2Based.springLength = 50;
                network.physics.options.forceAtlas2Based.springConstant = 0.3;
                
                // Update physics and stabilize
                network.setOptions({physics: network.physics.options});
                
                // Make sure physics is enabled
                if (!isPhysicsEnabled) {
                    network.setOptions({ physics: true });
                    isPhysicsEnabled = true;
                    updatePhysicsStatus(true);
                    releaseFixedPositions();
                }
                
                network.stabilize(100);
            } catch (e) {
                console.error("Error clustering nodes:", e);
            }
        }
        
        // Expand nodes
        function expandNodes() {
            if (!network) {
                console.warn("Cannot expand nodes: network not available");
                return;
            }
            
            try {
                // Decrease gravity to push nodes apart
                network.physics.options.forceAtlas2Based.gravitationalConstant = -50;
                network.physics.options.forceAtlas2Based.centralGravity = 0.05;
                network.physics.options.forceAtlas2Based.springLength = 200;
                network.physics.options.forceAtlas2Based.springConstant = 0.05;
                
                // Update physics and stabilize
                network.setOptions({physics: network.physics.options});
                
                // Make sure physics is enabled
                if (!isPhysicsEnabled) {
                    network.setOptions({ physics: true });
                    isPhysicsEnabled = true;
                    updatePhysicsStatus(true);
                    releaseFixedPositions();
                }
                
                network.stabilize(100);
            } catch (e) {
                console.error("Error expanding nodes:", e);
            }
        }
        
        // Reset layout - 確保所有節點還原到初始位置但保持固定
        function resetLayout() {
            if (!network || !nodesDataset) {
                console.warn("Cannot reset layout: network or nodes not available");
                return;
            }
            
            try {
                // 使用初始保存的位置重置所有節點
                if (Object.keys(savedPositions).length > 0) {
                    Object.keys(savedPositions).forEach(function(nodeId) {
                        const pos = savedPositions[nodeId];
                        nodesDataset.update({
                            id: nodeId,
                            x: pos.x,
                            y: pos.y,
                            fixed: {
                                x: true,
                                y: true
                            }
                        });
                    });
                }
                
                // 重置縮放
                network.fit();
                
                // 確保物理引擎關閉
                network.setOptions({ physics: false });
                isPhysicsEnabled = false;
                updatePhysicsStatus(false);
                
                console.log("Reset complete. All nodes restored to original positions.");
            } catch (e) {
                console.error("Error resetting layout:", e);
            }
        }
        
        // Freeze node positions
        function freezePositions() {
            if (!network || !nodesDataset) {
                console.warn("Cannot freeze positions: network or nodes not available");
                return;
            }
            
            try {
                // 儲存所有節點當前位置
                saveCurrentNodePositions();
                
                // 將所有基因節點設為固定位置
                const ids = nodesDataset.getIds();
                ids.forEach(id => {
                    // 修改所有節點為固定位置
                    const pos = savedPositions[id];
                    if (pos && typeof pos.x === 'number' && typeof pos.y === 'number') {
                        nodesDataset.update({
                            id: id,
                            x: pos.x,
                            y: pos.y,
                            fixed: {
                                x: true,
                                y: true
                            }
                        });
                    }
                });
                
                // 物理引擎已關閉，更新狀態顯示
                updatePhysicsStatus(false);
                console.log("All nodes frozen in place");
            } catch (e) {
                console.error("Error freezing positions:", e);
            }
        }
        
        // Unfreeze node positions
        function unfreezePositions() {
            if (!network || !nodesDataset) {
                console.warn("Cannot unfreeze positions: network or nodes not available");
                return;
            }
            
            try {
                const ids = nodesDataset.getIds();
                
                ids.forEach(id => {
                    // 所有節點都可以自由移動
                    nodesDataset.update({
                        id: id,
                        fixed: {
                            x: false,
                            y: false
                        }
                    });
                });
                
                updatePhysicsStatus(true);
                console.log("All nodes can now be dragged freely");
            } catch (e) {
                console.error("Error unfreezing positions:", e);
            }
        }
        </script>
        '''
        
        # Add html2canvas library and download buttons with script
        html2canvas_lib = '''
        <script src="https://html2canvas.hertzen.com/dist/html2canvas.min.js"></script>
        '''
        
        # Create download buttons
        download_buttons = '''
        <div class="download-buttons">
            <button id="download-btn">Download PNG</button>
            <button id="download-hires-btn">Download Hi-Res PNG</button>
        </div>
        '''
        
        # Add JavaScript for downloading PNG with different resolutions
        download_script = '''
        <script>
        // Regular PNG download
        document.getElementById('download-btn').addEventListener('click', function() {
            downloadNetworkImage();
        });
        
        // High-resolution PNG download
        document.getElementById('download-hires-btn').addEventListener('click', function() {
            downloadNetworkImage(3.0); // 3x higher resolution
        });
        
        function downloadNetworkImage(scaleFactor = 1.0) {
            // Get the entire document for capturing, including legend and title
            const fullPage = document.querySelector('body');
            
            if (fullPage) {
                // Hide the download buttons during capture
                const downloadBtns = document.querySelector('.download-buttons');
                const controlBtns = document.querySelector('.control-buttons');
                const statusIndicator = document.querySelector('.status-indicator');
                
                if (downloadBtns) downloadBtns.style.display = 'none';
                if (controlBtns) controlBtns.style.display = 'none';
                if (statusIndicator) statusIndicator.style.display = 'none';
                
                // Use html2canvas to capture the visualization including legend and title
                html2canvas(fullPage, {
                    scale: scaleFactor,
                    backgroundColor: '#f9f9f9',
                    logging: false,
                    allowTaint: true,
                    useCORS: true
                }).then(canvas => {
                    // Create a temporary link element
                    const link = document.createElement('a');
                    link.href = canvas.toDataURL('image/png');
                    
                    // Set filename based on resolution
                    const filename = scaleFactor > 1.0 ? 'multi_tumor_network_hires.png' : 'multi_tumor_network.png';
                    link.download = filename;
                    
                    // Append, click, and remove the link
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    
                    // Restore visibility
                    if (downloadBtns) downloadBtns.style.display = 'flex';
                    if (controlBtns) controlBtns.style.display = 'flex';
                    if (statusIndicator) statusIndicator.style.display = 'block';
                });
            } else {
                alert('Could not find the visualization. Please wait until it is fully loaded.');
            }
        }
        </script>
        '''
        
        # Find the position of </head> to inject CSS
        head_end = html_content.find('</head>')
        if head_end != -1:
            html_content = html_content[:head_end] + css_styles + html2canvas_lib + html_content[head_end:]
        
        # Find the body tag to wrap network in our custom container
        body_start = html_content.find('<body')
        if body_start != -1:
            body_tag_end = html_content.find('>', body_start)
            if body_tag_end != -1:
                # Insert the start of container and title after opening body tag
                html_content = html_content[:body_tag_end + 1] + '\n<div class="network-container">\n' + title_html + html_content[body_tag_end + 1:]
        
        # Find where to insert the closing container tag
        body_end = html_content.find('</body>')
        if body_end != -1:
            # Add close of container div before body end
            html_content = html_content[:body_end] + '\n</div>\n' + legend_html + status_indicator + control_buttons + download_buttons + control_script + download_script + html_content[body_end:]
        
        # Write modified HTML back to file
        with open(html_path, 'w', encoding='utf-8', errors='ignore') as file:
            file.write(html_content)
        
        return html_path

    def run(self, data_dir: str = 'data') -> str:
        """Run the visualization process.
        
        Args:
            data_dir: Directory containing tumor data files
            
        Returns:
            Path to the generated HTML file
        """
        start_time = time.time()
        
        # Load data for all tumor types
        tumor_data = self.load_all_data(data_dir)
        
        # Create network
        net = self.create_network(tumor_data)
        
        # Create HTML file with title, legend, and enhanced download buttons
        html_path = os.path.join(self.output_dir, f'{self.output_filename}.html')
        self.create_html_with_title_and_legend(html_path, net)
        
        elapsed_time = time.time() - start_time
        logger.info(f"Visualization created in {elapsed_time:.2f} seconds")
        logger.info(f"Interactive visualization with enhanced features created at: {html_path}")
        
        print(f"\nInteractive visualization created at: {html_path}")
        print("Enhanced Features:")
        print("- Multi-tumor network visualization with GCH1 as central node")
        print("- Regular and high-resolution PNG download options")
        print("- Node clustering controls (bottom-left)")
        print("- Interactive exploration tools")
        print(f"- Visualization includes {sum(self.gene_counts.values())} genes across {len(self.tumor_types)} tumor types")
        
        return html_path


if __name__ == "__main__":
    visualizer = MultiTumorNetworkVisualizer()
    html_file = visualizer.run()
    
    # Open the file in the default browser
    import webbrowser
    try:
        webbrowser.open('file://' + os.path.abspath(html_file))
        print("\nOpening visualization in web browser...")
    except Exception as e:
        print(f"\nCould not open browser automatically: {e}")
        print(f"Please open the file manually: {html_file}") 