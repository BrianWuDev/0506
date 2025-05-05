#!/usr/bin/env python
"""
Simple script to open the bad survival tree visualization in a web browser.
Use this if the visualization did not open automatically when running
the main script.
"""

import os
import webbrowser
from pathlib import Path

def open_visualization(file_name="bad_survival_tree.html"):
    """Open the visualization in the default web browser.
    
    Args:
        file_name: Name of the HTML file to open
    """
    output_dir = Path("output")
    html_path = output_dir / file_name
    
    if html_path.exists():
        file_url = f"file://{os.path.abspath(html_path)}"
        print(f"Opening {file_url} in your default web browser...")
        webbrowser.open(file_url)
    else:
        print(f"Error: Could not find {html_path}")
        print("Please run bad_survival_tree.py first to generate the visualization.")

if __name__ == "__main__":
    open_visualization()
    print("\nIf the browser did not open automatically, please manually navigate to:")
    print(os.path.abspath(Path("output") / "bad_survival_tree.html")) 