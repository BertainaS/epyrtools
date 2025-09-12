#!/usr/bin/env python3
"""
Script to fix voigtian function calls in the notebook to use correct parameter format.
The EPyR voigtian function uses widths=(gaussian_fwhm, lorentzian_fwhm) instead of sigma=X, gamma=Y.
"""

import json
import re

def fix_voigt_calls_in_cell(cell_source):
    """Fix voigtian calls in a single cell's source code."""
    
    # Pattern to match voigtian calls with sigma and gamma parameters
    # This handles various formats: sigma=X, gamma=Y or sigma=var1, gamma=var2
    pattern = r'voigtian\(([^)]*?),\s*sigma=([^,)]+),\s*gamma=([^,)]+)([^)]*)\)'
    
    def replacement(match):
        prefix = match.group(1)  # Everything before sigma
        sigma_val = match.group(2).strip()
        gamma_val = match.group(3).strip()
        suffix = match.group(4)  # Everything after gamma (like derivative=X)
        
        # Build the replacement with widths tuple
        if suffix:
            return f'voigtian({prefix}, widths=({sigma_val}, {gamma_val}){suffix})'
        else:
            return f'voigtian({prefix}, widths=({sigma_val}, {gamma_val}))'
    
    # Apply the replacement
    fixed_source = re.sub(pattern, replacement, cell_source)
    
    return fixed_source

def fix_notebook():
    """Fix the notebook file."""
    
    notebook_path = '/Users/sylvainbertaina/Documents/Cloud_CNRS/GitHub/epyrtools/examples/notebooks/07_Lineshapes_Functions_Complete.ipynb'
    
    # Read the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Track changes
    changes_made = 0
    
    # Fix each cell
    for cell in notebook['cells']:
        if cell['cell_type'] == 'code' and 'source' in cell:
            # Get the source as a string
            if isinstance(cell['source'], list):
                original_source = ''.join(cell['source'])
            else:
                original_source = cell['source']
            
            # Apply fixes
            fixed_source = fix_voigt_calls_in_cell(original_source)
            
            if fixed_source != original_source:
                # Convert back to list format (preserving line structure)
                cell['source'] = fixed_source.splitlines(keepends=True)
                changes_made += 1
                
                print(f"Fixed cell with changes:")
                print(f"  Original: {repr(original_source[:100])}...")
                print(f"  Fixed:    {repr(fixed_source[:100])}...")
                print()
    
    # Write the fixed notebook back
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1, ensure_ascii=False)
    
    print(f"✅ Fixed {changes_made} cells with voigtian parameter issues")
    return changes_made

if __name__ == "__main__":
    fix_notebook()