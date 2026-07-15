#!/usr/bin/env python3
"""
Simple HDF5 visualization script for microsim_gp phase-field simulations.
Requires: h5py, numpy, matplotlib

CRITICAL: Run this script from /app/test/ where the DATA/ directory is located.
Example:
    cd /app/test
    python /app/microsim_gp/scripts/postprocessing/visualize.py DATA/output_1000.h5
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Verify running from /app/test
if not os.path.exists('/app/test'):
    print("ERROR: /app/test/ does not exist")
    sys.exit(1)

# Change to /app/test/ if not already there
if os.getcwd() != '/app/test':
    os.chdir('/app/test')
    print(f"Changed working directory to: {os.getcwd()}")

def print_structure(fname):
    """Print contents of HDF5 file."""
    print(f"\n=== Contents of {fname} ===")
    with h5py.File(fname, 'r') as f:
        def print_attrs(name, obj):
            print(f"  {name}: {type(obj).__name__}")
            if isinstance(obj, h5py.Dataset):
                print(f"    shape: {obj.shape}, dtype: {obj.dtype}")
        f.visititems(print_attrs)

def plot_phase_fields(fname, output_prefix=None):
    """Plot phase fields from HDF5 file."""
    with h5py.File(fname, 'r') as f:
        phase_names = [k for k in f.keys() if not k.startswith('Mu_') 
                       and not k.startswith('Composition') and k not in ['T', 'x', 'y', 'z']]
        
        n_phases = len(phase_names)
        fig, axes = plt.subplots(1, n_phases, figsize=(4*n_phases, 4))
        if n_phases == 1:
            axes = [axes]
        
        for i, phase_name in enumerate(phase_names):
            data = f[phase_name][:]
            im = axes[i].imshow(data.T, origin='lower', cmap='viridis')
            axes[i].set_title(phase_name)
            axes[i].set_xlabel('X')
            axes[i].set_ylabel('Y')
            plt.colorbar(im, ax=axes[i])
        
        plt.tight_layout()
        if output_prefix:
            plt.savefig(f"{output_prefix}_phases.png", dpi=150)
            print(f"Saved {output_prefix}_phases.png")
        else:
            plt.show()
        plt.close()

def plot_composition(fname, output_prefix=None):
    """Plot composition fields from HDF5 file."""
    with h5py.File(fname, 'r') as f:
        comp_names = [k for k in f.keys() if k.startswith('Composition_')]
        
        if not comp_names:
            print("No composition fields found")
            return
            
        n_comps = len(comp_names)
        fig, axes = plt.subplots(1, n_comps, figsize=(4*n_comps, 4))
        if n_comps == 1:
            axes = [axes]
        
        for i, comp_name in enumerate(comp_names):
            data = f[comp_name][:]
            im = axes[i].imshow(data.T, origin='lower', cmap='plasma')
            axes[i].set_title(comp_name)
            axes[i].set_xlabel('X')
            axes[i].set_ylabel('Y')
            plt.colorbar(im, ax=axes[i])
        
        plt.tight_layout()
        if output_prefix:
            plt.savefig(f"{output_prefix}_composition.png", dpi=150)
            print(f"Saved {output_prefix}_composition.png")
        else:
            plt.show()
        plt.close()

def plot_chemical_potential(fname, output_prefix=None):
    """Plot chemical potential fields from HDF5 file."""
    with h5py.File(fname, 'r') as f:
        mu_names = [k for k in f.keys() if k.startswith('Mu_')]
        
        if not mu_names:
            print("No chemical potential fields found")
            return
            
        n_mus = len(mu_names)
        fig, axes = plt.subplots(1, n_mus, figsize=(4*n_mus, 4))
        if n_mus == 1:
            axes = [axes]
        
        for i, mu_name in enumerate(mu_names):
            data = f[mu_name][:]
            im = axes[i].imshow(data.T, origin='lower', cmap='coolwarm')
            axes[i].set_title(mu_name)
            axes[i].set_xlabel('X')
            axes[i].set_ylabel('Y')
            plt.colorbar(im, ax=axes[i])
        
        plt.tight_layout()
        if output_prefix:
            plt.savefig(f"{output_prefix}_mu.png", dpi=150)
            print(f"Saved {output_prefix}_mu.png")
        else:
            plt.show()
        plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize.py <hdf5_file> [output_prefix]")
        print("       python visualize.py --list <hdf5_file>")
        sys.exit(1)
    
    fname = sys.argv[1]
    
    if not os.path.exists(fname):
        print(f"Error: File not found: {fname}")
        sys.exit(1)
    
    if '--list' in sys.argv:
        print_structure(fname)
        return
    
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else None
    
    plot_phase_fields(fname, output_prefix)
    plot_composition(fname, output_prefix)
    plot_chemical_potential(fname, output_prefix)
    
    print(f"\nVisualization complete!")

if __name__ == '__main__':
    main()
