#!/usr/bin/env python3
"""
Combined Validator for Microsim_GP

Runs all validators together.
Run: python validate_all.py <input_file.in> [fill_file.in] [--tdb-dir <tdbs_directory>]
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from validate_input import InputValidator
from validate_filling import FillingValidator
from validate_thermo import ThermoValidator


def main():
    if len(sys.argv) < 2:
        print("Usage: python validate_all.py <input_file.in> [fill_file.in] [--tdb-dir <tdbs_directory>]")
        print("Example: python validate_all.py Input.in Fill.in --tdb-dir ../../tdbs")
        sys.exit(1)
        
    input_file = sys.argv[1]
    fill_file = None
    tdb_dir = None
    
    # Parse arguments
    args = sys.argv[2:]
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == '--tdb-dir' and i + 1 < len(args):
            tdb_dir = args[i + 1]
            i += 2
        elif not arg.startswith('--'):
            fill_file = arg
            i += 1
        else:
            i += 1
            
    print("\n" + "=" * 70)
    print("MICROSIM_GP COMBINED VALIDATOR")
    print("=" * 70)
    
    all_valid = True
    
    # Validate input file
    print("\n[1/3] Validating input file...")
    if os.path.exists(input_file):
        input_validator = InputValidator(input_file)
        if input_validator.parse():
            if input_validator.validate():
                print("  ✓ Input file valid")
            else:
                all_valid = False
                input_validator.print_results()
        else:
            all_valid = False
            input_validator.print_results()
            
        # Extract mesh size for filling validation
        mesh_size = (
            input_validator.params.get('MESH_X'),
            input_validator.params.get('MESH_Y'),
            input_validator.params.get('MESH_Z')
        )
    else:
        print(f"  ✗ Input file not found: {input_file}")
        all_valid = False
        mesh_size = (None, None, None)
        
    # Validate filling file
    if fill_file:
        print("\n[2/3] Validating filling file...")
        if os.path.exists(fill_file):
            fill_validator = FillingValidator(fill_file, mesh_size)
            if fill_validator.parse():
                if fill_validator.validate():
                    print("  ✓ Filling file valid")
                else:
                    all_valid = False
                    fill_validator.print_results()
            else:
                all_valid = False
                fill_validator.print_results()
        else:
            print(f"  ✗ Filling file not found: {fill_file}")
            all_valid = False
    else:
        print("\n[2/3] No filling file to validate")
        
    # Validate thermodynamic parameters
    print("\n[3/3] Validating thermodynamic parameters...")
    if os.path.exists(input_file):
        thermo_validator = ThermoValidator(input_file, tdb_dir)
        if thermo_validator.parse_input():
            if thermo_validator.validate():
                print("  ✓ Thermodynamic parameters valid")
            else:
                all_valid = False
                thermo_validator.print_results()
        else:
            all_valid = False
            thermo_validator.print_results()
    else:
        print("  ✗ Cannot validate thermodynamics without input file")
        all_valid = False
        
    # Summary
    print("\n" + "=" * 70)
    if all_valid:
        print("✓ ALL VALIDATIONS PASSED")
        print("\nThe input and filling files appear to be correctly formatted.")
        print("You can run the simulation with:")
        print(f"  mpirun -np <procs> ./microsim_gp {input_file} {fill_file or '<fill_file>'} <output_name>")
    else:
        print("✗ VALIDATION FAILED")
        print("\nPlease fix the errors listed above before running the simulation.")
    print("=" * 70)
    
    sys.exit(0 if all_valid else 1)


if __name__ == '__main__':
    main()
