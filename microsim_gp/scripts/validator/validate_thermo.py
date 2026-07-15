#!/usr/bin/env python3
"""
Thermodynamic Parameters Validator for Microsim_GP

Validates thermodynamic parameters consistency.
Run: python validate_thermo.py <input_file.in> [--tdb-dir <tdbs_directory>]
"""

import sys
import re
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set


class ThermoValidator:
    def __init__(self, input_filepath: str, tdb_dir: str = None):
        self.input_filepath = input_filepath
        self.tdb_dir = tdb_dir or 'tdbs'
        self.params: Dict[str, any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []
        
    def parse_input(self) -> bool:
        """Parse input file to get thermodynamic parameters."""
        if not os.path.exists(self.input_filepath):
            self.errors.append(f"Input file not found: {self.input_filepath}")
            return False
            
        with open(self.input_filepath, 'r') as f:
            content = f.read()
            
        # Remove comments
        content = re.sub(r'#.*$', '', content, flags=re.MULTILINE)
        
        # Parse parameters
        for match in re.finditer(r'(\w+)\s*=\s*([^;{}]+);', content):
            name = match.group(1).strip()
            value = match.group(2).strip()
            self.params[name] = self._parse_value(value)
            
        for match in re.finditer(r'(\w+)\s*=\s*\{([^}]*)\};', content):
            name = match.group(1).strip()
            values_str = match.group(2).strip()
            if values_str:
                values = [v.strip() for v in values_str.split(',')]
                self.params[name] = values
            else:
                self.params[name] = []
                
        return True
        
    def _parse_value(self, value: str) -> any:
        """Parse a scalar value."""
        value = value.strip()
        try:
            if '.' in value:
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value
            
    def validate(self) -> bool:
        """Run thermodynamic validation checks."""
        func_f = self.params.get('Function_F')
        
        if func_f == 1:
            self._validate_parabolic()
        elif func_f == 2:
            self._validate_tdb_based()
        elif func_f == 3:
            self._validate_linearized()
        elif func_f in [4, 6]:
            self._validate_tdb_based()
            self._validate_custom_a()
        elif func_f == 5:
            self._validate_thermal()
            
        self._check_matrix_consistency()
        
        return len(self.errors) == 0
        
    def _validate_parabolic(self):
        """Validate FUNCTION_F=1 (parabolic) parameters."""
        required = ['A', 'ceq', 'cfill', 'slopes']
        for param in required:
            if param not in self.params:
                self.errors.append(f"FUNCTION_F=1 requires parameter: {param}")
                
        # Check A matrix format
        if 'A' in self.params:
            a_values = self.params['A']
            n_phases = self.params.get('NUMPHASES', 2)
            n_comp = self.params.get('NUMCOMPONENTS', 2)
            expected_count = n_phases * (1 + (n_comp - 1)**2)
            
            if len(a_values) < expected_count:
                self.warnings.append(
                    f"A matrix has {len(a_values)} values, expected at least {expected_count}"
                )
                
        # Check composition matrices
        for param in ['ceq', 'cfill', 'slopes']:
            if param in self.params:
                self._validate_thermo_matrix(param)
                
    def _validate_tdb_based(self):
        """Validate FUNCTION_F=2 (TDB-based) parameters."""
        required = ['num_thermo_phases', 'tdbfname', 'tdb_phases', 'phase_map', 'ceq', 'cfill']
        for param in required:
            if param not in self.params:
                self.errors.append(f"FUNCTION_F=2 requires parameter: {param}")
                
        # Check TDB file exists
        if 'tdbfname' in self.params:
            tdb_file = self.params['tdbfname']
            tdb_path = os.path.join(self.tdb_dir, tdb_file)
            if not os.path.exists(tdb_path):
                self.errors.append(f"TDB file not found: {tdb_path}")
                
        # Check phase mapping
        if 'tdb_phases' in self.params and 'phase_map' in self.params:
            tdb_phases = self.params['tdb_phases']
            phase_map = self.params['phase_map']
            
            if len(tdb_phases) != len(phase_map):
                self.errors.append(
                    f"tdb_phases ({len(tdb_phases)}) must have same length as phase_map ({len(phase_map)})"
                )
                
    def _validate_linearized(self):
        """Validate FUNCTION_F=3 (linearized) parameters."""
        required = ['A', 'ceq', 'cfill', 'slopes']
        for param in required:
            if param not in self.params:
                self.errors.append(f"FUNCTION_F=3 requires parameter: {param}")
                
        self._validate_thermo_matrix('slopes')
        
    def _validate_custom_a(self):
        """Validate custom A matrix for FUNCTION_F=4,6."""
        if 'A' not in self.params:
            self.errors.append(f"FUNCTION_F={self.params.get('Function_F')} requires A matrix")
            
    def _validate_thermal(self):
        """Validate FUNCTION_F=5 (thermal) parameters."""
        required = ['Latent_heat', 'Thermal_conductivity']
        for param in required:
            if param not in self.params:
                self.errors.append(f"FUNCTION_F=5 requires parameter: {param}")
                
        if 'Latent_heat' in self.params:
            lf = self.params['Latent_heat']
            if lf <= 0:
                self.errors.append(f"Latent_heat must be positive, got: {lf}")
                
        if 'Thermal_conductivity' in self.params:
            k = self.params['Thermal_conductivity']
            if k <= 0:
                self.errors.append(f"Thermal_conductivity must be positive, got: {k}")
                
    def _validate_thermo_matrix(self, param_name: str):
        """Validate thermodynamic matrix format."""
        values = self.params.get(param_name, [])
        n_phases = self.params.get('NUMPHASES', 2)
        n_comp = self.params.get('NUMCOMPONENTS', 2)
        
        # Each phase needs 2 + (NUMCOMPONENTS-1) values
        expected_per_phase = 2 + (n_comp - 1)
        expected_total = n_phases * n_phases * expected_per_phase
        
        if len(values) < expected_total:
            self.warnings.append(
                f"{param_name} has {len(values)} values, expected at least {expected_total} "
                f"for {n_phases} phases × {n_phases} phases × {expected_per_phase} values"
            )
            
    def _check_matrix_consistency(self):
        """Check consistency between phases/components and matrices."""
        n_phases = self.params.get('NUMPHASES')
        n_comp = self.params.get('NUMCOMPONENTS')
        
        if not n_phases or not n_comp:
            self.warnings.append("Cannot validate matrix sizes without NUMPHASES and NUMCOMPONENTS")
            return
            
        # GAMMA matrix: n_phases * (n_phases - 1) / 2
        if 'GAMMA' in self.params:
            gamma_count = len(self.params['GAMMA'])
            expected = n_phases * (n_phases - 1) // 2
            if gamma_count != expected:
                self.errors.append(
                    f"GAMMA has {gamma_count} values, expected {expected} "
                    f"(NUMPHASES={n_phases} → N×(N-1)/2)"
                )
                
        # Tau matrix: same as GAMMA
        if 'Tau' in self.params:
            tau_count = len(self.params['Tau'])
            expected = n_phases * (n_phases - 1) // 2
            if tau_count != expected:
                self.errors.append(
                    f"Tau has {tau_count} values, expected {expected}"
                )
                
        # dab/fab matrices for anisotropy
        for mat_name in ['dab', 'fab']:
            if mat_name in self.params:
                mat_count = len(self.params[mat_name])
                expected = n_phases * (n_phases - 1) // 2
                if mat_count != expected:
                    self.warnings.append(
                        f"{mat_name} has {mat_count} values, expected {expected}"
                    )
                    
    def print_results(self):
        """Print validation results."""
        print("=" * 60)
        print(f"Thermodynamic Validator: {self.input_filepath}")
        print("=" * 60)
        
        func_f = self.params.get('Function_F', 'Not specified')
        print(f"\nThermodynamic model: FUNCTION_F = {func_f}")
        
        if self.errors:
            print(f"\nERRORS ({len(self.errors)}):")
            for err in self.errors:
                print(f"  ✗ {err}")
                
        if self.warnings:
            print(f"\nWARNINGS ({len(self.warnings)}):")
            for warn in self.warnings:
                print(f"  ⚠ {warn}")
                
        if not self.errors and not self.warnings:
            print("\n✓ All thermodynamic validations passed!")
            
        # Print key thermodynamic parameters
        print("\nThermodynamic Parameters:")
        for key in ['Function_F', 'num_thermo_phases', 'tdbfname', 'A', 'ceq', 'cfill', 
                    'slopes', 'Latent_heat', 'Thermal_conductivity']:
            if key in self.params:
                val = self.params[key]
                if isinstance(val, list):
                    print(f"  {key} = {{{len(val)} values}}")
                else:
                    print(f"  {key} = {val}")
                    
        print("=" * 60)
        return len(self.errors) == 0


def main():
    if len(sys.argv) < 2:
        print("Usage: python validate_thermo.py <input_file.in> [--tdb-dir <tdbs_directory>]")
        print("Example: python validate_thermo.py Input.in --tdb-dir ../../tdbs")
        sys.exit(1)
        
    input_file = sys.argv[1]
    
    # Parse optional tdb directory
    tdb_dir = None
    if '--tdb-dir' in sys.argv:
        idx = sys.argv.index('--tdb-dir')
        if idx + 1 < len(sys.argv):
            tdb_dir = sys.argv[idx + 1]
            
    validator = ThermoValidator(input_file, tdb_dir)
    if validator.parse_input():
        valid = validator.validate()
        validator.print_results()
        sys.exit(0 if valid else 1)
    else:
        validator.print_results()
        sys.exit(1)


if __name__ == '__main__':
    main()
