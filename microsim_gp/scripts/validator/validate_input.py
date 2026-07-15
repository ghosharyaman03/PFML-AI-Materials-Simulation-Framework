#!/usr/bin/env python3
"""
Input File Validator for Microsim_GP

Validates input file syntax and checks for required parameters.
Run: python validate_input.py <input_file.in>
"""

import sys
import re
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Required parameters for all simulations
REQUIRED_SCALARS = {
    'DIMENSION': ['int', [2, 3]],
    'MESH_X': ['int', None],
    'MESH_Y': ['int', None],
    'MESH_Z': ['int', None],
    'DELTA_X': ['float', None],
    'DELTA_Y': ['float', None],
    'DELTA_Z': ['float', None],
    'DELTA_t': ['float', None],
    'NUMPHASES': ['int', None],
    'NUMCOMPONENTS': ['int', None],
    'NTIMESTEPS': ['int', None],
    'SAVET': ['int', None],
    'Function_F': ['int', [1, 2, 3, 4, 5, 6]],
    # Phase-field parameters (required for all simulations)
    'epsilon': ['float', None],
    'tau': ['float', None],
    'GAMMA': ['matrix', None],
    'Tau': ['matrix', None],
    # Thermodynamic/function parameters
    'R': ['float', None],
    'V': ['float', None],
    # Temperature control
    'ISOTHERMAL': ['int', [0, 1]],
}

# Parameters required for specific FUNCTION_F values
FUNCTION_F_PARAMS = {
    1: ['A', 'ceq', 'cfill', 'slopes'],
    2: ['num_thermo_phases', 'tdbfname', 'tdb_phases', 'phase_map', 'ceq', 'cfill', 'c_guess'],
    3: ['A', 'ceq', 'cfill', 'slopes'],
    4: ['num_thermo_phases', 'tdbfname', 'tdb_phases', 'phase_map', 'ceq', 'cfill', 'c_guess'],
    5: ['Latent_heat', 'Thermal_conductivity'],
    6: ['num_thermo_phases', 'tdbfname', 'tdb_phases', 'phase_map', 'ceq', 'cfill', 'c_guess'],
}

# Optional feature flags
FEATURE_FLAGS = {
    'Function_anisotropy': [0, 1],
    'Function_W': [1, 2],
    'ELASTICITY': [0, 1],
    'LBM': [0, 1],
    'ISOTHERMAL': [0, 1],
    'GRAIN_GROWTH': [0, 1],
    'CONSTRAINED': [0, 1],
    'Shift': [0, 1],
    'DILUTE': [0, 1],
}

# Required when corresponding feature is enabled
CONDITIONAL_REQUIRED = {
    'Function_anisotropy': ['dab'],
    'ELASTICITY': ['EIGEN_STRAIN', 'VOIGT_ISOTROPIC', 'rho'],
    'LBM': ['NU_LBM', 'rho_LBM'],
    'Shift': ['Shiftj'],
    'CONSTRAINED': ['lambda'],
}

# EIDT parameters - requires eidt_mode > 0
EIDT_PARAMS = {
    'comp_ff': 'Far-field compositions',
    'comp_ff_rate': 'Composition rate (mode 1)',
    'rad_coeff': 'Radius coefficient (mode 2)',
    'step_coeff': 'Step coefficient (mode 2)',
    'step_zero': 'Step zero (mode 2)',
}

# Special conditional: ISOTHERMAL=1 requires T, ISOTHERMAL=0 requires Tempgrady
ISOTHERMAL_CONDITIONAL = {
    1: ['T'],
    0: ['Tempgrady'],
}

# Parameters required for all (not conditional)
ALWAYS_REQUIRED_MATRIX = ['DIFFUSIVITY']


class InputValidator:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.params: Dict[str, any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []
        
    def parse(self) -> bool:
        """Parse input file and extract parameters."""
        if not os.path.exists(self.filepath):
            self.errors.append(f"File not found: {self.filepath}")
            return False
            
        with open(self.filepath, 'r') as f:
            content = f.read()
            
        # Remove comments
        content = re.sub(r'#.*$', '', content, flags=re.MULTILINE)
        
        # Parse scalar parameters: PARAM = value;
        scalar_pattern = r'(\w+)\s*=\s*([^;{}]+);'
        for match in re.finditer(scalar_pattern, content):
            name = match.group(1).strip()
            value = match.group(2).strip()
            self.params[name] = self._parse_value(value)
            
        # Parse array parameters: PARAM = {value1, value2, ...};
        array_pattern = r'(\w+)\s*=\s*\{([^}]*)\};'
        for match in re.finditer(array_pattern, content):
            name = match.group(1).strip()
            values_str = match.group(2).strip()
            if values_str:
                values = [v.strip() for v in values_str.split(',')]
                if name in self.params:
                    if isinstance(self.params[name], list) and len(self.params[name]) > 0 and isinstance(self.params[name][0], list):
                        self.params[name].append(values)
                    else:
                        self.params[name] = [self.params[name], values]
                else:
                    self.params[name] = [values]
            else:
                self.params[name] = []
                
        return True
        
    def _parse_value(self, value: str) -> any:
        """Parse a scalar value to int or float."""
        value = value.strip()
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value
            
    def validate(self) -> bool:
        """Run all validation checks."""
        self._check_required_scalars()
        self._check_dimension_consistency()
        self._check_function_f_params()
        self._check_conditional_params()
        self._check_matrix_params()
        self._check_boundary_conditions()
        self._check_tempgrady_format()
        self._check_deprecated_params()
        
        return len(self.errors) == 0
        
    def _check_required_scalars(self):
        """Check all required parameters are present."""
        for param, (expected_type, valid_values) in REQUIRED_SCALARS.items():
            if param not in self.params:
                self.errors.append(f"Missing required parameter: {param}")
                continue
                
            value = self.params[param]
            
            # Check type
            if expected_type == 'int' and not isinstance(value, int):
                self.errors.append(f"{param} must be an integer, got: {type(value).__name__}")
            elif expected_type == 'float' and not isinstance(value, (int, float)):
                self.errors.append(f"{param} must be a number, got: {type(value).__name__}")
                
            # Check valid values if specified
            if valid_values and value not in valid_values:
                self.errors.append(f"{param} must be one of {valid_values}, got: {value}")
                
    def _check_dimension_consistency(self):
        """Check dimension-related consistency."""
        dim = self.params.get('DIMENSION')
        mesh_z = self.params.get('MESH_Z', 1)
        
        if dim == 2 and mesh_z != 1:
            self.warnings.append(f"DIMENSION=2 but MESH_Z={mesh_z}. For 2D, MESH_Z should be 1.")
            
    def _check_function_f_params(self):
        """Check FUNCTION_F specific parameters."""
        func_f = self.params.get('Function_F')
        if func_f not in FUNCTION_F_PARAMS:
            return
            
        required = FUNCTION_F_PARAMS[func_f]
        for param in required:
            if param not in self.params:
                self.errors.append(f"Function_F={func_f} requires parameter: {param}")
                
    def _check_conditional_params(self):
        """Check conditional parameter requirements."""
        for feature, required_params in CONDITIONAL_REQUIRED.items():
            if self.params.get(feature):
                for param in required_params:
                    if param not in self.params:
                        self.errors.append(f"When {feature}=1, parameter required: {param}")
        
        isothermal = self.params.get('ISOTHERMAL')
        if isothermal is not None:
            required_params = ISOTHERMAL_CONDITIONAL.get(isothermal, [])
            for param in required_params:
                if param not in self.params:
                    if isothermal == 1:
                        self.errors.append(f"ISOTHERMAL=1 requires parameter: {param}")
                    else:
                        self.errors.append(f"ISOTHERMAL=0 (temperature gradient) requires parameter: {param}")
                        
        eidt_mode = self.params.get('eidt_mode')
        if eidt_mode is not None and eidt_mode > 0:
            for param in EIDT_PARAMS:
                if param not in self.params:
                    self.errors.append(f"eidt_mode={eidt_mode} requires parameter: {param}")
                        
    def _check_matrix_params(self):
        """Check matrix parameters that are always required."""
        for param in ALWAYS_REQUIRED_MATRIX:
            if param not in self.params:
                self.errors.append(f"Missing required parameter: {param}")
                
    def _check_boundary_conditions(self):
        """Validate BOUNDARY parameters."""
        boundary_list = self.params.get('BOUNDARY', [])
        
        if not boundary_list:
            self.warnings.append("No BOUNDARY specified - defaults to periodic (type 3)")
            return
            
        # Check format: {field, X+, X-, Y+, Y-, Z+, Z-}
        if boundary_list and isinstance(boundary_list[0], list):
            boundaries = boundary_list
        else:
            boundaries = [boundary_list] if boundary_list else []
            
        for i, bc in enumerate(boundaries):
            if len(bc) != 7:
                self.errors.append(f"BOUNDARY entry {i}: expected 7 elements, got {len(bc)}")
                
            # Check field type
            valid_fields = ['phi', 'mu', 'c', 'T', 'u', 'F']
            if bc[0] not in valid_fields:
                self.warnings.append(f"BOUNDARY field '{bc[0]}' may not be valid. Valid: {valid_fields}")
                
            # Check boundary types (0-3)
            for j in range(1, 7):
                try:
                    bc_type = int(bc[j])
                    if bc_type not in [0, 1, 2, 3]:
                        self.errors.append(f"BOUNDARY entry {i}, position {j}: type must be 0-3, got {bc_type}")
                except (ValueError, IndexError):
                    if len(bc) <= j:
                        self.errors.append(f"BOUNDARY entry {i}: missing value at position {j}")
                    else:
                        self.errors.append(f"BOUNDARY entry {i}, position {j}: must be integer")
                    
        # Check BOUNDARY_VALUE for Dirichlet (type 2)
        for bc in boundaries:
            if '2' in bc[1:]:  # Dirichlet on any face
                if 'BOUNDARY_VALUE' not in self.params:
                    self.errors.append("BOUNDARY has Dirichlet (2) but no BOUNDARY_VALUE specified")
                    break
                    
    def _check_tempgrady_format(self):
        """Validate Tempgrady parameter format if present."""
        tempgrady = self.params.get('Tempgrady')
        if tempgrady:
            if isinstance(tempgrady, list) and len(tempgrady) > 0:
                if isinstance(tempgrady[0], list):
                    if len(tempgrady[0]) != 5:
                        self.errors.append(f"Tempgrady: expected 5 elements, got {len(tempgrady[0])}")
                else:
                    if len(tempgrady) != 5:
                        self.errors.append(f"Tempgrady: expected 5 elements, got {len(tempgrady)}")
            else:
                self.errors.append("Tempgrady: must have 5 values {base_temp, DeltaT, Distance, offset, velocity}")
                
    def _check_deprecated_params(self):
        """Check for deprecated parameters."""
        deprecated = ['BINARY', 'TERNARY']
        for param in deprecated:
            if param in self.params:
                self.warnings.append(f"{param} is deprecated - auto-set from NUMCOMPONENTS")
                
    def print_results(self):
        """Print validation results."""
        print("=" * 60)
        print(f"Input File Validator: {self.filepath}")
        print("=" * 60)
        
        if self.errors:
            print(f"\nERRORS ({len(self.errors)}):")
            for err in self.errors:
                print(f"  ✗ {err}")
                
        if self.warnings:
            print(f"\nWARNINGS ({len(self.warnings)}):")
            for warn in self.warnings:
                print(f"  ⚠ {warn}")
                
        if not self.errors and not self.warnings:
            print("\n✓ All validations passed!")
            
        print("\nParsed Parameters:")
        for k, v in sorted(self.params.items()):
            if isinstance(v, list):
                print(f"  {k} = {{{', '.join(str(x) for x in v)}}}")
            else:
                print(f"  {k} = {v}")
                
        print("=" * 60)
        return len(self.errors) == 0


def main():
    if len(sys.argv) < 2:
        print("Usage: python validate_input.py <input_file.in>")
        sys.exit(1)
        
    filepath = sys.argv[1]
    
    validator = InputValidator(filepath)
    if validator.parse():
        valid = validator.validate()
        validator.print_results()
        sys.exit(0 if valid else 1)
    else:
        validator.print_results()
        sys.exit(1)


if __name__ == '__main__':
    main()
