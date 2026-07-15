#!/usr/bin/env python3
"""
Filling File Validator for Microsim_GP

Validates filling file syntax and geometry.
Run: python validate_filling.py <fill_file.in> <mesh_x> <mesh_y> <mesh_z>
"""

import sys
import re
import os
from typing import Dict, List, Tuple, Optional

# Valid fill keywords
VALID_FILL_TYPES = [
    'FILLSPHERE',
    'FILLSPHERERANDOM',
    'FILLCUBE',
    'FILLCYLINDER',
    'FILLCYLINDERRANDOM',
    'FILLVORONOI2D',
    'FILLVORONOI3D',
    'FILLCIRCLE',
    'FILLVELOCITYCUBE',
]

# Required parameters per fill type (position-based)
FILL_PARAMS = {
    'FILLSPHERE': 5,           # x, y, z, radius, phase
    'FILLSPHERERANDOM': 4,     # count, r_min, r_max, phase
    'FILLCUBE': 7,             # x_start, x_end, y_start, y_end, z_start, z_end, phase
    'FILLCYLINDER': 5,         # x, y, radius, height, phase
    'FILLCYLINDERRANDOM': 4,   # count, r_min, r_max, phase
    'FILLVORONOI2D': 1,        # num_grains
    'FILLVORONOI3D': 1,        # num_grains
    'FILLCIRCLE': 4,           # x, y, radius, phase
    'FILLVELOCITYCUBE': 7,     # x_start, x_end, y_start, y_end, z_start, z_end
}


class FillingValidator:
    def __init__(self, filepath: str, mesh_size: Tuple[int, int, int] = None):
        self.filepath = filepath
        self.mesh_size = mesh_size or (None, None, None)
        self.fills: List[Dict] = []
        self.errors: List[str] = []
        self.warnings: List[str] = []
        
    def parse(self) -> bool:
        """Parse filling file and extract fill operations."""
        if not os.path.exists(self.filepath):
            self.errors.append(f"File not found: {self.filepath}")
            return False
            
        with open(self.filepath, 'r') as f:
            lines = f.readlines()
            
        for line_num, line in enumerate(lines, 1):
            # Remove comments
            line = re.sub(r'#.*$', '', line).strip()
            if not line:
                continue
                
            # Parse fill command: FILLTYPE = {params};
            match = re.match(r'(\w+)\s*=\s*\{([^}]*)\};', line)
            if not match:
                if '=' in line:
                    self.warnings.append(f"Line {line_num}: Unrecognized format, ignoring: {line[:50]}")
                continue
                
            fill_type = match.group(1).strip()
            params_str = match.group(2).strip()
            
            if fill_type not in VALID_FILL_TYPES:
                self.errors.append(f"Line {line_num}: Unknown fill type '{fill_type}'")
                continue
                
            # Parse parameters
            if params_str:
                params = [p.strip() for p in params_str.split(',')]
            else:
                params = []
                
            self.fills.append({
                'type': fill_type,
                'params': params,
                'line': line_num
            })
            
        return True
        
    def validate(self) -> bool:
        """Run all validation checks."""
        self._check_fill_param_counts()
        self._check_phase_indices()
        self._check_geometry_bounds()
        self._check_overlaps()
        
        return len(self.errors) == 0
        
    def _check_fill_param_counts(self):
        """Check each fill has correct number of parameters."""
        for fill in self.fills:
            fill_type = fill['type']
            params = fill['params']
            expected = FILL_PARAMS.get(fill_type)
            
            if expected and len(params) != expected:
                self.errors.append(
                    f"Line {fill['line']}: {fill_type} requires {expected} parameters, got {len(params)}"
                )
                
    def _check_phase_indices(self):
        """Check phase indices are valid."""
        # Need NUMPHASES to validate - use placeholder if not provided
        max_phase = 10  # Assume reasonable max for warning
        
        for fill in self.fills:
            fill_type = fill['type']
            params = fill['params']
            
            # Find phase index position for each fill type
            phase_pos = {
                'FILLSPHERE': 4,
                'FILLCUBE': 6,
                'FILLCYLINDER': 4,
                'FILLCIRCLE': 3,
            }
            
            if fill_type in phase_pos:
                pos = phase_pos[fill_type]
                if pos < len(params):
                    try:
                        phase = int(params[pos])
                        if phase < 0:
                            self.errors.append(
                                f"Line {fill['line']}: Phase index must be >= 0, got {phase}"
                            )
                        # Note: Can't fully validate without NUMPHASES
                    except ValueError:
                        self.errors.append(
                            f"Line {fill['line']}: Phase must be integer, got '{params[pos]}'"
                        )
                        
    def _check_geometry_bounds(self):
        """Check geometry is within mesh bounds."""
        mesh_x, mesh_y, mesh_z = self.mesh_size
        
        for fill in self.fills:
            params = fill['params']
            fill_type = fill['type']
            
            try:
                if fill_type == 'FILLSPHERE':
                    x, y, z, r, phase = [float(p) for p in params]
                    if mesh_x and (x + r > mesh_x or x - r < 0):
                        self.warnings.append(
                            f"Line {fill['line']}: Sphere may extend beyond mesh X bounds"
                        )
                    if mesh_y and (y + r > mesh_y or y - r < 0):
                        self.warnings.append(
                            f"Line {fill['line']}: Sphere may extend beyond mesh Y bounds"
                        )
                        
                elif fill_type == 'FILLCUBE':
                    xs, xe, ys, ye, zs, ze, phase = [float(p) for p in params]
                    if mesh_x and (xe > mesh_x or xs < 0):
                        self.errors.append(
                            f"Line {fill['line']}: Cube X bounds [{xs}, {xe}] exceed mesh [0, {mesh_x}]"
                        )
                    if mesh_y and (ye > mesh_y or ys < 0):
                        self.errors.append(
                            f"Line {fill['line']}: Cube Y bounds [{ys}, {ye}] exceed mesh [0, {mesh_y}]"
                        )
                    if mesh_z and (ze > mesh_z or zs < 0):
                        self.errors.append(
                            f"Line {fill['line']}: Cube Z bounds [{zs}, {ze}] exceed mesh [0, {mesh_z}]"
                        )
                        
                elif fill_type == 'FILLSPHERERANDOM':
                    count, r_min, r_max = int(params[0]), float(params[1]), float(params[2])
                    if r_min >= r_max:
                        self.errors.append(
                            f"Line {fill['line']}: r_min ({r_min}) must be < r_max ({r_max})"
                        )
                        
            except (ValueError, IndexError) as e:
                self.warnings.append(f"Line {fill['line']}: Could not validate geometry: {e}")
                
    def _check_overlaps(self):
        """Check for obvious overlaps (basic check)."""
        # This is a simplified check - full overlap detection is complex
        spheres = []
        
        for fill in self.fills:
            if fill['type'] == 'FILLSPHERE':
                try:
                    x, y, z, r, phase = [float(p) for p in fill['params']]
                    spheres.append((x, y, z, r, phase, fill['line']))
                except:
                    pass
                    
        # Check sphere-sphere overlaps
        for i in range(len(spheres)):
            for j in range(i + 1, len(spheres)):
                x1, y1, z1, r1, p1, line1 = spheres[i]
                x2, y2, z2, r2, p2, line2 = spheres[j]
                
                # Distance between centers
                dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ** 0.5
                
                if dist < (r1 + r2) and p1 == p2:
                    self.warnings.append(
                        f"Lines {line1}, {line2}: Potential overlap between spheres "
                        f"(distance: {dist:.1f}, sum of radii: {r1+r2:.1f})"
                    )
                    
    def print_results(self):
        """Print validation results."""
        mesh_str = f"{self.mesh_size[0]}x{self.mesh_size[1]}x{self.mesh_size[2]}" if self.mesh_size[0] else "not specified"
        
        print("=" * 60)
        print(f"Filling File Validator: {self.filepath}")
        print(f"Mesh size: {mesh_str}")
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
            
        print(f"\nFound {len(self.fills)} fill operations:")
        for fill in self.fills:
            print(f"  Line {fill['line']}: {fill['type']} = {{{', '.join(fill['params'])}}}")
            
        print("=" * 60)
        return len(self.errors) == 0


def main():
    if len(sys.argv) < 2:
        print("Usage: python validate_filling.py <fill_file.in> [mesh_x] [mesh_y] [mesh_z]")
        print("Example: python validate_filling.py Fill.in 200 200 1")
        sys.exit(1)
        
    filepath = sys.argv[1]
    
    # Parse optional mesh size
    mesh_size = None
    if len(sys.argv) >= 5:
        try:
            mesh_size = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
        except ValueError:
            print("Error: Mesh dimensions must be integers")
            sys.exit(1)
            
    validator = FillingValidator(filepath, mesh_size)
    if validator.parse():
        valid = validator.validate()
        validator.print_results()
        sys.exit(0 if valid else 1)
    else:
        validator.print_results()
        sys.exit(1)


if __name__ == '__main__':
    main()
