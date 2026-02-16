#!/usr/bin/env python3
"""
Validation script for checking dune volume consistency in YAML configuration files.

This script verifies that the vramp value in dune_volumes matches the geometric
parameters in dune_parameters.

Usage:
    python validate_yaml_volumes.py <config_file.yaml>
"""

import sys
import yaml
import argparse


def calculate_vramp_from_geometry(ysprim, ys, betas):
    """
    Calculate ramp volume using the simplified triangular approximation.
    
    The ramp is approximated as a triangle with:
    - Base length: (ys - ysprim)
    - Height: betas × base (from slope)
    - Volume per unit width: 0.5 × base² × betas
    
    Parameters:
    -----------
    ysprim : float
        X-position of ramp crest (m)
    ys : float
        X-position of ramp toe (m)
    betas : float
        Ramp slope (dimensionless, dz/dx)
    
    Returns:
    --------
    float
        Ramp volume per unit width (m³/m)
    """
    base_length = ys - ysprim
    vramp = 0.5 * base_length**2 * betas
    return vramp


def validate_yaml_config(config_file):
    """
    Validate volume consistency in a YAML configuration file.
    
    Parameters:
    -----------
    config_file : str
        Path to YAML configuration file
    
    Returns:
    --------
    dict
        Dictionary with validation results
    """
    # Load YAML file
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract parameters
    dune_params = config.get('dune_parameters', {})
    dune_volumes = config.get('dune_volumes', {})
    
    ysprim = dune_params.get('ysprim')
    ys = dune_params.get('ys')
    betas = dune_params.get('betas')
    vramp_config = dune_volumes.get('vramp')
    
    # Check if all required parameters exist
    missing_params = []
    if ysprim is None:
        missing_params.append('dune_parameters.ysprim')
    if ys is None:
        missing_params.append('dune_parameters.ys')
    if betas is None:
        missing_params.append('dune_parameters.betas')
    if vramp_config is None:
        missing_params.append('dune_volumes.vramp')
    
    if missing_params:
        return {
            'valid': False,
            'error': f"Missing required parameters: {', '.join(missing_params)}"
        }
    
    # Calculate expected vramp
    vramp_calculated = calculate_vramp_from_geometry(ysprim, ys, betas)
    
    # Calculate error
    absolute_error = vramp_config - vramp_calculated
    relative_error = (absolute_error / vramp_calculated) * 100 if vramp_calculated != 0 else float('inf')
    
    # Determine if valid (within 1% tolerance)
    is_valid = abs(relative_error) < 1.0
    
    return {
        'valid': is_valid,
        'ysprim': ysprim,
        'ys': ys,
        'betas': betas,
        'base_length': ys - ysprim,
        'vramp_config': vramp_config,
        'vramp_calculated': vramp_calculated,
        'absolute_error': absolute_error,
        'relative_error': relative_error
    }


def print_validation_report(results):
    """Print a formatted validation report."""
    print("=" * 70)
    print("DUNE VOLUME VALIDATION REPORT")
    print("=" * 70)
    print()
    
    if 'error' in results:
        print(f"ERROR: {results['error']}")
        return
    
    print("Configuration Parameters:")
    print(f"  ysprim (ramp crest x-position): {results['ysprim']:.2f} m")
    print(f"  ys (ramp toe x-position):       {results['ys']:.2f} m")
    print(f"  betas (ramp slope):             {results['betas']:.4f}")
    print(f"  Base length (ys - ysprim):      {results['base_length']:.2f} m")
    print()
    
    print("Volume Calculation:")
    print(f"  Formula: vramp = 0.5 × (ys - ysprim)² × betas")
    print(f"         = 0.5 × {results['base_length']:.2f}² × {results['betas']:.4f}")
    print(f"         = 0.5 × {results['base_length']**2:.2f} × {results['betas']:.4f}")
    print(f"         = {results['vramp_calculated']:.2f} m³/m")
    print()
    
    print("Comparison:")
    print(f"  vramp in config file:     {results['vramp_config']:.2f} m³/m")
    print(f"  vramp from geometry:      {results['vramp_calculated']:.2f} m³/m")
    print(f"  Absolute error:           {results['absolute_error']:+.2f} m³/m")
    print(f"  Relative error:           {results['relative_error']:+.1f}%")
    print()
    
    if results['valid']:
        print("✓ VALIDATION PASSED: Volume is consistent with geometry")
    else:
        print("✗ VALIDATION FAILED: Volume is inconsistent with geometry")
        print()
        print("RECOMMENDATION:")
        print(f"  Update dune_volumes.vramp from {results['vramp_config']:.2f} to {results['vramp_calculated']:.2f}")
        print()
        print("  In your YAML file, change:")
        print(f"    dune_volumes:")
        print(f"      vramp: {results['vramp_calculated']:.2f}")
    
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='Validate dune volume consistency in YAML configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python validate_yaml_volumes.py config.yaml
    
This will check if the vramp value matches the geometric parameters.
        """
    )
    parser.add_argument('config_file', help='Path to YAML configuration file')
    
    args = parser.parse_args()
    
    try:
        results = validate_yaml_config(args.config_file)
        print_validation_report(results)
        
        # Exit with error code if validation failed
        sys.exit(0 if results.get('valid', False) else 1)
        
    except FileNotFoundError:
        print(f"ERROR: Configuration file not found: {args.config_file}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR: Invalid YAML file: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
