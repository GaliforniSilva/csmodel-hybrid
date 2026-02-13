"""
Comprehensive mass conservation report for CS Model.
This script runs a simulation and generates a detailed mass conservation report.
"""
import sys
import yaml

sys.path.insert(0, 'src')
from cs_model.main import run_model
from cs_model.mass_conservation import check_mass_conservation


def run_with_mass_conservation_report(config_path):
    """
    Run the model and generate a comprehensive mass conservation report.
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file
        
    Returns
    -------
    bool
        True if mass is conserved, False otherwise
    """
    print("=" * 80)
    print("CS MODEL MASS CONSERVATION REPORT")
    print("=" * 80)
    print()
    
    # Load config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    output_file = config.get('output', {}).get('filename', 'cs_model_results.pkl')
    
    print(f"Configuration file: {config_path}")
    print(f"Output file: {output_file}")
    print()
    
    # Run the model
    print("Running model...")
    results = run_model(config_path, check_mass_conservation=False)
    print("Model run completed.")
    print()
    
    # Perform mass conservation check
    mass_check_result = check_mass_conservation(results, config, verbose=True)
    
    return mass_check_result['conserved']


if __name__ == "__main__":
    if len(sys.argv) < 2:
        config_path = 'test_mass_conservation.yaml'
    else:
        config_path = sys.argv[1]
    
    run_with_mass_conservation_report(config_path)
