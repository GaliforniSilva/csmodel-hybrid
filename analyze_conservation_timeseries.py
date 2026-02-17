#!/usr/bin/env python
"""
Analyze conservation time-series from simulation results.

This script demonstrates how to use the conservation time-series data
to track where and when conservation issues occur during a simulation.

Usage:
    python analyze_conservation_timeseries.py <results_file.pkl>
"""

import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt


def analyze_conservation_timeseries(results_file):
    """
    Analyze and visualize conservation time-series from simulation results.
    
    Parameters
    ----------
    results_file : str
        Path to the pickle file containing simulation results
    """
    
    # Load results
    print(f"Loading results from: {results_file}")
    with open(results_file, 'rb') as f:
        results = pickle.load(f)
    
    # Check if conservation data exists
    if 'conservation' not in results:
        print("ERROR: No conservation data found in results file.")
        print("Make sure the simulation was run with the updated code that tracks conservation.")
        return
    
    print("✓ Conservation data found!")
    print()
    
    # Extract time series data
    # Conservation data includes initial state (t=0) plus all timesteps
    # So we need to create a time array that includes t=0
    time_array = np.array(results['time'])
    dt = time_array[1] - time_array[0] if len(time_array) > 1 else 0
    # Prepend t=0 to time array to match conservation data
    time_with_initial = np.concatenate([[0], time_array])
    time_hours = time_with_initial / 3600  # Convert to hours
    conservation = results['conservation']
    
    # Print summary statistics
    print("=" * 80)
    print("CONSERVATION TIME-SERIES SUMMARY")
    print("=" * 80)
    print(f"Total simulation time: {time_hours[-1]:.2f} hours ({time_hours[-1]/24:.2f} days)")
    print(f"Number of timesteps: {len(time_hours)}")
    print()
    
    # Find maximum errors for each compartment
    max_dune_error = np.max(np.abs(conservation['dune_error']))
    max_beach_error = np.max(np.abs(conservation['beach_error']))
    max_bar_error = np.max(np.abs(conservation['bar_error']))
    max_total_error = np.max(np.abs(conservation['total_error']))
    
    # Find when maximum errors occurred
    idx_max_dune = np.argmax(np.abs(conservation['dune_error']))
    idx_max_beach = np.argmax(np.abs(conservation['beach_error']))
    idx_max_bar = np.argmax(np.abs(conservation['bar_error']))
    idx_max_total = np.argmax(np.abs(conservation['total_error']))
    
    print("Maximum conservation errors:")
    print(f"  Dune:  {max_dune_error:.6e} m³/m at t={time_hours[idx_max_dune]:.2f} hours (step {idx_max_dune})")
    print(f"  Beach: {max_beach_error:.6e} m³/m at t={time_hours[idx_max_beach]:.2f} hours (step {idx_max_beach})")
    print(f"  Bar:   {max_bar_error:.6e} m³/m at t={time_hours[idx_max_bar]:.2f} hours (step {idx_max_bar})")
    print(f"  Total: {max_total_error:.6e} m³/m at t={time_hours[idx_max_total]:.2f} hours (step {idx_max_total})")
    print()
    
    # Calculate cumulative errors
    cumulative_dune_error = np.cumsum(np.abs(conservation['dune_error']))
    cumulative_beach_error = np.cumsum(np.abs(conservation['beach_error']))
    cumulative_bar_error = np.cumsum(np.abs(conservation['bar_error']))
    cumulative_total_error = np.cumsum(np.abs(conservation['total_error']))
    
    print("Cumulative absolute errors (end of simulation):")
    print(f"  Dune:  {cumulative_dune_error[-1]:.6e} m³/m")
    print(f"  Beach: {cumulative_beach_error[-1]:.6e} m³/m")
    print(f"  Bar:   {cumulative_bar_error[-1]:.6e} m³/m")
    print(f"  Total: {cumulative_total_error[-1]:.6e} m³/m")
    print()
    
    # Identify problematic timesteps (error > 0.001% of initial volume)
    threshold_pct = 0.001
    problematic_steps = []
    for i in range(len(time_hours)):
        if (conservation['dune_error_pct'][i] > threshold_pct or 
            conservation['beach_error_pct'][i] > threshold_pct or 
            conservation['bar_error_pct'][i] > threshold_pct):
            problematic_steps.append(i)
    
    if problematic_steps:
        print(f"Timesteps with errors > {threshold_pct}% of initial volume:")
        for step in problematic_steps[:10]:  # Show first 10
            print(f"  Step {step:4d} (t={time_hours[step]:7.2f}h): "
                  f"Dune={conservation['dune_error_pct'][step]:.4f}%, "
                  f"Beach={conservation['beach_error_pct'][step]:.4f}%, "
                  f"Bar={conservation['bar_error_pct'][step]:.4f}%")
        if len(problematic_steps) > 10:
            print(f"  ... and {len(problematic_steps)-10} more")
    else:
        print(f"✓ No timesteps with errors > {threshold_pct}% of initial volume")
    
    print()
    print("=" * 80)
    
    # Create visualization
    print("\nGenerating conservation time-series plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Conservation Time-Series Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Volume time-series
    ax = axes[0, 0]
    ax.plot(time_hours, conservation['dune_volume'], label='Dune', linewidth=2)
    ax.plot(time_hours, conservation['beach_volume'], label='Beach', linewidth=2)
    ax.plot(time_hours, conservation['bar_volume'], label='Bar', linewidth=2)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Volume (m³/m)')
    ax.set_title('Compartment Volumes Over Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Absolute errors time-series
    ax = axes[0, 1]
    ax.plot(time_hours, np.abs(conservation['dune_error']), label='Dune', linewidth=2)
    ax.plot(time_hours, np.abs(conservation['beach_error']), label='Beach', linewidth=2)
    ax.plot(time_hours, np.abs(conservation['bar_error']), label='Bar', linewidth=2)
    ax.plot(time_hours, np.abs(conservation['total_error']), label='Total', 
            linewidth=2, linestyle='--', color='black')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Absolute Error (m³/m)')
    ax.set_title('Conservation Errors Over Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # Plot 3: Percentage errors time-series
    ax = axes[1, 0]
    ax.plot(time_hours, conservation['dune_error_pct'], label='Dune', linewidth=2)
    ax.plot(time_hours, conservation['beach_error_pct'], label='Beach', linewidth=2)
    ax.plot(time_hours, conservation['bar_error_pct'], label='Bar', linewidth=2)
    ax.plot(time_hours, conservation['total_error_pct'], label='Total', 
            linewidth=2, linestyle='--', color='black')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Error (% of initial volume)')
    ax.set_title('Conservation Errors (% of Initial Volume)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Cumulative absolute errors
    ax = axes[1, 1]
    ax.plot(time_hours, cumulative_dune_error, label='Dune', linewidth=2)
    ax.plot(time_hours, cumulative_beach_error, label='Beach', linewidth=2)
    ax.plot(time_hours, cumulative_bar_error, label='Bar', linewidth=2)
    ax.plot(time_hours, cumulative_total_error, label='Total', 
            linewidth=2, linestyle='--', color='black')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Cumulative Absolute Error (m³/m)')
    ax.set_title('Cumulative Conservation Errors')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = results_file.replace('.pkl', '_conservation_analysis.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Conservation analysis plot saved to: {output_file}")
    
    return results


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_conservation_timeseries.py <results_file.pkl>")
        print("\nExample:")
        print("  python analyze_conservation_timeseries.py test_short.pkl")
        sys.exit(1)
    
    results_file = sys.argv[1]
    analyze_conservation_timeseries(results_file)
