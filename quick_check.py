"""
Quick test to check if ramp_crest height is reasonable after fix.
"""
import sys
import pickle
import numpy as np

sys.path.insert(0, '/home/runner/work/csmodel-hybrid/csmodel-hybrid/src')

# Load results from short test
with open('cs_model_test_short_ramp.pkl', 'rb') as f:
    results = pickle.load(f)

print("=" * 80)
print("RAMP_CREST HEIGHT CHECK")
print("=" * 80)

# Get initial and final ramp_crest heights
initial_height = results['dune_z'][0]['ramp_crest']
final_height = results['dune_z'][-1]['ramp_crest']
dfoot = results['substrate_z'][4]

print(f"\nBase elevation (dfoot): {dfoot:.2f} m")
print(f"Initial ramp_crest height: {initial_height:.2f} m")
print(f"Final ramp_crest height: {final_height:.2f} m")
print(f"Height change: {final_height - initial_height:.2f} m")

# Check if heights are reasonable (should be above dfoot)
if initial_height < dfoot:
    print(f"\n❌ ERROR: Initial ramp_crest ({initial_height:.2f} m) is below dfoot ({dfoot:.2f} m)!")
else:
    print(f"\n✓ Initial ramp_crest is above dfoot ({initial_height - dfoot:.2f} m above)")

if final_height < dfoot:
    print(f"❌ ERROR: Final ramp_crest ({final_height:.2f} m) is below dfoot ({dfoot:.2f} m)!")
else:
    print(f"✓ Final ramp_crest is above dfoot ({final_height - dfoot:.2f} m above)")

# Check for any NaN or invalid values
heights = [results['dune_z'][i]['ramp_crest'] for i in range(len(results['dune_z']))]
if any(np.isnan(heights)):
    print(f"\n❌ ERROR: Found NaN values in ramp_crest heights!")
elif any(h < dfoot - 0.1 for h in heights):
    problematic_steps = [i for i, h in enumerate(heights) if h < dfoot - 0.1]
    print(f"\n❌ ERROR: Found {len(problematic_steps)} time steps where ramp_crest is below dfoot!")
    print(f"   First problematic step: {problematic_steps[0]}, height: {heights[problematic_steps[0]]:.2f} m")
else:
    print(f"\n✓ All ramp_crest heights are valid throughout the simulation")

# Calculate total sediment input
total_aeolian = sum(results['qs_ramp']) + sum(results['qs_crest']) + sum(results['qs_lee'])
print(f"\nTotal aeolian sediment input: {total_aeolian:.2f} m³/m")

print("\n" + "=" * 80)
