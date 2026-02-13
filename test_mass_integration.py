#!/usr/bin/env python
"""
Test script for mass conservation integration.
This runs a short simulation and checks mass conservation.
"""

import sys
import os

# Add src to path
sys.path.insert(0, 'src')

from cs_model.main import run_model

print("=" * 80)
print("Testing Mass Conservation Integration")
print("=" * 80)
print()

# Run a short test simulation with mass conservation check enabled
print("Running short test simulation with mass conservation check...")
print()

results = run_model('test_short.yaml', check_mass_conservation=True)

print()
print("=" * 80)
print("Test completed successfully!")
print("=" * 80)
