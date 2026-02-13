#!/usr/bin/env python
"""
Test script to verify mass conservation terminal output.
This runs a very short simulation to demonstrate the terminal printing functionality.
"""

import sys
import os

# Add src to path
sys.path.insert(0, 'src')

print("=" * 80)
print("TESTING MASS CONSERVATION TERMINAL OUTPUT")
print("=" * 80)
print()
print("This test demonstrates that mass conservation results are now printed")
print("to the terminal after simulation completion, in addition to being saved")
print("to a text file.")
print()
print("Running short simulation with mass conservation check enabled...")
print("=" * 80)
print()

# Import here to ensure sys.path is set first
from cs_model.main import run_model

# Run a short test simulation with mass conservation check enabled
results = run_model('test_short.yaml', check_mass_conservation=True)

print()
print("=" * 80)
print("TEST COMPLETED SUCCESSFULLY!")
print("=" * 80)
print()
print("Notice that:")
print("1. The detailed mass conservation report was printed above")
print("2. A detailed report file was saved to disk")
print("3. A final summary was printed at the end for quick reference")
print()
