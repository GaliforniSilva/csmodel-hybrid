"""
Trace exactly what happens when the dune gets stuck.
"""
import pickle
import numpy as np

with open('cs_model_test9b.pkl', 'rb') as f:
    results = pickle.load(f)

# Find the timestep where volume stops growing
volumes = []
for i in range(len(results['dune_x'])):
    import sys
    sys.path.insert(0, 'src')
    from cs_model.functions import total_volume
    
    dune_x = results['dune_x'][i]
    dune_z = results['dune_z'][i]
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    
    vol = total_volume(
        [dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                              'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                              'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        [dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                              'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                              'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        substrate_x, substrate_z
    )
    volumes.append(vol)
    
    if i > 100 and abs(vol - volumes[i-1]) < 0.0001:  # Volume stopped changing
        print(f"Volume stopped changing at step {i}, t={results['time'][i]/3600:.1f}h")
        print(f"  Volume: {vol:.4f} m³/m")
        print(f"  ramp_crest: z={dune_z['ramp_crest']:.6f} m, x={dune_x['ramp_crest']:.6f} m")
        print(f"  seaward_crest: z={dune_z['seaward_crest']:.6f} m")
        print(f"  qs_ramp: {results['qs_ramp'][i]:.8f} m³/m")
        print(f"  substrate_z[3] (dike top): {substrate_z[3]:.6f} m")
        print(f"  dfoot: {substrate_z[4]:.6f} m")
        
        # Check the condition on line 1378
        s_hs = substrate_z[3] - substrate_z[4]
        ysprim_hs = substrate_x[3]
        ys_hs = substrate_x[4]
        betas = 0.1
        
        vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5)
        
        # Estimate vrampnew (cumulative volume so far)
        vrampnew = vol  # Approximate
        vol_extra = vrampnew - vol_to_fill_face
        
        vol_to_fill_top = (ysprim_hs - (2*(dune_z['seaward_crest']-dune_z['ramp_crest'])/betas))*((dune_z['seaward_crest']-dune_z['ramp_crest'])/2)
        
        print(f"\n  vol_to_fill_face: {vol_to_fill_face:.4f} m³/m")
        print(f"  vol_extra (vrampnew - vol_to_fill_face): {vol_extra:.4f} m³/m")
        print(f"  vol_to_fill_top: {vol_to_fill_top:.4f} m³/m")
        print(f"\n  Condition: vol_extra < vol_to_fill_top and vol_extra > 0")
        print(f"             {vol_extra:.4f} < {vol_to_fill_top:.4f} and {vol_extra:.4f} > 0")
        print(f"             = {vol_extra < vol_to_fill_top and vol_extra > 0}")
        
        if vol_extra <= 0:
            print(f"\n  ❌ PROBLEM: vol_extra <= 0, so line 1378 condition FAILS!")
            print(f"     The code should enter f_ramp_dike_top, but my added check prevents it!")
        
        break

