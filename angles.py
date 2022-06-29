import numpy as np
import pandas as pd
import matplotlib.pyplot as pp
import analyse as ana
import argparse
from scipy.stats import gaussian_kde as kde 

resid_offset = 1 # output resid minus input resid

### Parameter sets

parsets = [
   'CTPOL', 'OPLS-AA', 'OPLS-AA', 'opt-CTPOL'
]

ref = 'trajectories/input.pdb'

### Inputs for angle calculating function

inputs = [
    dict(
        ang_label = "Angle of S-Zn-S",
        selections = [
            dict(name='SG', resid=4),
            dict(name='ZN'),
            dict(name='SG', resid=7)
        ]
    ),
    dict(
        ang_label = "Angle of N-Zn-N",
        selections = [
            dict(name='NE2', resid=20),
            dict(name='ZN'),
            dict(name='NE2', resid=24)
        ]
    )
]

### Get Universes


unis = dict(NMR=ana.Universe(ref, ref))
for par in parsets:
    dcd = f"trajectories/MD-1ZNF/{par}/centered_output.dcd"
    unis[par] = ana.Universe('output.pdb', dcd)

### Plot


for inp in inputs:
    fig, ax = pp.subplots()
    for par, uni in unis.items():
        sels = []
        for sel in inp['selections']:
            if 'resid' in sel:
                sel['resid'] -= (par == 'NMR')*resid_cutoff
            seltxt = ' and '.join(
                [' '.join([k, v] for k, v in sel.items()]
            )
            sels.append(uni.select_atoms(seltxt))
        
        angles = uni.get_angles(*sels)
        _, bins = np.histogram(angles)
        
        ang_range = np.linspace(bins[0], bins[-1], 200)

        ax.plot(ang_range, kde(angles)(x)) 
