import numpy as np
import pandas as pd
import matplotlib.pyplot as pp
import analyse as ana
import argparse


### Argument handler
parser = argparse.ArgumentParser(
    description = '''
Script to find and plot distribution of angles and distances of:
    Zn - S (distance)
    Zn - N (distance)
    S - Zn - S (angle)
    N - Zn - N (angle)
Modify the SELECTIONS in script if you want to change the selections above.'''
)

parser.add_argument(
    'top',
    help="Topology file (*.psf, *.pdb, etc) with atom names, connectivity, etc."
)

parser.add_argument(
    'coor',
    help="Coordinate files (*.dcd, *.pdb, etc)."
)

args = parser.parse_args()

### Get Universe
uni = ana.Universe(args.top, args.coor)

### SELECTIONS: each selection should be a single atom selection
s1 = uni.select_atoms('name SG and resid 4')
s2 = uni.select_atoms('name ZN')
s3 = uni.select_atoms('name SG and resid 7')

### Get list of angles and list of distances
angles = uni.get_angles(s1, s2, s3)  # angle of s1-s2-s3, where s2 is at the vertex of angle
dists12 = uni.get_distances(s1, s2)   # order doesn't matter
dists23 = uni.get_distances(s3, s2)   # order doesn't matter
dists13 = uni.get_distances(s3, s1)   # order doesn't matter

### Plot
ang_label = "Angle of S-Zn-S"
dist_labels_data = {
    "Zn-C4.Sulfur distance": dists12,
    "Zn-C7.Sulfur distance": dists23,
    "Sulfur-Sulfur distance": dists13
}

bins = uni.trajectory.n_frames//50
kwargs = dict(bins=bins, histtype='step', linewidth=2) # common style for all plots

f, a = pp.subplots()
a.hist(angles, **kwargs)
a.set_xlabel('Angle $^\circ$')
a.set_ylabel('Relative frequency')
f.suptitle(ang_label)

f, a = pp.subplots()
for label, dist in dist_labels_data.items():
    a.hist(dist, label=label, **kwargs)
a.set_xlabel('Distance $\AA$')
a.set_ylabel('Relative frequency')
f.suptitle('Distances')

a.legend()
pp.show()
