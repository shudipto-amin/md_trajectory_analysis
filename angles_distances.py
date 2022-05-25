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


args = parser.parse_args()

def plot_all(unis):
    ang_label = "Angle of S-Zn-S"
    f_ang, a_ang = pp.subplots()
    f_dist, a_dist = pp.subplots(3, sharex=True)
    f_dist.subplots_adjust(hspace=0.0)
    a_dist[-1].set_xlabel('Distance $\AA$')
    a_dist[1].set_ylabel('Density $(\AA^{-1})$')
    for a in a_dist[:-1]:
        a.get_xaxis().set_visible(False)

    a_ang.set_xlabel('Angle $^\circ$')
    a_ang.set_ylabel('Relative frequency')
    for label, uni in unis.items():
        if label == 'NMR':
            s1_resid = 3
            s2_resid = 6
        else:
            s1_resid = 4
            s2_resid = 7
        ### SELECTIONS: each selection should be a single atom selection
        s1 = uni.select_atoms(f'name SG and resid {s1_resid}')
        s2 = uni.select_atoms('name ZN')
        s3 = uni.select_atoms(f'name SG and resid {s2_resid}')

        ### Get list of angles and list of distances
        angles = uni.get_angles(s1, s2, s3)  # angle of s1-s2-s3, where s2 is at the vertex of angle
        dists12 = uni.get_distances(s1, s2)   # order doesn't matter
        dists23 = uni.get_distances(s3, s2)   # order doesn't matter
        dists13 = uni.get_distances(s3, s1)   # order doesn't matter
        dist_labels_data = {
            "Zn-C4.Sulfur distance": dists12,
            "Zn-C7.Sulfur distance": dists23,
            "Sulfur-Sulfur distance": dists13
        }

        ### Plot

        bins = 20 # uni.trajectory.n_frames//50
        kwargs = dict(
            bins=bins, histtype='step', density=True,
            label=label, linewidth=2) # common style for all plots

        if label == 'NMR':
            kwargs['color'] = 'k'
            kwargs['linewidth'] = 1

        a_ang.hist(angles, **kwargs)

        for n, (lab, dist) in enumerate(dist_labels_data.items()):
            a_dist[n].hist(dist, **kwargs)
            a_dist[n].text(0.5, 0.9, lab, 
                horizontalalignment='center', transform=a_dist[n].transAxes)
    
        
        f_dist.suptitle(f'Distances')
        f_ang.suptitle(ang_label)

    a_ang.legend()
    a_dist[0].legend(loc=1)
    #for a in a_dist:
    #    a.legend()
    return f_ang, f_dist
### Parameter sets

parsets = [
    'CTPOL', 'OPLS-AA', 'OPLS-AA', 'opt-CTPOL'
]

ref = 'trajectories/input.pdb'


### Get Universes
unis = dict(NMR = ana.Universe(ref, ref))
for par in parsets:
    dcd = f"trajectories/MD-1ZNF/{par}/centered_output.dcd"
    unis[par] = ana.Universe('output.pdb', dcd)

#for uni in [NMR]:
f_ang, f = plot_all(unis)
f_ang.savefig("CTPOL_ang.png")
f.savefig("CTPOL_dist.png")
