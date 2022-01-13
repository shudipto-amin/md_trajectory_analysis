import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf as RDF
from MDAnalysis import transformations
import json

import matplotlib.pyplot as pp

def get_data(psf, traj, sysinfo):
    '''
    Constructs Universe object from files.

    Arguments:
        psf     : Protein Structure File for topological information
        traj    : Trajectory (pdb, xtc, dcd, etc)
        sysinfo : Json-format file with dimensions of system

    returns: MDAnalysis.Universe object
    '''

    uni = mda.Universe(psf, traj)
    
    with open(sysinfo, 'r') as inp:
        data = json.load(inp)
            
    trans = transformations.boxdimensions.set_dimensions(data['dimensions'])
    uni.trajectory.add_transformations(trans)

    return uni
       
def get_rdf(uni, atom1, atom2, **kwargs):
    '''
    Calculates rdf for atom selections in a Universe object

    returns: MDAnalysis.analysis.rdf object, and atom selections
    '''
    g1 = uni.select_atoms(f'name {atom1}')
    g2 = uni.select_atoms(f'name {atom2}')

    rdf = RDF.InterRDF(g1, g2, **kwargs)
    rdf.run()
    return rdf, g1, g2

 
if __name__ == "__main__":
    psf = 'example.psf'
    traj = 'example.pdb'
    sysinfo = 'sysinfo.dat'

    uni = get_data(psf, traj, sysinfo)

    # Get distances (r), rdf (gab)
    rdf, g1, g2 = get_rdf(uni, 'OH2', 'OH2', nbins=200, range=(1,5))
    r = rdf.results['bins']
    gab = rdf.results['rdf'] 

    # Get average no. of particles within r, as function of r (Nab)
    Nab = np.cumsum(rdf.results['count']) / g1.n_atoms
    
    # Plot gab and cdf vs r
    fig, ax = pp.subplots()
    ax2 = ax.twinx()

    ax.plot(r, gab, color='grey')
    ax2.plot(r, Nab, color='black')

    ax.set_xlabel("$r$", fontsize=16)

    ax.set_ylabel("$g_{ab}(r)$", color='grey', fontsize=16)
    ax2.set_ylabel("$N_{ab}(r)$", color='black', fontsize=16)

    fig.savefig("rdf_and_coordination.png")
    

        
