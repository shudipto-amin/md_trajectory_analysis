import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf as RDF
from MDAnalysis import transformations
import json

import matplotlib.pyplot as pp

def get_data(psf, traj, sysinfo=None):
    '''
    Constructs Universe object from files.

    Arguments:
        psf     : Protein Structure File for topological information
        traj    : Trajectory (pdb, xtc, dcd, etc)
        sysinfo : Json-format file with dimensions of system

    returns: MDAnalysis.Universe object
    '''

    uni = mda.Universe(psf, traj)
    
    if sysinfo is not None:
        with open(sysinfo, 'r') as inp:
            data = json.load(inp)
                
        trans = transformations.boxdimensions.set_dimensions(data['dimensions'])
        uni.trajectory.add_transformations(trans)

    return uni
       
def get_rdf(uni, sel1, sel2, **kwargs):
    '''
    Calculates rdf for atom selections in a Universe object

    returns: MDAnalysis.analysis.rdf object, and atom selections
    '''
    g1 = uni.select_atoms(sel1)
    g2 = uni.select_atoms(sel2)

    rdf = RDF.InterRDF(g1, g2, **kwargs)
    rdf.run()
    return rdf, g1, g2

def plot_rdf_cdf(r, gab, Nab):
    fig, ax = pp.subplots()
    ax2 = ax.twinx()

    ax.plot(r, gab, color='grey')
    ax2.plot(r, Nab, color='black')

    ax.set_xlabel("$r$", fontsize=16)

    ax.set_ylabel("$g_{ab}(r)$", color='grey', fontsize=16)
    ax2.set_ylabel("$N_{ab}(r)$", color='black', fontsize=16)

    return fig, ax, ax2

def get_rdf_peaks(rdf):
    pass

def get_contact(uni, seltxt):
    sel = uni.select_atoms(seltxt, updating=True)
    print(uni)
    print(uni.trajectory)
    counts = []
    for ts in uni.trajectory:
        counts.append(sel.n_atoms)
    
    bin_edges = np.arange(min(counts) - 0.5, max(counts) + 1.5)
    
    return counts, bin_edges
    
 
if __name__ == "__main__":
    psf = 'output.pdb'
    traj = 'output.dcd'
    sysinfo = None

    uni = get_data(psf, traj, sysinfo)

    # Atom selections
    names = {
        'sulfurs' : 'name S*',
        'nitrogens' : 'name N*',
        'oxygens' : 'name O*',
        'all' : 'name S* or name N* or name O*'
    }

    # RDFs of each selection around Zn:
    for sel, seltxt in names.items():
        # Get distances (r), rdf (gab)
        rdf, g1, g2 = get_rdf(uni, 'name ZN', f'{seltxt}', nbins=200, range=(1,5))
        r = rdf.results['bins']
        gab = rdf.results['rdf'] 

        # Get average no. of particles within r, as function of r (Nab)
        Nab = np.cumsum(rdf.results['count']) / (g1.n_atoms * uni.trajectory.n_frames)
        
        # Plot gab and cdf vs r
        fig, ax, ax2 = plot_rdf_cdf(r, gab, Nab)
        fig.suptitle(sel)

        #pp.show()
        fig.savefig(f"rdf_and_coordination_{sel}.png")
    

    # Distribution of number of atoms around Zn within cutoff
    cutoff = 2.1  # choose based on RDF
    atom = "ZN"
    seltxt = f"around {cutoff} name {atom}"
    print(seltxt)
    counts, bin_edges = get_contact(uni, seltxt)

    fig, ax = pp.subplots()
    ax.set_xlabel("Coordination number")
    ax.set_ylabel("Frequency")
    ax.hist(counts, bin_edges)
    fig.savefig(f"contact_around_{atom}_{cutoff}.png")

    fig,ax = pp.subplots()
    ax.set_xlabel("Frame no.")
    ax.set_ylabel("Coordination number")
    ax.plot(counts)
    fig.savefig(f"cn_vs_frame_{atom}_{cutoff}.png")
    pp.show()
