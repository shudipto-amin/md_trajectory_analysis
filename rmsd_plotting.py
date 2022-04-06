import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as pp

from matplotlib import gridspec
from MDAnalysis.analysis.rms import RMSD
from scipy import stats
from matplotlib.ticker import AutoLocator
psf = 'output.pdb'
ref = 'trajectories/input.pdb'

trajes = {
    'OPLS-AA' : "trajectories/MD-1ZNF/OPLS-AA/centered_output.dcd",
    'opt-OPLS-AA': "trajectories/MD-1ZNF/opt-OPLS-AA/centered_output.dcd",
    'CTPOL' : "trajectories/MD-1ZNF/CTPOL/centered_output.dcd",
    'opt-CTPOL' : "trajectories/MD-1ZNF/opt-CTPOL/centered_output.dcd",
}

TimeStep = 0.02

selections = {
    "backbone" : "backbone",
    r"Zn-Cys$_2$His$_2$" : "name ZN or resname CYS or resname HIS",
    r"Zn-S(Cys$^-$)/N(His)" : "name ZN or name SG or (name NE2 and resname HIS) "
}

def get_density(vals):
    kernel = stats.gaussian_kde(vals)

    x = np.linspace(min(vals), max(vals))
    y = kernel(x)

    return x, y

unis = {param : mda.Universe(psf, dcd) for param, dcd in trajes.items()}
refuni = mda.Universe(ref)


fig = pp.figure()
fig.set_size_inches(8,8)
gs = gridspec.GridSpec(3,3)
gs.update(hspace=0.01, wspace=0.01)
for n, sel in enumerate(selections.items()):
    print(sel)
    if n == 0:
        ax_left = fig.add_subplot(gs[n, 0:2])
        ticks = ax_left.get_xticks
    else:
        ax_left = fig.add_subplot(gs[n, 0:2], sharex=ax_left, sharey=ax_left)
        
    ax_right = fig.add_subplot(gs[n, 2], sharey=ax_left)

    if n < len(selections) - 1 :
        ax_left.xaxis.set_visible(False)
    else:
        ax_left.set_xlabel(r"Time (ns)", fontsize=16)

    if n == (len(selections)) // 2:
        ax_left.set_ylabel(r"RMSD ($\AA$)", fontsize=16)

    

    ax_right.yaxis.set_visible(False)
    ax_right.xaxis.set_visible(False)

    ax_left.text(
        0.5, 0.7, sel[0], 
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax_left.transAxes,
        fontsize=14
    )
    for param, uni in unis.items():
        rmsd = RMSD(uni, select=sel[1])
        rmsd.run()
        time = rmsd.results['rmsd'][:,0]*TimeStep
        ax_left.plot(
            time,
            rmsd.results['rmsd'][:,2],
            label = param
        )

        x, y = get_density(rmsd.results['rmsd'][:,2])
        ax_right.plot(y, x)
    
    ref_rmsd = RMSD(refuni, select=sel[1])
    ref_rmsd.run()
    x, y = get_density(ref_rmsd.results['rmsd'][:,2])
    ax_right.plot(y, x, linestyle='--', color='k')
   
    ax_right.set_xlim(left=0)
    ax_left.set_xlim(min(time), max(time)) 
    #ax_left.legend()

fig.tight_layout()
pp.savefig("rmsd_density.png")
pp.show()
 

        


