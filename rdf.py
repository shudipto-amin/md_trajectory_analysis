import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp
import analyse as ana


# Color-blind friendly pallette
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle)


def plot_all(trajes, save=False, show=True):
    psf = 'output.pdb'
    sysinfo = None
    outfile = 'data.csv'

    # Atom selections
    names = {
        #'sulfurs' : 'name S*',
        #'nitrogens' : 'name N*',
        'oxygens' : 'name O*',
        'water_O' : 'resname HOH and name O*',
        'backbone_O' : 'backbone and name O*',
        'acid_O' : '(not backbone) and (resname ASP or resname GLU) and name O*' 
    }
    unis = {name:ana.Universe(psf, traj, outfile=outfile, sysinfo=None) for name, traj in trajes.items()}

    for n, (sel, seltxt) in enumerate(names.items()):
        fig, ax = pp.subplots()
        fig.suptitle(sel)
        ax2 = ax.twinx()
        ax.set_xlabel("$r$", fontsize=16)

        ax.set_ylabel("$g_{ab}(r)$", fontsize=16)
        ax2.set_ylabel("$N_{ab}(r)$", fontsize=16)
        # RDFs of each selection around Zn:
        for name, traj in trajes.items():
            uni = unis[name] 
            # Get distances (r), rdf (gab)
            rdf, g1, g2 = uni.get_rdf('name ZN', f'{seltxt}', nbins=300, range=(1,7))
            r = rdf.results['bins']
            gab = rdf.results['rdf'] 

            # Get average no. of particles within r, as function of r (Nab)
            Nab = np.cumsum(rdf.results['count']) / (g1.n_atoms * uni.trajectory.n_frames)
           
            # Convolve rdf to get smoother curve
            width = 5 # no. of points to use in convolution
            conv = np.ones(width)/width

            smooth_gab = np.convolve(conv, gab, mode='valid')
            smooth_r = np.convolve(conv, r, mode='valid')
 
            # Plot gab and cdf vs r

            #ax.plot(r, gab)
            ax.plot(smooth_r, smooth_gab, linewidth=2, label=name)
            ax2.plot(r, Nab, alpha=0.6, linestyle='--', label=name)

        ax.legend(loc='upper left')
        ax2.legend(loc='lower right')
            #pp.show()
        if save:
            fig.savefig(f"rdf_and_coordination_{sel}.png")
    
    if show:
        pp.show()
    

parsets = [
    'CTPOL', 'opt-CTPOL'
]
trajes = {
    par:f"trajectories/MD-1ZNF/{par}/centered_output.dcd" for par in parsets
}

plot_all(trajes, save=True, show=False)
