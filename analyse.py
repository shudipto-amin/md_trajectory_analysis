import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf as RDF, rms as RMS
from MDAnalysis import transformations

import json
import pandas as pd
import matplotlib.pyplot as pp

class InterRDF(RDF.InterRDF):
    '''
    Modifies RDF.InterRDF() to accept option to not normalize by volume of box.
    This is becase the original version fails if this volume is not provided.
    This new method returns rdf which is essentially No. of particles per volume, per frame.
    '''
    def __init__(self, *args, normalize_volume=False, **kwargs):
        self.normalize_volume = normalize_volume
        RDF.InterRDF.__init__(self, *args, **kwargs)
        
    def _conclude(self):
        # Number of each selection
        nA = len(self.g1)
        nB = len(self.g2)
        N = nA * nB

        # If we had exclusions, take these into account
        if self._exclusion_block:
            xA, xB = self._exclusion_block
            nblocks = nA / xA
            N -= xA * xB * nblocks

        # Volume in each radial shell
        vols = np.power(self.results.edges, 3)
        vol = 4/3 * np.pi * np.diff(vols)

        # Average number density
        if self.normalize_volume:
            box_vol = self.volume / self.n_frames
            density = N / box_vol
        else:
            density = 1

        rdf = self.results.count / (density * vol * self.n_frames)
        self.results.rdf = rdf

        
class InterRDF_s(RDF.InterRDF_s):
        
    def _conclude(self):
        # Volume in each radial shell
        vols = np.power(self.results.edges, 3)
        vol = 4/3 * np.pi * np.diff(vols)

        # Empty lists to restore indices, RDF
        indices = []
        rdf = []

        for i, (ag1, ag2) in enumerate(self.ags):
            # Number of each selection
            indices.append([ag1.indices, ag2.indices])

            # Average number density
            

            if self._density:
                rdf.append(self.results.count[i] / (vol * self.n_frames))
            else:
                box_vol = self.volume / self.n_frames
                density = 1 / box_vol
                rdf.append(
                    self.results.count[i] / (density * vol * self.n_frames))

        self.results.rdf = rdf
        self.results.indices = indices
        
        
class Universe(mda.Universe):
    '''
    Class based on mda.Universe with additional methods for analysis
    '''
    def __init__(self, psf, traj, *args, outfile=None,sysinfo=None, **kwargs):
        '''
        Constructs Universe object from files.

        Arguments:
            psf     : <str> Protein Structure File for topological information
            traj    : <str> Trajectory (pdb, xtc, dcd, etc)
            outfile : <str> OpenMM output file as formatted by Xiaojuan
            sysinfo : <str> Json-format file with dimensions of system

        '''
        mda.Universe.__init__(self, psf, traj, *args, **kwargs)
        self.psffile = psf
        self.trajfile = traj

        if sysinfo is not None:
            with open(sysinfo, 'r') as inp:
                data = json.load(inp)
                    
            trans = transformations.boxdimensions.set_dimensions(data['dimensions'])
            self.trajectory.add_transformations(trans)

        if outfile is not None:
            self.outfile = outfile
            self.get_outputs()
         
    def get_contact(self, seltxt):
        '''
        Get no. of atoms which satisfy selection criteria, as a function of frame number.

        Arguments:
            seltxt  : <str> Selection criteria

        Returns:
            counts  : <list> No. of atoms per frame which match seltxt criterion
            bin_edges : <np.1darray> bin edges for plotting a histogram of counts
        '''

        sel = self.select_atoms(seltxt, updating=True)
        counts = []
        for ts in self.trajectory:
            counts.append(sel.n_atoms)

        bin_edges = np.arange(min(counts) - 0.5, max(counts) + 1.5)

        return counts, bin_edges       

    def get_macro(self):
        '''Calcualte/get macroscopic quantities (energy, pressure, volume) vs time.'''
        ### Must get energy and pressure from OpenMM output.
        self.volume = []
        for ts in self.trajectory:
            pos = self.atoms.positions
            sides = np.max(pos, axis=0) - np.min(pos, axis=0)
            V = np.product(sides)
            self.volume.append(V)

        self.rmsd = self._get_rmsd('protein')

    def get_outputs(self):
        '''Get outputs from output file as pd.DataFrame'''
        df = pd.read_csv(self.outfile)
        self.outputs = df
        
    def get_rdf(self, sel1, sel2, **kwargs):
        '''
        Calculates rdf for atom selections in a Universe object

        Arguments:
            sel1    : <str> Selection string
            sel2    : <str> Selection string
            kwargs  : keyword arguments to pass to InterRDF()

        returns: MDAnalysis.analysis.rdf object, and atom selection objects
        '''

        g1 = self.select_atoms(sel1)
        g2 = self.select_atoms(sel2)

        rdf = InterRDF(g1, g2, **kwargs)
        rdf.run()
        return rdf, g1, g2

    def write_traj(self, seltxt, fname, start=0, stop=-1, step=1):
        sel = self.select_atoms(seltxt)
        with mda.Writer(fname, sel.n_atoms) as out:
            for ts in self.trajectory[start : stop : step]:
                out.write(sel)

    def center_protein_write(self, fname, **kwargs):
        '''
        Adds transformation to center protein and saves to file.
        This slows down iteration, so use only for writing new trajectory to file, 
        then load that as new universe.
        '''
        protein = self.select_atoms('protein')
        not_protein = self.select_atoms('not protein')

        trans = [
            transformations.unwrap(protein),
            transformations.center_in_box(protein),
            transformations.wrap(not_protein)
        ]
        self.trajectory.add_transformations(*trans)
        
        self.write_traj('all', fname, **kwargs)

    def _get_rmsd(self, seltxt, ref_frame=0):
        '''Calculates rmsd vs frame no.'''
        sel = self.select_atoms(seltxt)
        rmsd = RMS.RMSD(sel)
        rmsd.run()
        return rmsd


    def get_angles(self, sel1, sel2, sel3):
        '''
        Each sel* must be a selection for ONE atom only.
        sel2 should be the selection for the atom in the middle.
        The angle will be calculated for sel1-sel2-sel3.
        '''

        norm = np.linalg.norm
        arccos = np.arccos
        angles = []
        for n, ts in enumerate(self.trajectory):
            p1 = sel1.positions[0]
            p2 = sel2.positions[0]
            p3 = sel3.positions[0]

            A = p1 - p2
            B = p3 - p2
            a = norm(A)
            b = norm(B)

            angle = arccos(np.dot(A, B) / (a*b)) * 180/np.pi
            angles.append(angle)

        return angles
        
    def get_distances(self, sel1, sel2):
        '''
        Get distances between two selections.
        If there is more than one atom in the selection, 
        '''

        dists = []
        
        for n, ts in enumerate(self.trajectory):
            dists.append([])
            for p1 in sel1.positions:
                for p2 in sel2.positions:
                    d = np.linalg.norm(p1 - p2)
                    dists[n].append(d)

        return np.array(dists)
    
    @staticmethod
    def get_normal(sel):
        '''Get the normal direction for a sel containing 3 atoms'''
        if sel.n_atoms != 3:
            raise ValueError(f"{sel} does not contain 3 atoms")

        a = sel[1].position - sel[0].position
        b = sel[2].position - sel[0].position

        return np.cross(a, b)

    def get_dihedral(self, sel1 ,sel2):
        '''Get dihedral between two planes defined by two atomgroups (3 atoms each)'''
        
        dihedrals = []
        for n, ts, in enumerate(self.trajectory):
            normal1 = self.get_normal(sel1)
            normal2 = self.get_normal(sel2)
            l1 = np.linalg.norm(normal1)
            l2 = np.linalg.norm(normal2)
            dihed = np.arccos(np.dot(normal1, normal2) / (l1*l2)) * 180/np.pi        

            dihedrals.append(dihed)
        
        return np.array(dihedrals)
    
    @staticmethod
    def get_bisector(sel):
        '''Get the bisector of angle between 3 atoms'''
        if sel.n_atoms != 3:
            raise ValueError(f"{sel} does not contain 3 atoms")
        
        a = sel[0].position - sel[1].position
        b = sel[2].position - sel[1].position
        
        a = a / np.linalg.norm(a)
        b = b / np.linalg.norm(b)
        
        return (a + b) / 2
    
    def get_bisector_angles(self, sel1, sel2):
        bi_angles = []
        for n, ts, in enumerate(self.trajectory):
            b1 = self.get_bisector(sel1)
            b2 = self.get_bisector(sel2)
            l1 = np.linalg.norm(b1)
            l2 = np.linalg.norm(b2)
            bi_ang = np.arccos(np.dot(b1, b2) / (l1*l2)) * 180/np.pi        

            bi_angles.append(bi_ang)
        
        return np.array(bi_angles)        
        
## Plotting functions
        
def plot_rdf_cdf(r, gab, Nab):
    '''
    Plots rdf (gab) and coordination (Nab) vs radial distance (r)
    '''
    fig, ax = pp.subplots()
    ax2 = ax.twinx()

    ax.plot(r, gab, color='grey')
    ax2.plot(r, Nab, color='black')

    ax.set_xlabel("$r$", fontsize=16)

    ax.set_ylabel("$g_{ab}(r)$", color='grey', fontsize=16)
    ax2.set_ylabel("$N_{ab}(r)$", color='black', fontsize=16)

    return fig, ax, ax2

def plot_macros(uni):
    fig, axes = pp.subplots(2,1, sharex=True)
    fig.set_size_inches(10,10)
    fig.subplots_adjust(hspace=0)
    # Volume and RMSD plot (top)
    ax = axes[0]
    ax2 = ax.twinx()

    ax.set_ylabel('Volume $\AA^3$')
    ax2.set_ylabel('RMSD $\AA$')

    ax.plot(uni.volume, color='k', linestyle='--',label="Volume")
    ax2.plot(uni.rmsd.results.rmsd[:,-1], color='k', alpha=0.6, linewidth=1.5,label="RMSD")

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()

    ax2.legend(lines + lines2, labels + labels2, loc=0)

    pot_ener =  uni.outputs['Potential Energy (kJ/mole)']
    kin_ener =  uni.outputs['Kinetic Energy (kJ/mole)']
    tot_ener = pot_ener + kin_ener

    # Energy plot (bottom)
    ax1 = axes[1]
    ax12 = ax1.twinx()

    ax1.plot(tot_ener, label="Total Energy", color='k', linestyle='--')
    ax12.plot(pot_ener, label="Potential Energy", color='k', alpha=0.6, linewidth=2 )
    ax1.set_xlabel('Frame no.')
    ax1.set_ylabel('Total Energy (kJ/mole)')
    ax12.set_ylabel('Potential Energy (kJ/mole)')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax12.get_legend_handles_labels()

    ax12.legend(lines + lines2, labels + labels2, loc=0)
    

    return fig, axes

def main(save=False, show=True):
    psf = 'output.pdb'
    traj = 'output.dcd'
    sysinfo = None
    outfile = 'data.csv'
    uni = Universe(psf, traj, outfile=outfile, sysinfo=sysinfo)

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
        rdf, g1, g2 = uni.get_rdf('name ZN', f'{seltxt}', nbins=200, range=(1,5))
        r = rdf.results['bins']
        gab = rdf.results['rdf'] 

        # Get average no. of particles within r, as function of r (Nab)
        Nab = np.cumsum(rdf.results['count']) / (g1.n_atoms * uni.trajectory.n_frames)
        
        # Plot gab and cdf vs r
        fig, ax, ax2 = plot_rdf_cdf(r, gab, Nab)
        fig.suptitle(sel)

        #pp.show()
        if save:
            fig.savefig(f"rdf_and_coordination_{sel}.png")
    

    # Distribution of number of atoms around Zn within cutoff
    cutoff = 2.1  # choose based on RDF
    atom = "ZN"
    seltxt = f"around {cutoff} name {atom}"
    counts, bin_edges = uni.get_contact(seltxt)

    fig, ax = pp.subplots()
    ax.set_xlabel("Coordination number")
    ax.set_ylabel("Frequency")
    ax.hist(counts, bin_edges)
    if save:
        fig.savefig(f"contact_around_{atom}_{cutoff}.png")

    fig,ax = pp.subplots()
    ax.set_xlabel("Frame no.")
    ax.set_ylabel("Coordination number")
    ax.plot(counts)
    if save:
        fig.savefig(f"cn_vs_frame_{atom}_{cutoff}.png")

    if show:
        pp.show()
    

if __name__ == "__main__":
    main(save=True, show=True)
