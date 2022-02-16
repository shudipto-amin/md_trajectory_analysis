import MDAnalysis as mda  # make sure you have latest version > 2.0.0
import numpy as np
import matplotlib.pyplot as pp

from MDAnalysis.analysis.rms import RMSD

# Create two universes for each trajectory
uni1 = mda.Universe('output.pdb', 'output1.dcd')
uni2 = mda.Universe('output.pdb', 'output2.dcd')

# Create two rmsd objects for each trajectory
rmsd1 = RMSD(uni1, select='protein') 
rmsd2 = RMSD(uni2, select='protein', reference=uni1, ref_frame=0)

# Run first one on first 100, second one on last 100
rmsd1.run(stop=100)
rmsd2.run(start = uni2.trajectory.n_frames - 100)

# Combine the two separate rmsds into one array
rmsd_only1 = rmsd1.results['rmsd'][:,2]
rmsd_only2 = rmsd2.results['rmsd'][:,2]

rmsd_ext = np.concatenate([rmsd_only1, rmsd_only2])

# Plot
fig, ax = pp.subplots()
ax.plot(rmsd_ext)
ax.set_xlabel("Combined frame no.")
ax.set_ylabel("RMSD of protein ($\AA$)")
pp.savefig("combined_rmsd.pdf")
