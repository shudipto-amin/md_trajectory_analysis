# README 

This is an example of how to use MDAnalysis to analyse protein trajectories - 
calculating things like RDF, coordination numbers, average angles, etc. 

## How it works

The main script here is [analyse.py](https://github.com/shudipto-amin/md_trajectory_analysis/blob/master/analyse.py). 
It defines a class `Universe` which works the same way as 
[MDAnalysis.Universe](https://userguide.mdanalysis.org/stable/universe.html),
but with added functionality.

This script is meant to be used as like a module,
i.e., it is meant to be imported in another script or interactive shell:

```
import analyse as ana
uni = ana.Universe(<topology_file>, <coordinate_file>, *args, **kwargs )
```

Now you have a `Universe` object
which inherits everything from `MDAnalysis.Universe`
and adds some extra functionality, like centering the protein,
getting angles and distances, etc.

## Centering the protein

One important feature implemented in `analyse.Universe`
is the function to center the protein within the box
and write it to a new file. 

For example,
```
uni = ana.Univese('output.pdb', 'output.dcd')
uni.center_protein_write('centered_output.dcd')
```
will  center the protein and save it to `centered_output.dcd`.

> __It is really important that you center the protein first,
> then load the new centered trajectory for further analysis.__
> Otherwise, the distance based calculations (RDF, RMSD, etc) might 
> not be right.

## Example: Calculating angles and distances with `analyse.py` as a module.

[angles_distances.py](https://github.com/shudipto-amin/md_trajectory_analysis/blob/master/angles_distances.py)
shows how to import `analyse.py` as a module 
and calculate angles and distances of particular atom selections per frame.
This data is then plotted using [matplotlib](https://matplotlib.org/).

## Running `analyse.py` directly.

If importing `analyse.py` you really don't need anything else here, it is stand alone. 
However, running it directly,
```
python analyse.py
```
will require at least the `output.pdb` and `output.pdb`,
unless they are changed within the script. 
Running it directly will also plot and save some RDF, RMSD, and Coordination
figures, as done within the `main` function of `analyse.py`.

