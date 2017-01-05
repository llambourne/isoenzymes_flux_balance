This is the code used in the paper:
**Upon accounting for the impact of isoenzyme loss, gene deletion costs anticorrelate with their evolutionary rates**, published in PLOS ONE.

### Reproducing figures

There are four scripts to run to generate the FBA model predictions for the cost of knocking
out genes, in different media, with different
assumptions in the cost calculations.

`download_and_process_data.py`: Run once to obtain the input data.

`send_*_knockouts_to_cluster.py`: Calculates growth for all single/double gene knockouts using the Open Grid Scheduler batch system on a computing cluster.

`collate_cluster_output.py`

The figures related to the single knockouts are produced in the jupyter notebook `examine_correlations.ipynb` and the figures related to the double knockouts are produced in `epistasis.ipynb`.

### Dependancies

The FBA requires a linear programming solver.
I've used gurobi, which requires signing up for an
account on their website before you can install the
software.
The other python package dependancies are listed in `requirements.txt`.

### Boston University cluster

To get gurobi working in the scc1 cluster:

```
module load gurobi
export PYTHONPATH=/share/pkg/gurobi/6.5.2/install/python/python2.7/site-packages:$PYTHONPATH
```
