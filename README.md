[![DOI](https://zenodo.org/badge/78135650.svg)](https://zenodo.org/badge/latestdoi/78135650)

This is the code used in the [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0170164):

>Jacobs C, Lambourne L, Xia Y, Segrè D. Upon accounting for the impact of isoenzyme loss, gene deletion costs anticorrelate with their evolutionary rates. PLoS ONE. 2017 Jan. 12(1): e0170164.

### Reproducing figures

There are four scripts to run to generate the FBA model predictions for the cost of knocking
out genes, in different media, with different
assumptions in the cost calculations.

`download_and_process_data.py`: Run once to obtain the input data.

`send_*_knockouts_to_cluster.py`: Calculates growth for all single/double gene knockouts using the Open Grid Scheduler batch system on a computing cluster.

`collate_cluster_output.py`: Run once after all batch jobs have finished to combine their outputs.

The figures related to the single knockouts are produced in the jupyter notebook `examine_correlations.ipynb` and the figures related to the double knockouts are produced in `epistasis.ipynb`.

### Dependancies

Uses Python 2.7. Uses the [cobrapy](opencobra.github.io/cobrapy) package. The FBA requires a linear programming solver.
I've used [Gurobi](www.gurobi.com), which requires signing up for an
account on their website before you can install the
software.
The other python package dependancies are listed in `requirements.txt`.

### Using other FBA models

The code is written to analyze the [yeast 7.6](https://sourceforge.net/projects/yeast/files/) FBA model. Since there are no standardized reaction IDs across different models, you need to edit the `carbon_sources.txt` and `nitrogen_sources.txt` files by hand to have the correct names for the exchange reactions for the new model.

### Boston University cluster

To get gurobi working in the scc1 cluster:

```
module load gurobi
export PYTHONPATH=/share/pkg/gurobi/6.5.2/install/python/python2.7/site-packages:$PYTHONPATH
```
