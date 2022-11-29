# Repo for 'Joint distribution of wave height, wave period, and wind speed'

This is the repository for the Chapter "Joint distribution of wave height, wave
period, and wind speed" in the doctoral thesis of Andreas Haselsteiner. https://doi.org/10.26092/elib/1615

The chapter is based on the paper paper 'Global hierarchical models for wind and
wave contours: Physical interpretations of the dependence functions' by A. F.
Haselsteiner, A. Sander, J.-H. Ohlendorf and K.-D. Thoben. https://doi.org/10.1115/OMAE2020-18668

The results presented in this paper were also submitted to the benchmarking
exercise on estimating extreme environmental conditions (https://doi.org/10.1016/j.oceaneng.2021.109504).

To reproduce the results from the OMAE 2020 paper you can navigate to the commit from March 25, 2020:
https://github.com/ahaselsteiner/2020-paper-omae-hierarchical-models/commit/e5781fcec447fbdd30f1248f6d834342326bb9ae

At that time it contained the following things:
 * The used datasets (folder 'datasets')
 * Python files to reproduce the analysis
   * compute_contours_datasets_abc.py for the sea state models
   * compute_contours_datasets_def.py for the wind-wave state models
 * The coordinates of the computed environmental contours (folder 'contour-coordinates')

The newest version also contains scripts to produce additional figure that are presented in the PhD thesis.

## Download and use the repository
To download this repository type
```console
git clone https://github.com/ahaselsteiner/2020-paper-omae-hierarchical-models.git
```

Then install the required Python packages by typing
```console
pip install -r requirements.txt
```
(if you are using pip for installing packages)
