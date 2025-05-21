# Characterization of peripheral nerve recordings

This is the repository for the article *Biophysical characterization of the recording of unmyelinated and myelinated fiber activity with peripheral interfaces*. iScience (2025). DOI: [10.1016/j.isci.2025.112495](https://doi.org/10.1016/j.isci.2025.112495).

**Authors**: Claudio Verardo, Simone Romeni, Silvestro Micera

## Code overview
* `hm`: Matlab classes to instatiate hybrid models and simulate SUAPs
* `neuron`: Python functions to generate the fiber excitation templates in Neuron
* `run_simulations`: Matlab scripts to simulate the data used in the paper
* `run_analyses`: Matlab scripts to reproduce the figures of the paper

## Requirements
* Matlab R2023a
* Python 3.10.11
* Comsol Multiphysics 6.1 (with LiveLink for Matlab)
* Neuron 8.2.2

## Installation
* Add Python dependencies (from terminal)
```
pip install -r requirements.txt
```
* Compile the Neuron mechanisms (from terminal)
```
nrnivmodl neuron/fibers/mod_mrg02
nrnivmodl neuron/fibers/mod_sundt15
```
* Add paths (from Matlab): run `setup.m`

## Getting started

### Generate fiber excitation templates
```
cd neuron
python excitation_template_fibers_parallel.py internode-CV DIR_DATA N_CORES
```
* `DIR_DATA`: path to the directory where to store the simulated data
* `N_CORES`: number of CPU cores to use for the simulations

### Instantiate hybrid models and simulate SUAPs
Execute the following Matlab script(s) in `run_simulations` depending of the figure(s) data you want to reproduce. The figure panels related to each script are listed at the first lines of its code. `DIR_DATA` must be specified at the beginning of each script before executing it.

* `sim1_SUAPs.m`
* `sim2_SUAPs_along_fiber.m`
* `sim3_SUAPs_along_fiber_perineurium.m`
* `sim4_SUAPs_along_fiber_cuff_coverage.m`
* `sim5_SUAPs_along_fiber_cuff_design.m`
* `sim6_SUAPs_along_fiber_chronic.m`

NOTE: the Comsol LiveLink for Matlab must be active (`mphstart` command). Please, refer to the official documentation for further details.

### Reproduce figures of the paper
Execute the following Matlab script(s) in `run_analyses` depending of the figure plot(s) you want to reproduce. The figure panels related to each script are listed at the first lines of the script. `DIR_DATA` must be specified at the beginning of each script before executing it.

* `fig1_SUAPs_overview.m`
* `fig2_S1_SUAPs_spatial_maps.m`
* `fig3_4_S3_S4_SUAPs_fixed_distance.m`
* `fig3_4_S3_S4_SUAPs_overlap_win.m`
* `fig5_6_SUAPs_vs_distance_violin.m`
* `fig5_SUAPs_perineurium.m`
* `fig6_SUAPs_perineurium.m`
* `fig6_SUAPs_cuff_coverage.m`
* `fig6_SUAPs_cuff_design.m`
* `fig7_SUAPs_overlap_win_chronic.m`

NOTE: the SUAPs data must be simulated beforehand (see previous section).