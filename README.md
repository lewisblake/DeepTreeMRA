# Deep Tree MRA - A Parallel Implementation of the Multi-Resolution Approximation for High Performance Computing Environments

### Software authors:
*	Lewis Blake, Colorado School of Mines (lblake@mines.edu).
*	Dorit Hammerling, Colorado School of Mines (hammerling@mines.edu).

This code is based on the MRA model described in
["A multi-resolution approximation for massive spatial datasets" by Matthias Katzfuss, 2017](https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2015.1123632) in the Journal of American Statistical Association (DOI: 10.1080/01621459.2015.1123632). Also available at [arXiv](https://arxiv.org/abs/1507.04789).
Throughout the codebase there are references to this manuscript.

Designed and implemented with Matlab 2020b
Previous versions may not be supported.
Required toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox

## Getting Started

This Matlab codebase allows users to apply the multi-resolution approximation model to 2D spatial data sets.

The `user_input.m` script is where much of the user input can be modified (see [User Input](#user_input) and [Additional User Input](#additional_user_input) below).
If desired, users should modify model parameters within `user_input.m`.
The `main.m` script runs the model. 
Within the Matlab Editor Tab, selecting the 'Run' button from `main.m` will execute the code.


The repository is structured as follows: 
The repositiory contains all files required to run the code.
Within the `Distributed MRA` folder, there are four other folders: `Data`, `Plots`, `Results`, and `subroutines`.
The `Data` folder contains example data sets. 
The `Results` folder is the default folder for results to be saved. Initially empty. 
The `Plots` folder is the default folder for spatial prediction plots to be saved. Initially empty.
The `subroutines` folder contains the functions used to execute the model which are called from `main.m`.

## Example Data

Two datasets are included in this distribution: `satelliteData.mat` and `simulatedData.mat`. 
These files are contained within the `Data` folder.
Both of these data sets were originally used in [Heaton, M.J., Datta, A., Finley, A.O. et al. JABES (2018).](https://doi.org/10.1007/s13253-018-00348-w)


## <a name = "parallelization"></a> Paralleization:

A significant benefit of the MRA is that it lends itself to execution in parallel. 
For this reason, the portions of the creation of the prior, portions of the posterior inference, and spatial prediction in this codebase were designed to run in parallel using `spmd. 
Within `user_input.m`, users can specify the number of workers in the parallel pool by setting `NUM_WORKERS`. 
SPMD blocks throughout the code will execute in parallel using `NUM_WORKERS` workers.
Within the Matlab Cluster Profile Manager, the user can specify the desired cluster settings.
These settings include the number of nodes, nmber of cores, number of MPI processes, wall times, and many others.
For further reference please see https://www.mathworks.com/help/distcomp/discover-clusters-and-use-cluster-profiles.html .

## Preliminaries:

### `find_num_levels_suggested.m` 

This is stand-alone function estimates `NUM_LEVELS_M` for a given dataset size (`NUM_DATA_POINTS_n`), number of knots (`NUM_KNOTS_r`), and number of partitions (`NUM_PARTITIONS_J`).
Note that `NUM_LEVELS_M` is a positive integer.

## <a name="user_input"></a> User Input

### user_input.m 

In user_input.m, the areas requiring user input are as follows:

dataSource: | 'satellite' | 'simulated' |
    - These are the dataSource's for the data provided. Default is 'satellite'.
    - In order to use a different data set, see the section of load_data.m below and feed the case string for the data used to dataSource in user_input.m.

calculationType: | 'prediction' | 'optimize' | 'likelihood' | 'build_structure' |
Default is 'likelihood'.
calculationType can be set to any of the following calculation modes:
	- prediction: Uses given values for the parameters (theta and varEps) and just conducts spatial prediction. Parameters can be changed in load_data.m	
	- optimize: Optimizes over the range, variance and measurement error. The range and variance parameters are stored as a vector: theta. The measurement error is stored as a double: varEps.	
	- likelihood: Calculates the log-likelihood.
	- build_structure: Builds the multi-resolution structure. Reports summary statistics and produces a histogram of the number of observations assigned to regions at the finest resolution.

#### User Input relevant for any calculationType:

`NUM_LEVELS_M`: Total number of levels in the hierarchical domain-partitioning. By default set to 9.

`NUM_PARTITIONS_J`: Number of partitions for each region at each level. Only implemented for J = 2 or J = 4. By default set to 2.

`NUM_KNOTS_r`: Number of knots per partition. By default set to 64.

* `offsetPercentage`: Offset percentage from partition boundaries. Must be between 0 and 1.
This quantity determines the buffer between the boundaries of a region where knots can be placed.
`offsetPercentage` is also used at the coarsest resolution to extend the maximal x and y domain boundaries as to include data points that may be exactly on the boundary within a region. The domain boundaries define a rectangular region determined by the minimal and maximal x and y coordinate locations. Preferably set `offsetPercentage` to be a small number (e.g. 0.01).

* `NUM_WORKERS`: Number of workers in the parallel pool. Must be set to be a power of J <= nRegions(nLevelsInSerial). See [Parallelization](#parallelization) above. Default is 4.

* `NUM_LEVEL_ASSIGN_REGIONS_P`: The level at which regions are assigned across workers. This determines how much for the hierarchical domain paritioning each worker is assigned. Default is 3.

* `verbose`: Boolean variable indicating whether to produce progress indicators. Default is true.

* `resultsFilePath`: Optional file path to save results for each `calculationType`. 
Set to be a string (e.g. `resultsFilesPath = '/Users/JerryGarcia/Desktop/';`). By default results are saved in the `Results` folder.

#### User inputs relevant if `calculationType = "prediction"` or `"build_structure"`

* `displayPlots`: Boolean variable indicating whether to display plots if predicting or building the multi-resolution structure.
* `savePlots`: Boolean variable indicating whether to save plots if predicting or building the multi-resolution structure.
(Note: If not executing the `"prediction"` or `"build_structure"` `calculationType`, these booleans are not relevant.)

* `nXGrid`: Number of prediction grid points in x-direction. By default set to 200.
* `nYGrid`: Number of prediction gridpoints in y-direction. By default set to 200.
(Note: These parameters define a `nXGrid` x `nYGrid` prediction grid of spatial prediction locations if predicting.
The prediction grid is only defined within rectangular region given by the domain boundaries discussed above.)

* `plotsFilePath`: Optional file path to save prediction plots if plotting.
Set to be a string (e.g.`plotsFilesPath = '/Users/JerryGarcia/Pictures/';`).
By default plots are saved in the `Plots` folder.

#### User inputs relevant if calculationType = 'optimize'

* `lowerBound`: Vector of lower-bound values required for the optimization search. Default is [0,0,0].

* `upperBound`: Vector of upper-bound values required for the optimization search. Default is [10,1,5].

* `initialEstimate`: Vector of inital estimates of parameteres required for the optimization search. Default is [5, 0.3, 0.1].



## <a name = "additional_user_input"></a> Additional User Input



