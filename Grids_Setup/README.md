The procedure for producing interpolation grids using NNLOJET and FastNLO is outlined below.

The directory structure where the grids will be produced should follow this layout:

```
main folder
 ├── LO
 ├── R
 ├── V
 ├── RRa
 ├── RRb
 ├── RV
 ├── VV
 └── combine
```
Here, LO, R, V, RRa, RRb, RV, and VV are folders containing results corresponding to different contributions to the final result.

This repository provides four examples in the following directories:

 - ATLAS8TeV/Wplus

 - ATLAS8TeV/Wminus

 - ATLAS13TeV/Wplus

 - ATLAS13TeV/Wminus

These correspond to two different setups (ATLAS8TeV and ATLAS13TeV) and the production of W⁺ and W⁻ bosons in each case.

Inside each contribution folder, the process begins with preparing a runcard (`.run` file), which specifies the process, observables, binning, scale choice, and other relevant parameters for the NNLOJET calculation.

The first computational step is the VEGAS warm-up, which explores the phase space and optimises the integration grid. This step can be run in a multithreaded environment to minimise computing time.

The VEGAS warm-up is performed by running NNLOJET in warm-up mode. In the runcard (with the production line commented out), this looks as follows:

```
RUN  WpJ_LO
  PDF = NNPDF31_nnlo_as_0118[0]
  tcut = 1d-7
  iseed   = 1
  imethod = 2
  iplot   = 0
  angular_average  = .true.
  print_max_weight = .false.
  cache_kinematics = .false.
  pole_check       = .false.
  scale_coefficients = .true.
  warmup = 10000000[6]
!  production = 50000000[1]
END_RUN
```

Once the VEGAS warm-up is complete and the corresponding files are generated, a FastNLO warm-up is performed to initialise the interpolation grid structure. This FastNLO warm-up run is automatically executed if `.wrm` files are not found in the directory when running NNLOJET in production mode and the `grid` argument is specified in the `HISTOGRAMS` section of the runcard. For example:

```
RUN  WmJ_LO
  PDF = NNPDF31_nnlo_as_0118[0]
  tcut = 1d-7
  iseed   = 1
  imethod = 2
  iplot   = 0
  angular_average  = .true.
  print_max_weight = .false.
  cache_kinematics = .false.
  pole_check       = .false.
  scale_coefficients = .true.
!  warmup = 10000000[6]
  production = 50000000[1]
END_RUN
```

```
HISTOGRAMS
	ht_full > ht_full_var [50,100,150,200,250,300,350,400,450,500,550,600,650,700,800,900,1000,1100,1200,1400,1600,2500] grid=ht_full.fast
	ptw > ptw_var [0,25,50,75,100,125,150,175,200,250,300,350,400,450,500,600,800] grid=ptw.fast
	ptj1 > ptj1_var [30,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,450,500,550,600,700,1000] grid=ptj1.fast
	abs_yj1 > abs_yj1_var [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4] grid=abs_yj1.fast
END_HISTOGRAMS
```


During the FastNLO warmup, events are generated without integration in order to evaluate their distribution across the interpolation dimensions: x1 , x2 , and $\mu$. Given the specified number of interpolation nodes per dimension and bin, FastNLO applies a procedure to optimally allocate them. The details on the grid configuration are specified in the `fastnlo.str` file. Since the grid structure must be identical for all contributions of the calculation, the `.wrm` files must be produced only once, and copied to every contribution folder (LO, R, V, RRa, RRb, RV and VV). 

The full event generation and integration (production run) can then be submitted to a computing cluster, where the different cross section contributions (LO, R,
V, RRa , RRb , RV, and VV) are evaluated. The integration is parallelised by submitting certain number of randomly initialised jobs, whose results are subsequently combined. The submission to the cluster can be handled using the scripts `submit.sh` and `job_node.sh`, assuming the cluster is managed by SLURM.

The results are combined in two stages within the `combine` folder:

 - Intra-contribution combination:
    Run `combine_all_1.sh` to combine results from different jobs within the same contribution folder.

 - Final combination:
    Run `combine_all_2.sh` to combine results from all contributions, producing the final result.

Some configuration settings for the combination process are specified in the `combine.ini` file. After this step, both the "exact" NNLOJET results computed using the PDF set specified in the runcard (here NNPDF31_nnlo_as_0118) and the corresponding interpolation grids are available for each requested histogram.

The resulting grids can then be evaluated using the same PDF set that was used for the "exact" NNLOJET prediction: NNPDF31_nnlo_as_0118. This can be done by running the `read_tab_gz.sh` script. The comparison between the grid-based prediction and the corresponding result from NNLOJET constitutes the closure test of the grid, through which the accuracy and reliability of the interpolation grid is assessed.

The script `read_tab_gz.sh` also re-evaluates the grids with scale variations by factors of 2 and 1/2 relative to the central scale choice.
This allows for the estimation of theoretical (scale) uncertainties.

The `analysis.py` and `analysis_thesis.py` scripts allow to plot closure tests and comparison with experimental data including scale uncertainties.

After validation, the grids are ready to be used in PDF fits.