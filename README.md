# HawkesTrajectoryManuscript

R Code to conduct simulations and fit results from _Hawkes Trajectories: A Framework for Leadership Inference in Coordinated Movement_ by Clark and Mateu.

The file contains code to fit all results given in the manuscript. It requires the pigeon flight data available in the Dryad Digital Repository (Watts et al. 2016) at  https://doi.org/10.5061/dryad.508j2

Each block is named for what it produces in the manuscript.

   Tables1to6      Tables 1-6   simulation studies: recovery with fixed
                                 kernels, estimated kernels, robustness,
                                 head-to-head against directional delay
      Table7          Table 7      bin width and working likelihood
      Table8          Table 8      sandwich coverage under a binary response
      Sec6Edges       Section 6    threshold-free edge selection, AUC and
                                   average precision, bootstrap rule
      Sec6Driver      Section 6    shared group-level driver simulation
      Sec6Calib       Section 6    calibrated coupling sweep
      Sec7Fits        Section 7    intensity decomposition, rho(R) distribution
      Sec7Boot        Section 7    cluster bootstrap intervals
      Sec7Resid       Section 7    residual diagnostics
      Sec7Lag         Section 7    2 s oscillation and its shared component
      Sec7Origin      Section 7    raw-signal spectrum before binning
      Table9          Table 9      null battery
      Sec7Blocks      Section 7    block-restricted timing nulls
      Sec7Thresh      Section 7    null battery across event thresholds
      Sec7Events      Section 7    alternative event definitions
      Sec7RhoACF      Section 7    rho(R) against residual structure
      Table10         Table 10     group-activity decomposition
      Sec7FlockNull   Section 7    null battery under the group-activity model
      Figure1         Figure 1     representative flight
 
    USAGE
      Rscript HawkesTrajAnalysis.R            # everything
      Rscript HawkesTrajAnalysis.R Table9     # one block
 
    SEEDS. SEED_BASE is set once below and each block adds a fixed offset, so
    any block reproduces identically whether run alone or as part of a full
    run. 
