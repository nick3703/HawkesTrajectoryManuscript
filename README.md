# HawkesTrajectoryManuscript

R code to conduct the simulations and fit the results reported in *Hawkes
Trajectories: A Framework for Leadership Inference in Coordinated Movement*, by
Clark and Mateu.

`HawkesTrajAnalysis.R` regenerates every numbered table and figure in the
manuscript. Each block of the script is named for what it produces.

## Data

The pigeon flight data are available from the Dryad Digital Repository
(Watts et al. 2016) at <https://doi.org/10.5061/dryad.508j2>. Place the
unzipped CSV files in a `Data/` directory alongside the script.

## Requirements

R (>= 4.0) with the `data.table` and `sf` packages.

```r
install.packages(c("data.table", "sf"))
```

## Usage

```bash
Rscript HawkesTrajAnalysis.R                 # regenerate everything
Rscript HawkesTrajAnalysis.R Table9          # a single block
Rscript HawkesTrajAnalysis.R Table9 Figure1  # several blocks
```

Results are written to CSV in the working directory, and the figure to PNG.

## Blocks

| Block | Produces | Description |
|---|---|---|
| `Tables1to6` | Tables 1–6 | Simulation studies: recovery with fixed kernels, estimated kernels, robustness, and head-to-head comparison against directional delay |
| `Table7` | Table 7 | Bin width and working likelihood |
| `Table8` | Table 8 | Sandwich coverage under a binary response |
| `Sec6Edges` | Section 6 | Threshold-free edge selection: AUC, average precision, bootstrap rule |
| `Sec6Driver` | Section 6 | Shared group-level driver simulation |
| `Sec6Calib` | Section 6 | Calibrated coupling sweep |
| `Sec7Fits` | Section 7 | Intensity decomposition and the distribution of ρ(R) |
| `Sec7Boot` | Section 7 | Cluster bootstrap intervals |
| `Sec7Resid` | Section 7 | Residual diagnostics |
| `Sec7Lag` | Section 7 | The 2 s oscillation and its shared component |
| `Sec7Origin` | Section 7 | Raw-signal spectrum before binning |
| `Table9` | Table 9 | Null battery |
| `Sec7Blocks` | Section 7 | Block-restricted timing nulls |
| `Sec7Thresh` | Section 7 | Null battery across event thresholds |
| `Sec7Events` | Section 7 | Alternative event definitions |
| `Sec7RhoACF` | Section 7 | ρ(R) against residual structure |
| `Table10` | Table 10 | Group-activity decomposition |
| `Sec7FlockNull` | Section 7 | Null battery under the group-activity model |
| `Figure1` | Figure 1 | Representative flight |

## Reproducibility

`SEED_BASE` is set once at the top of the script and each block adds a fixed
offset, so any block reproduces identically whether it is run on its own or as
part of a full run. Setting the seed only once would not have this property,
since the state on entering a block would depend on which blocks preceded it.

The null-battery blocks (`Table9`, `Sec7Blocks`, `Sec7Thresh`,
`Sec7FlockNull`) are checkpointed per flock-release in a `ckpt/` directory and
resume if interrupted. These take several hours each; the remaining blocks run
in minutes.

## Reference

Watts, I., Nagy, M., Burt de Perera, T. and Biro, D. (2016). Misinformed
leaders lose influence over pigeon flocks. *Biology Letters* 12(9), 20160544.
<https://doi.org/10.1098/rsbl.2016.0544>
