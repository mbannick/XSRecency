[![DOI](https://zenodo.org/badge/327072490.svg)](https://zenodo.org/badge/latestdoi/327072490)

# Cross-Sectional Incidence Estimation of HIV Using Recent Infection Testing Algorithms

## Authors
[Marlena Bannick](https://marlenabannick.com/) and [Fei Gao](https://www.fredhutch.org/en/faculty-lab-directory/gao-fei.html)

## Description

For background on cross-sectional incidence estimators, please see [Gao and Bannick (2022)](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9296).

### Estimation of Cross-Sectional Incidence
This package provides functions to estimate cross-sectional incidence of HIV using recent infection testing algorithms (RITA). Available estimators are:

- Snapshot estimator ([Kaplan and Brookmeyer 1999](https://doi.org/10.1287/opre.47.1.29))
- Adjusted estimator ([Kassanjee et al. 2012](https://doi.org/10.1097/EDE.0b013e3182576c07))
- Enhanced estimator (Bannick et al. 2023+)

### Estimate Properties of Recent Infection Testing Algorithms
We include functions to interface with the CEPHIA Public Use Dataset ([Grebe et al. 2021](https://doi.org/10.5281/zenodo.4900634)) which provides data on recent infection testing algorithms. The functions in this package use CEPHIA data to estimate properties of recent infection testing algorithms, like the mean duration of recent infection (MDRI). They are also crucial for the enhanced estimator, as the enhanced estimator requires an estimate of the entire test-recent function for a RITA (which is not typically available in the literature).

### Simulate Data
We provide functions to simulate cross-sectional data samples based on different epidemiological dynamics. These are useful for demonstrating how to use cross-sectional incidence estimation techniques and checking the validity of estimators.

## Installation

You can install the package with `devtools`:
```{R}
require(devtools)
require(formula.tools)
require(zen4R)
require(geepack)
devtools::install_github("mbannick/XSRecency")
```

Or to download the package, you may clone the repository:
```{R}
git clone https://github.com/mbannick/XSReceny.git
```

## References

### Related Publications

- Bannick, M. S., Donnell, D., Hayes, R., Laeyendecker, O., Gao F. (2023+). Enhanced Cross-Sectional HIV Incidence Estimators that Incorporate Prior HIV Test Results.
- Gao, F., & Bannick, M.S. (2022). Statistical considerations for cross-sectional HIV incidence estimation based on recency test. *Statistics in Medicine*. 41(8): 1446–1461. [doi:10.1002/sim.9296](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9296)

To reproduce the simulations in Bannick et al. 2023, please see [this repository](https://github.com/mbannick/RITA-plus-sims), which
uses functions from this package.

To reproduce the simulations in [Gao and Bannick (2022)](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9296), please download a previous release (0.1.0) of this repository [here](https://github.com/mbannick/XSRecency/releases/tag/0.1.0).

### Additional Materials

- Kaplan, E. H., & Brookmeyer, R. (1999). Snapshot Estimators of Recent HIV Incidence Rates. *Operations Research*. 47(1): 29–37. [doi:10.1287/opre.47.1.29](https://doi.org/10.1287/opre.47.1.29)
- Kassanjee, R., McWalter, T. A., Bärnighausen, T., & Welte, A. (2012). A New General Biomarker-based Incidence Estimator. *Epidemiology*. 23(5): 721–728. [doi:10.1097/EDE.0b013e3182576c07](https://doi.org/10.1097/EDE.0b013e3182576c07)
- Grebe, E., et al. (2021). CEPHIA public use data (1.0) Data set. Zenodo. [doi:10.5281/zenodo.4900634](https://doi.org/10.5281/zenodo.4900634)
