# The Effect of Federal Delisting on US Wolf Population Dynamics

Data and R code for reproducing the analyses and figures in:

> Ramírez-Loza, J. P., López-Bao, J. V., Epstein, Y., & Chapron, G. (2026).
> **The Effect of Federal Delisting on US Wolf Population Dynamics.**
> *Conservation Letters*, 19(1). [doi:10.1111/con4.70028](https://doi.org/10.1111/con4.70028)

## Repository contents

| File or folder | Purpose |
| --- | --- |
| [all_timeseries.xlsx](all_timeseries.xlsx) | Population, monitoring, ESA listing, and harvest/culling data. |
| [data_prep.r](data_prep.r) | Prepares inputs for density-independent models and records the local repository path. |
| [data_prep_dd.r](data_prep_dd.r) | Prepares inputs for density-dependent models. |
| `model_*.r` scripts | Fit models, run diagnostics, and generate figures; see the model guide below. |
| [Figures/](Figures/) | Generated PDF figures. |
| [JAGS_Outputs/](JAGS_Outputs/) | Generated MCMC outputs, including inputs for the informative priors. |

This Markdown README is the maintained guide displayed by GitHub. Scripts and
the source workbook remain at the repository root so their relative paths and
existing output names stay consistent.

## Requirements

- R.
- JAGS, installed separately from R.
- R packages: `readxl`, `boot`, `rjags`, `coda`, `scales`, `lubridate`, and `dplyr`.

Install missing R packages from an R console:

```r
install.packages(c("readxl", "boot", "rjags", "coda", "scales", "lubridate", "dplyr"))
```

JAGS must be available for `rjags` to load. Package versions are not locked in
this repository.

## Run the analyses

### 1. Prepare the data

Start a dedicated R session with the repository root as its working directory
(the folder containing `all_timeseries.xlsx`). The scripts clear objects from
the current R environment when sourced.

```r
setwd("/path/to/US-wolf-delisting")  # Replace with your local repository path.
source("data_prep.r")
source("data_prep_dd.r")
```

`data_prep.r` uses the current working directory, creates the output folders if
needed, and writes `root_dir.rds` and `bugs_data.rds`. `data_prep_dd.r` writes
`bugs_data_DD.rds`, excluding West Washington and West Oregon.

Run all subsequent scripts from the repository root: they read `root_dir.rds`
before changing directories. If you move the repository, rerun both preparation
scripts to update the saved local path and inputs.

### 2. Fit models and generate figures

Run the scripts individually as needed, in the order below when reproducing all
analyses. Models 1–3 use `bugs_data.rds`; models 4–6 and the final three scripts
use `bugs_data_DD.rds`.

```r
# Density-independent models
source("model_1.r")
source("model_2.r")
source("model_3.r")

# Density-dependent models with full mortality compensation
source("model_4.r")
source("model_5.r")
source("model_6.r")

# Generate informative priors before fitting the corresponding additive models
source("model_7_1.r")
source("model_8.r")
source("model_7_2.r")
source("model_9.r")
source("model_7_3.r")
source("model_10.r")
```

Model fits use long MCMC runs and can take substantial time. Run them
sequentially within a checkout: every model writes the same temporary
`model.bug` file. Rerunning a model overwrites its existing outputs. Random seeds
are not explicitly fixed in the scripts, so repeated MCMC runs can differ.

### Models and figures

**Script numbering differs from paper numbering after Model 6.** The three
`model_7_*.r` scripts are auxiliary models used to estimate informative priors;
their results do not appear in the paper. Scripts `model_8.r`, `model_9.r`, and
`model_10.r` implement paper Models 7, 8, and 9, respectively.

| Script | Paper model | Description | Figures |
| --- | --- | --- | --- |
| [model_1.r](model_1.r) | 1 | General, density independent | Figure 3, both panels |
| [model_2.r](model_2.r) | 2 | Region specific, density independent | Figure S4 |
| [model_3.r](model_3.r) | 3 | State specific, density independent | Figure 4, columns 1–3 |
| [model_4.r](model_4.r) | 4 | General, density dependent, full harvest/culling mortality compensation | Figure S3a |
| [model_5.r](model_5.r) | 5 | Region specific, density dependent, full harvest/culling mortality compensation | Figure S5a |
| [model_6.r](model_6.r) | 6 | State specific, density dependent, full harvest/culling mortality compensation | Figure S6 |
| [model_7_1.r](model_7_1.r) | Auxiliary | General model incorporating harvest/culling; priors for paper Model 7 | — |
| [model_7_2.r](model_7_2.r) | Auxiliary | Region-specific model incorporating harvest/culling; priors for paper Model 8 | — |
| [model_7_3.r](model_7_3.r) | Auxiliary | State-specific model incorporating harvest/culling; priors for paper Model 9 | — |
| [model_8.r](model_8.r) | 7 | General, density dependent, full harvest mortality additivity | Figure S3b |
| [model_9.r](model_9.r) | 8 | Region specific, density dependent, full harvest mortality additivity | Figure S5b |
| [model_10.r](model_10.r) | 9 | State specific, density dependent, full harvest mortality additivity | Figure 4, column 4 |

Each model script saves `*_coda.RData` and `*_jags.RData` files under
`JAGS_Outputs/`. These outputs retain their original filenames so existing
results and informative priors remain compatible. For example, `model_1.r`
writes `1.Model1_coda.RData` and `1.Model1_jags.RData`. PDF filenames are also
unchanged. Use the new lowercase script names in `source()` and `Rscript`
commands.

The prior dependencies are:

| Run first | Required output in `JAGS_Outputs/` | Used by |
| --- | --- | --- |
| `model_7_1.r` | `7.1Model7_coda.RData` | `model_8.r` |
| `model_7_2.r` | `7.2Model7_coda.RData` | `model_9.r` |
| `model_7_3.r` | `7.3Model7_coda.RData` | `model_10.r` |

PDF files are saved under `Figures/`:

| Script | PDF filenames | Figure placement |
| --- | --- | --- |
| `model_1.r` | `M1_Beta1.pdf`, `M1_Preds.pdf` | Figure 3 |
| `model_2.r` | `M2_beta1_<Region>.pdf`, `M2_predGR_<Region>.pdf` | Figure S4 |
| `model_3.r` | `M3_tseries_<STATECODE>.pdf` | Figure 4, column 1 |
| `model_3.r` | `M3_fit_<State>.pdf` | Figure 4, column 2 |
| `model_3.r` | `M3_predGR_<State>.pdf`, `M3_beta1_<State>.pdf` | Figure 4, column 3 and its subpanels |
| `model_4.r` | `M4_Beta1.pdf` | Figure S3a |
| `model_5.r` | `M5_beta1_Great Lakes.pdf`, `M5_beta1_Northern Rocky Mountains.pdf` | Figure S5a |
| `model_6.r` | `M6_predGR_<State>.pdf`, `M6_beta1_<State>.pdf` | Figure S6 |
| `model_8.r` | `M8_Beta1.pdf` | Figure S3b |
| `model_9.r` | `M9_beta1_<Region>.pdf` | Figure S5b |
| `model_10.r` | `M10_predGR_<State>.pdf`, `M10_beta1_<State>.pdf` | Figure 4, column 4 |

`<State>` and `<Region>` are placeholders for the names used by the scripts;
`<STATECODE>` is one of the workbook sheet codes below.

## Data dictionary

The workbook contains one sheet per state or state segment:

| Sheet | State or segment |
| --- | --- |
| `ID` | Idaho |
| `MT` | Montana |
| `WY` | Wyoming |
| `WA_NRM` | Washington, Northern Rocky Mountains |
| `OR_NRM` | Oregon, Northern Rocky Mountains |
| `WA_W` | West Washington |
| `OR_W` | West Oregon |
| `WI` | Wisconsin |
| `MI` | Michigan |
| `MN` | Minnesota |

Each sheet contains these columns (names are case sensitive):

| Column | Description |
| --- | --- |
| `Year` | Year of the population estimate. |
| `Popsize` | Estimated population size for that year. |
| `MP_Start` | Start date of the monitoring/reporting period to which the estimate applies. |
| `MP_End` | End date of that monitoring/reporting period. |
| `MP_Total` | Number of days in the monitoring/reporting period. |
| `Date_Listed` | Date of listing under the Endangered Species Act (ESA) during the period, if applicable. |
| `Date_Delisted` | Date of delisting under the ESA during the period, if applicable. |
| `Days_Listed` | Number of days in the period during which wolves were listed under the ESA. |
| `Days_Delisted` | Number of days in the period during which wolves were not listed under the ESA. |
| `N_Killed` | Number of wolves killed in legal hunting or culling between `MP_Start` and `MP_End`. |
