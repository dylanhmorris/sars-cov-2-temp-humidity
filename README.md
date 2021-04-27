# [Mechanistic theory predicts the effects of temperature and humidity on inactivation of SARS-CoV-2 and other enveloped viruses](https://doi.org/10.7554/eLife.65902)
[Dylan H. Morris](https://dylanhmorris.com)(1\*), Kwe Claude Yinda (2\*), Amandine Gamble(3\*), Fernando W. Rossine(1), Qishen Huang(4), Trenton Bushmaker(2, 5), M. Jeremiah Matson(2, 6), Neeltje van Doremalen(2), Peter J. Vikesland(4), Linsey C. Marr(4), Vincent J. Munster(2), James O. Lloyd-Smith(3)

\* These authors contributed equally

1. Dept. of Ecology and Evolutionary Biology, Princeton University, Princeton, NJ, USA
2. Laboratory of Virology, Division of Intramural Research, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Hamilton, MT, USA
3. Dept. of Ecology and Evolutionary Biology, University of California, Los Angeles, Los Angeles, CA, USA
4. Department of Civil and Environmental Engineering, Virginia Tech, Blacksburg, VA, USA
5. Department of Microbiology and Immunology, Montana State University, Bozeman, MT, USA
6. Joan C. Edwards School of Medicine, Marshall University, Huntington, WV, USA


## Repository information
This repository accompanies the article "Mechanistic theory predicts the effects of temperature and humidity on inactivation of SARS-CoV-2 and other enveloped viruses" (DH Morris et al, *eLife*, 2021). It provides code for reproducing data analysis and theoretical modeling from the paper and recreating display figures.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work as below:

- DH Morris, KC Yinda, A Gamble et al. Mechanistic theory predicts the effects of temperature and humidity on inactivation of SARS-CoV-2 and other enveloped viruses. *eLife*. 2021-04-27. DOI:[10.7554/eLife.65902](https://doi.org/10.7554/eLife.65902).

Bibtex record:
```
@article {10.7554/eLife.65902,
article_type = {journal},
title = {Mechanistic theory predicts the effects of temperature and humidity on inactivation of SARS-CoV-2 and other enveloped viruses},
author = {Morris, Dylan H and Yinda, Kwe Claude and Gamble, Amandine and Rossine, Fernando W and Huang, Qishen and Bushmaker, Trenton and Fischer, Robert J and Matson, M Jeremiah and Van Doremalen, Neeltje and Vikesland, Peter J and Marr, Linsey C and Munster, Vincent J and Lloyd-Smith, James O},
editor = {Ogbunugafor, C. Brandon},
volume = 10,
year = 2021,
month = {04},
pub_date = {2021-04-27},
pages = {e65902},
citation = {eLife 2021;10:e65902},
doi = {10.7554/eLife.65902},
journal = {eLife},
}
```

## Article abstract 
Ambient temperature and humidity strongly affect inactivation rates of enveloped viruses, but a mechanistic, quantitative theory of these effects has been elusive. We measure the stability of SARS-CoV-2 on an inert surface at nine temperature and humidity conditions and develop a mechanistic model to explain and predict how temperature and humidity alter virus inactivation. We find SARS-CoV-2 survives longest at low temperatures and extreme relative humidities (RH); median estimated virus half-life is >24 hours at 10C and 40% RH, but ~1.5 hours at 27C and 65% RH. Our mechanistic model uses fundamental chemistry to explain why inactivation rate increases with increased temperature and shows a U-shaped dependence on RH. The model accurately predicts existing measurements of five different human coronaviruses, suggesting that shared mechanisms may affect stability for many viruses. The results indicate scenarios of high transmission risk, point to mitigation strategies, and advance the mechanistic study of virus transmission.


## Directories
- ``src``: all code, including data preprocessing, Bayesian model definition and fitting, and results post-processing and figure generation:
    - ``src/cleaning``: scripts for cleaning raw data and preparing it for analysis
    - ``src/figures``: scripts for generating figures
    - ``src/fitting``: scripts for fitting models to data
    - ``src/parameters``: specification of model parameters
    - ``src/stan``: Stan source files for Bayesian inference MCMC
    - ``src/tables``: scripts for autogenerating tables of results
- ``dat``: data files in comma-separated values (``.csv``) formats
    - ``dat/raw``: raw data files
    - ``dat/cleaned``: data files processed and prepared for model fitting
- ``out``: output files
    - ``out/mcmc_chains``: Markov Chain Monte Carlo (MCMC) output, as serialized R data (``.Rds``) files. 
    - ``out/figures``: figures generated from results
    - ``out/tables``: tables generated from results
    - ``out/chain_diagnostics.csv``: diagnostic tests for MCMC convergence.
- ``virusenv``: R package containing useful functions and styling called in the scripts used

## Reproducing analysis

A guide to reproducing the analysis from the paper follows. If you encounter issues, see the **Troubleshooting** section at the end of this README.

### Getting the code
First download this repository. The recommended way is to ``git clone`` it from the command line:

    git clone https://github.com/dylanhmorris/sars-cov-2-temp-humidity.git

Downloading it manually via Github's download button should also work.

### Dependency installation
The analysis can be auto-run from the project ``Makefile``, but you may need to install some external dependencies first. See the **Dependency installation guide** below for a complete walkthrough. In the first instance, you'll need a working installation of the statistical programming language R, a working C++ compiler, and a working installation of Gnu Make or similar. A few external R packages can then be installed from the command line by typing.

    make depend

from within the project directory. This will install external R packages as well as the project package ``virusenv``.

### Running the analysis

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures, tables, MCMC output (saved as R Dataset ``.Rds`` files in the ``out/mcmc-chains/`` directory as ``<model_name>_chains.Rds``), etc. The MCMC output can be loaded in any working R installation, as long as the package ``rstan`` is also installed.

If you want to do things piecewise, typing ``make <filename>`` for any of the files listed in the ``dat/cleaned`` or ``out`` directories below should run the steps needed to produce that file.

Some shortcuts are available:

- ``make depend`` installs dependencies
- ``make data`` produces cleaned data files.
- ``make chains`` produces all MCMC output
- ``make diagnostics`` extracts MCMC diagnostic statistics
- ``make figures`` produces all figures
- ``make tables`` produces all tables
- ``make macros`` produces macros used in LaTeX manuscript file
- ``make delcached`` removes cached compiled Stan MCMC programs, leading them to be recompiled from source
- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)

### Examining code

Examining the raw Stan code is the place to start to understand how Bayesian inference models have been specified. But note that parameters for the prior distributions are set at runtime rather than hard-coded into the ``.stan`` files, so that recompilation is not required when parameter choices are changed (this makes it easier to try the models using different priors, for sensitivity analysis).

Prior parameter choices are specified in the ``src/parameters`` directory.

Stan model code is modularized. Top level model files (found in ``src/stan``) are built by including code from module files found in subdirectories of ``src/stan``.

## Project structure when complete

Once the full analysis has been run, you should be able to find a full set of figures in ``out/figures`` and a set of results tables in ``out/tables``.

## Dependency installation guide
You will need a working R installation with the command line interpreter ``Rscript`` (macOS and Linux) or ``Rscript.exe`` (Windows). On mac and Linux, you can check that you have an accessible ``Rscript`` by typing ``which Rscript``at the command line and seeing if one is found.

If you do not have an R installation, you can install it from [the R project website](https://www.r-project.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

Once R is installed, you can automatically install all other dependencies (including the Hamiltonian Monte Carlo software Stan and its R interface rstan) on most systems using ``make``. In the top level project directory, type the following at the command line:

    make depend

Alternatively, you can run the script ``src/install_needed_packages.R`` manually. 

Note that installing Stan and RStan can be time-consuming, Stan is a large program that must be compiled from source. Some of the packages in the very valuable [tidyverse](https://www.tidyverse.org/) may also take some time to install.
