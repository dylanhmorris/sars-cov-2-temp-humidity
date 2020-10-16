# The effect of temperature and humidity on the stability of SARS-CoV-2 and other enveloped viruses
[Dylan H. Morris](https://dylanhmorris.com)(1\*), Kwe Claude Yinda (2\*), Amandine Gamble(3\*), Fernando W. Rossine(1), Qishen Huang(4), Trenton Bushmaker(2, 5), M. Jeremiah Matson(2, 6), Neeltje van Doremalen(2), Peter J. Vikesland(4), Linsey C. Marr(4), Vincent J. Munster(2), James O. Lloyd-Smith(3)

\* These authors contributed equally

1. Dept. of Ecology and Evolutionary Biology, Princeton University, Princeton, NJ, USA
2. Laboratory of Virology, Division of Intramural Research, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Hamilton, MT, USA
3. Dept. of Ecology and Evolutionary Biology, University of California, Los Angeles, Los Angeles, CA, USA
4. Department of Civil and Environmental Engineering, Virginia Tech, Blacksburg, VA, USA
5. Department of Microbiology and Immunology, Montana State University, Bozeman, MT, USA
6. Joan C. Edwards School of Medicine, Marshall University, Huntington, WV, USA


## Repository information
This repository accompanies the article "The effect of temperature and humidity on the stability of SARS-CoV-2 and other enveloped viruses" (DH Morris et al, 2020). It provides code for reproducing data analysis and theoretical modeling from the paper and recreating display figures.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work as below:

- DH Morris, KC Yinda, A Gamble et al. The effect of temperature and humidity on the stability of SARS-CoV-2 and other enveloped viruses. 2020.

Bibtex record:
```
@electronic{morris2020temperature,
    Author = {
      Dylan H. Morris AND 
      Kwe Claude Yinda AND 
      Amandine Gamble AND 
      Fernando W. Rossine AND
      Qishen Huang AND
      Trenton Bushmaker AND
      M. Jeremiah Matson AND 
      Neeltje van Doremalen AND 
      Peter J. Vikesland AND 
      Linsey C. Marr AND
      Vincent J. Munster AND
      James O. Lloyd-Smith},
      Title = {The effect of temperature and humidity on the stability of SARS-CoV-2 and other enveloped viruses},
      Date = {2020},
      URL = {https://github.com/dylanhmorris/sars-cov-2-temp-humidity}
}
```

## Article abstract 
Since emerging in late 2019, SARS-CoV-2 has caused a global pandemic, and it may become an endemic human pathogen. Understanding the impact of environmental conditions on SARS-CoV-2 viability and its transmission potential is crucial to anticipating epidemic dynamics and designing mitigation strategies. Ambient temperature and humidity are known to have strong effects on the environmental stability of viruses, but there is little data for SARS-CoV-2, and a general quantitative understanding of how temperature and humidity affect virus stability has remained elusive. Here, we characterise the stability of SARS-CoV-2 on an inert surface at a variety of temperature and humidity conditions, and introduce a mechanistic model that enables accurate prediction of virus stability in unobserved conditions. We find that SARS-CoV-2 survives better at low temperatures and extreme relative humidities; median estimated virus half-life was more than 24 hours at 10C and 40\% RH, but less than an hour and a half at 27C and 65\% RH. Our results highlight scenarios of particular transmission risk, and provide a mechanistic explanation for observed superspreading events in cool indoor environments such as food processing plants. Moreover, our model predicts observations from other human coronaviruses and other studies of SARS-CoV-2, suggesting the existence of shared mechanisms that determine environmental stability across a number of enveloped viruses.

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
