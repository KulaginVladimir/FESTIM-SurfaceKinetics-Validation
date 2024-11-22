# FESTIM-SurfaceKinetics-Validation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14036908.svg)](https://doi.org/10.5281/zenodo.14036908)

## Overview

This repository contains scripts for V&V of the kinetic surface model implemented in FESTIM.

## Cases

The repository includes several validation cases and one verification test on the MMS basis.

[1] [H absorption in Ti](./H_Ti): The case reproduces simulation results of [Y. Shimohata et al.](https://www.sciencedirect.com/science/article/pii/S0920379621006098) on the H absorption in Ti. The simulations are based on experiments performed by [Y. Harooka et al.](https://www.sciencedirect.com/science/article/abs/pii/0022311581905663?via%3Dihub)

[2] [D adsorption on oxidised W](./D_WO): The case reproduces simulation results of [E. A. Hodille et al.](https://iopscience.iop.org/article/10.1088/1741-4326/ad2a29) on the D adsorption/desorption on/from oxidised W. The simulations are based on experiments performed by [A. Dunand et al.](https://iopscience.iop.org/article/10.1088/1741-4326/ac583a)

[3] [D atom exposure of self-damaged W](./D_damagedW): The case reproduces simulation results of [E. A. Hodille et al.](https://iopscience.iop.org/article/10.1088/1741-4326/aa5aa5/meta) on the isothermal D implantation in self-damaged W, followed by isothermal desorption. The simulations are based on experiments performed by [S. Markelj et al.](https://www.sciencedirect.com/science/article/pii/S0022311515303470?via%3Dihub)

[4] [D implantation in W-damaged EUROFER](./D_EUROFER): The case reproduces experimental results of [K. Schmid et al.](https://www.sciencedirect.com/science/article/pii/S2352179122002228) on the D implantation in damaged EUROFER, followed by TDS measurements. The FESTIM model is mainly based on the simulations of [K. Schmid et al.](https://www.sciencedirect.com/science/article/pii/S2352179123001333?via%3Dihub)

[5] [MMS](./MMS) test

## How to use

### Cloud

Jupyter notebooks can be run in the browser with [Binder](https://mybinder.org/v2/gh/KulaginVladimir/FESTIM-SurfaceKinetics-Validation/HEAD):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KulaginVladimir/FESTIM-SurfaceKinetics-Validation/HEAD)

### Locally

For a local use, clone this repository to your local machine.

```
git clone https://github.com/KulaginVladimir/FESTIM-SurfaceKinetics-Validation
```

Create and activate the correct conda environment with the required dependencies:

```
conda env create -f environment.yml
conda activate festim-surface-kinetics-vv-env
```

This will set up a Conda environment named `festim-surface-kinetics-vv-env` with all the required dependencies for running the FESTIM scripts.

Navigate to the desired case folder and run the Jupyter books using the activated Conda environment.

> [!NOTE]  
> LaTeX is required to reproduce paper-quality figures. To install required dependencies, run the following command in your terminal:
> ```
> sudo apt-get install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
> ```
