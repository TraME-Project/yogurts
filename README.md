# yogurts

Parallelized implementations of the MSA and IPFP algorithms for the Yogurts Project.

## Installation

The quickest way to install yogurts is via the devtools package:
``` R
install.packages("devtools")
devtools:::install_github("TraME-Project/yogurts")
```

Note that yogurts requires compilation, so an appropriate development environment is necessary to install the package.
* Windows users should get [Rtools](https://cran.r-project.org/bin/windows/Rtools/). Please ensure that R and Rtools are installed to `C:\` (and not `C:\Program Files`), and that the PATH variables are set correctly (described during the Rtools installation process).
* Mac users should install Xcode and then check [here](https://cran.r-project.org/bin/macosx/tools/) for additional tools (Clang6 and gfortran).
