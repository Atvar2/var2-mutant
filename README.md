# Arabidopsis-var2-mutant

Collection of R codes or horse perl scripts to perfrom analysis and data visualization in research article "Exploring the probabilistic fates of chloroplasts in Arabidopsis by single-cell RNA sequencing".

|         |                                                                  |
| ------- | ---------------------------------------------------------------- |
| Authors | Junhui Chen ([chenjunhui](https://github.com/Atvar2))         |
| Email   | <chenjhbio@163.com>                                           |
| License | [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)               |


## Citations
If you used the codes or reference part of codes for your analysis, please kindly cited:
_Chai et al. (2022) Exploring the probabilistic fates of chloroplasts in Arabidopsis by single-cell RNA sequencing

## Dependece packages
Following are a list of R packages that are used by the analysis pipeline. Before implement the codes, please comfirmed the packages have been installed on your R platform.

- Seaurat version 4.1.0
- DoubletFinder version 2.02
- clusterProfiler version 4.0.0
- monocle version: 2.20.0
- ggplot2 version: 3.3.5

The other R packages such as dplyr, should depend on R environment variable space to install the responding packages.

## Installation

Just download the codes and implement in R platform under linux or window system. For packages installation, mainly intallated through [conda](https://docs.conda.io/en/latest/) or build-in funciont of R install.packages() 
