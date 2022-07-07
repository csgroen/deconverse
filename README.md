## `deconverse`: bulk RNA-seq deconvolution benchmark using single-cell reference profiles

**Note**: This package is still very immature, API's are not yet well documented and may change

To install, run:

```{r}
remotes::install_github("csgroen/ggheatmapper")
remotes::install_github("jamesotto852/ggdensity")
remotes::install_github("csgroen/deconverse")
```

`docker` or `singularity` must be available to run some deconvolution methods. To install, see: <https://docs.docker.com/get-docker/>

To install CIBERSORTx docker, run:

```{r}
library(deconverse)
install_cibersortx()
```

Usage instructions TBA
