<<<<<<< HEAD
[![Build
Status](https://travis-ci.org/dsrobertson/onlineFDR.svg?branch=master)](https://travis-ci.org/dsrobertson/onlineFDR)
[![codecov](https://codecov.io/gh/dsrobertson/onlineFDR/branch/master/graph/badge.svg)](https://codecov.io/gh/dsrobertson/onlineFDR)

onlineFDR <img src="man/figures/logo.png" align="right" />
==========================================================

`onlineFDR` allows users to control the false discovery rate (FDR) or
familywise error rate (FWER) for online hypothesis testing, where
hypotheses arrive sequentially in a stream. In this framework, a null
hypothesis is rejected based only on the previous decisions, as the
future p-values and the number of hypotheses to be tested are unknown.

Installation
------------

To install the latest (development) version of the onlineFDR package
from Bioconductor, please run the following code:

=======

[![Build
Status](https://travis-ci.org/dsrobertson/onlineFDR.svg?branch=master)](https://travis-ci.org/dsrobertson/onlineFDR)
[![codecov](https://codecov.io/gh/dsrobertson/onlineFDR/branch/master/graph/badge.svg)](https://codecov.io/gh/dsrobertson/onlineFDR)

# onlineFDR <img src="man/figures/logo.png" align="right" />

`onlineFDR` allows users to control the false discovery rate (FDR) or
familywise error rate (FWER) for online hypothesis testing, where
hypotheses arrive sequentially in a stream. In this framework, a null
hypothesis is rejected based only on the previous decisions, as the
future p-values and the number of hypotheses to be tested are unknown.

## Installation

To install the latest (development) version of the onlineFDR package
from Bioconductor, please run the following code:

>>>>>>> master
``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

<<<<<<< HEAD
# The following initializes usage of Bioc devel
BiocManager::install(version='3.12')
#> Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.2 (2020-06-22)
#> Old packages: 'annotate', 'AnnotationDbi', 'AnnotationFilter',
#>   'AnnotationForge', 'apeglm', 'backports', 'Biobase', 'BiocCheck',
#>   'BiocFileCache', 'BiocParallel', 'BiocStyle', 'biocViews', 'biomaRt',
#>   'Biostrings', 'biovizBase', 'bookdown', 'broom', 'BSgenome', 'bumphunter',
#>   'callr', 'car', 'Category', 'cli', 'clipr', 'clusterProfiler', 'coda',
#>   'codetools', 'cpp11', 'DALEX', 'data.table', 'dbplyr', 'DEFormats',
#>   'DelayedArray', 'deldir', 'derfinder', 'derfinderHelper', 'derfinderPlot',
#>   'DESeq2', 'DNAcopy', 'doParallel', 'DOSE', 'downlit', 'DT', 'e1071', 'edgeR',
#>   'enrichplot', 'ensembldb', 'fgsea', 'foreach', 'Formula', 'furrr', 'future',
#>   'gdsfmt', 'genefilter', 'geneLenDataBase', 'geneplotter', 'generics',
#>   'GenomeInfoDb', 'GenomeInfoDbData', 'GenomicAlignments', 'GenomicFeatures',
#>   'GenomicFiles', 'GenomicRanges', 'GEOquery', 'ggbio', 'ggfortify', 'globals',
#>   'GO.db', 'GOSemSim', 'goseq', 'GOstats', 'graph', 'graphlayouts', 'GSEABase',
#>   'gtsummary', 'GWASTools', 'haplo.stats', 'hardhat', 'HardyWeinberg',
#>   'htmlwidgets', 'ideal', 'igraph', 'IHW', 'insight', 'iterators',
#>   'kableExtra', 'KernSmooth', 'knitr', 'labeling', 'lava', 'lhs', 'limma',
#>   'lme4', 'lpsymphony', 'lubridate', 'MASS', 'meta', 'mgcv', 'mice',
#>   'modeldata', 'multcomp', 'nlme', 'officer', 'onlineFDR', 'openxlsx',
#>   'org.Hs.eg.db', 'OrganismDbi', 'parsnip', 'patchwork', 'pcaExplorer',
#>   'pkgmaker', 'profvis', 'ProtGenerics', 'ps', 'quantreg', 'quantsmooth',
#>   'qvalue', 'R6', 'raster', 'RBGL', 'RcppArmadillo', 'reactable', 'readr',
#>   'recipes', 'recount', 'RefManageR', 'regionReport', 'rentrez', 'renv',
#>   'Rgraphviz', 'Rhtslib', 'rlang', 'rmarkdown', 'rsample', 'Rsamtools',
#>   'RSQLite', 'rstan', 'rstudioapi', 'rtracklayer', 'sandwich', 'seriation',
#>   'sever', 'shinyFeedback', 'shinyWidgets', 'sjlabelled', 'SNPassoc', 'sp',
#>   'SQUAREM', 'statmod', 'SummarizedExperiment', 'survival', 'systemfonts',
#>   'testthat', 'textshaping', 'tibble', 'tidypredict', 'tidytext', 'tinytex',
#>   'topGO', 'V8', 'VariantAnnotation', 'vroom', 'withr', 'workflows', 'xfun',
=======
# The following initializes usage of Bioc
BiocManager::install()
#> Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.2 (2020-06-22)
#> Old packages: 'annotate', 'AnnotationDbi', 'AnnotationFilter',
#>   'AnnotationForge', 'apeglm', 'Biobase', 'BiocCheck', 'BiocFileCache',
#>   'BiocParallel', 'BiocStyle', 'biocViews', 'biomaRt', 'Biostrings',
#>   'biovizBase', 'bookdown', 'broom', 'BSgenome', 'bumphunter', 'car',
#>   'Category', 'clipr', 'clusterProfiler', 'coda', 'codetools', 'DALEX',
#>   'data.table', 'dbplyr', 'DEFormats', 'DelayedArray', 'deldir', 'derfinder',
#>   'derfinderHelper', 'derfinderPlot', 'DESeq2', 'DNAcopy', 'doParallel',
#>   'DOSE', 'downlit', 'DT', 'e1071', 'edgeR', 'enrichplot', 'ensembldb',
#>   'fgsea', 'foreach', 'Formula', 'furrr', 'future', 'gdsfmt', 'genefilter',
#>   'geneLenDataBase', 'geneplotter', 'GenomeInfoDb', 'GenomeInfoDbData',
#>   'GenomicAlignments', 'GenomicFeatures', 'GenomicFiles', 'GenomicRanges',
#>   'GEOquery', 'ggbio', 'ggfortify', 'ggraph', 'globals', 'GO.db', 'GOSemSim',
#>   'goseq', 'GOstats', 'graph', 'graphlayouts', 'GSEABase', 'gtsummary',
#>   'GWASTools', 'haplo.stats', 'hardhat', 'HardyWeinberg', 'here',
#>   'htmlwidgets', 'ideal', 'igraph', 'IHW', 'insight', 'iterators',
#>   'kableExtra', 'KernSmooth', 'knitr', 'lava', 'lhs', 'limma', 'lme4',
#>   'lpsymphony', 'lubridate', 'MASS', 'meta', 'mgcv', 'mice', 'modeldata',
#>   'multcomp', 'nlme', 'officer', 'onlineFDR', 'openxlsx', 'org.Hs.eg.db',
#>   'OrganismDbi', 'parsnip', 'patchwork', 'pcaExplorer', 'pkgmaker', 'profvis',
#>   'ProtGenerics', 'quantreg', 'quantsmooth', 'qvalue', 'raster', 'RBGL',
#>   'RcppArmadillo', 'reactable', 'readr', 'recipes', 'recount', 'RefManageR',
#>   'regionReport', 'rentrez', 'renv', 'Rgraphviz', 'Rhtslib', 'rmarkdown',
#>   'rprojroot', 'rsample', 'Rsamtools', 'RSQLite', 'rstan', 'rtracklayer',
#>   'sandwich', 'seriation', 'sever', 'shinyFeedback', 'shinyWidgets',
#>   'sjlabelled', 'SNPassoc', 'sp', 'SQUAREM', 'statmod', 'SummarizedExperiment',
#>   'survival', 'systemfonts', 'textshaping', 'tidypredict', 'tidytext',
#>   'tinytex', 'topGO', 'V8', 'VariantAnnotation', 'vroom', 'workflows', 'xfun',
>>>>>>> master
#>   'XVector', 'zlibbioc'

BiocManager::install("onlineFDR")
#> Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.2 (2020-06-22)
#> Installing package(s) 'onlineFDR'
#> 
#> The downloaded binary packages are in
<<<<<<< HEAD
#>  /var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T//Rtmpthf8vu/downloaded_packages
#> Old packages: 'annotate', 'AnnotationDbi', 'AnnotationFilter',
#>   'AnnotationForge', 'apeglm', 'backports', 'Biobase', 'BiocCheck',
#>   'BiocFileCache', 'BiocParallel', 'BiocStyle', 'biocViews', 'biomaRt',
#>   'Biostrings', 'biovizBase', 'bookdown', 'broom', 'BSgenome', 'bumphunter',
#>   'callr', 'car', 'Category', 'cli', 'clipr', 'clusterProfiler', 'coda',
#>   'codetools', 'cpp11', 'DALEX', 'data.table', 'dbplyr', 'DEFormats',
#>   'DelayedArray', 'deldir', 'derfinder', 'derfinderHelper', 'derfinderPlot',
#>   'DESeq2', 'DNAcopy', 'doParallel', 'DOSE', 'downlit', 'DT', 'e1071', 'edgeR',
#>   'enrichplot', 'ensembldb', 'fgsea', 'foreach', 'Formula', 'furrr', 'future',
#>   'gdsfmt', 'genefilter', 'geneLenDataBase', 'geneplotter', 'generics',
#>   'GenomeInfoDb', 'GenomeInfoDbData', 'GenomicAlignments', 'GenomicFeatures',
#>   'GenomicFiles', 'GenomicRanges', 'GEOquery', 'ggbio', 'ggfortify', 'globals',
#>   'GO.db', 'GOSemSim', 'goseq', 'GOstats', 'graph', 'graphlayouts', 'GSEABase',
#>   'gtsummary', 'GWASTools', 'haplo.stats', 'hardhat', 'HardyWeinberg',
#>   'htmlwidgets', 'ideal', 'igraph', 'IHW', 'insight', 'iterators',
#>   'kableExtra', 'KernSmooth', 'knitr', 'labeling', 'lava', 'lhs', 'limma',
#>   'lme4', 'lpsymphony', 'lubridate', 'MASS', 'meta', 'mgcv', 'mice',
#>   'modeldata', 'multcomp', 'nlme', 'officer', 'openxlsx', 'org.Hs.eg.db',
#>   'OrganismDbi', 'parsnip', 'patchwork', 'pcaExplorer', 'pkgmaker', 'profvis',
#>   'ProtGenerics', 'ps', 'quantreg', 'quantsmooth', 'qvalue', 'R6', 'raster',
#>   'RBGL', 'RcppArmadillo', 'reactable', 'readr', 'recipes', 'recount',
#>   'RefManageR', 'regionReport', 'rentrez', 'renv', 'Rgraphviz', 'Rhtslib',
#>   'rlang', 'rmarkdown', 'rsample', 'Rsamtools', 'RSQLite', 'rstan',
#>   'rstudioapi', 'rtracklayer', 'sandwich', 'seriation', 'sever',
#>   'shinyFeedback', 'shinyWidgets', 'sjlabelled', 'SNPassoc', 'sp', 'SQUAREM',
#>   'statmod', 'SummarizedExperiment', 'survival', 'systemfonts', 'testthat',
#>   'textshaping', 'tibble', 'tidypredict', 'tidytext', 'tinytex', 'topGO', 'V8',
#>   'VariantAnnotation', 'vroom', 'withr', 'workflows', 'xfun', 'XVector',
=======
#>  /var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T//RtmpRHIxqk/downloaded_packages
#> Old packages: 'annotate', 'AnnotationDbi', 'AnnotationFilter',
#>   'AnnotationForge', 'apeglm', 'Biobase', 'BiocCheck', 'BiocFileCache',
#>   'BiocParallel', 'BiocStyle', 'biocViews', 'biomaRt', 'Biostrings',
#>   'biovizBase', 'bookdown', 'broom', 'BSgenome', 'bumphunter', 'car',
#>   'Category', 'clipr', 'clusterProfiler', 'coda', 'codetools', 'DALEX',
#>   'data.table', 'dbplyr', 'DEFormats', 'DelayedArray', 'deldir', 'derfinder',
#>   'derfinderHelper', 'derfinderPlot', 'DESeq2', 'DNAcopy', 'doParallel',
#>   'DOSE', 'downlit', 'DT', 'e1071', 'edgeR', 'enrichplot', 'ensembldb',
#>   'fgsea', 'foreach', 'Formula', 'furrr', 'future', 'gdsfmt', 'genefilter',
#>   'geneLenDataBase', 'geneplotter', 'GenomeInfoDb', 'GenomeInfoDbData',
#>   'GenomicAlignments', 'GenomicFeatures', 'GenomicFiles', 'GenomicRanges',
#>   'GEOquery', 'ggbio', 'ggfortify', 'ggraph', 'globals', 'GO.db', 'GOSemSim',
#>   'goseq', 'GOstats', 'graph', 'graphlayouts', 'GSEABase', 'gtsummary',
#>   'GWASTools', 'haplo.stats', 'hardhat', 'HardyWeinberg', 'here',
#>   'htmlwidgets', 'ideal', 'igraph', 'IHW', 'insight', 'iterators',
#>   'kableExtra', 'KernSmooth', 'knitr', 'lava', 'lhs', 'limma', 'lme4',
#>   'lpsymphony', 'lubridate', 'MASS', 'meta', 'mgcv', 'mice', 'modeldata',
#>   'multcomp', 'nlme', 'officer', 'openxlsx', 'org.Hs.eg.db', 'OrganismDbi',
#>   'parsnip', 'patchwork', 'pcaExplorer', 'pkgmaker', 'profvis', 'ProtGenerics',
#>   'quantreg', 'quantsmooth', 'qvalue', 'raster', 'RBGL', 'RcppArmadillo',
#>   'reactable', 'readr', 'recipes', 'recount', 'RefManageR', 'regionReport',
#>   'rentrez', 'renv', 'Rgraphviz', 'Rhtslib', 'rmarkdown', 'rprojroot',
#>   'rsample', 'Rsamtools', 'RSQLite', 'rstan', 'rtracklayer', 'sandwich',
#>   'seriation', 'sever', 'shinyFeedback', 'shinyWidgets', 'sjlabelled',
#>   'SNPassoc', 'sp', 'SQUAREM', 'statmod', 'SummarizedExperiment', 'survival',
#>   'systemfonts', 'textshaping', 'tidypredict', 'tidytext', 'tinytex', 'topGO',
#>   'V8', 'VariantAnnotation', 'vroom', 'workflows', 'xfun', 'XVector',
>>>>>>> master
#>   'zlibbioc'
```

Alternatively, you can install the package directly from GitHub:

``` r
# install.packages("devtools") # If devtools not installed

devtools::install_github("dsrobertson/onlineFDR")
#> Downloading GitHub repo dsrobertson/onlineFDR@HEAD
<<<<<<< HEAD
#> rlang      (0.4.7  -> 0.4.8) [CRAN]
#> R6         (2.4.1  -> 2.5.0) [CRAN]
#> labeling   (0.3    -> 0.4.2) [CRAN]
#> backports  (1.1.10 -> 1.2.0) [CRAN]
#> tibble     (3.0.3  -> 3.0.4) [CRAN]
#> diffobj    (NA     -> 0.3.2) [CRAN]
#> rstudioapi (0.11   -> 0.13 ) [CRAN]
#> withr      (2.2.0  -> 2.3.0) [CRAN]
#> waldo      (NA     -> 0.2.3) [CRAN]
#> ps         (1.3.4  -> 1.4.0) [CRAN]
#> cli        (2.0.2  -> 2.1.0) [CRAN]
#> callr      (3.4.4  -> 3.5.1) [CRAN]
#> testthat   (2.3.2  -> 3.0.0) [CRAN]
#> generics   (0.0.2  -> 0.1.0) [CRAN]
#> cpp11      (0.2.1  -> 0.2.4) [CRAN]
#> Installing 15 packages: rlang, R6, labeling, backports, tibble, diffobj, rstudioapi, withr, waldo, ps, cli, callr, testthat, generics, cpp11
#> 
#> The downloaded binary packages are in
#>  /var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T//Rtmpthf8vu/downloaded_packages
#>      checking for file ‘/private/var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T/Rtmpthf8vu/remotes14169461feba4/dsrobertson-onlineFDR-8b8ecc8/DESCRIPTION’ ...  ✓  checking for file ‘/private/var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T/Rtmpthf8vu/remotes14169461feba4/dsrobertson-onlineFDR-8b8ecc8/DESCRIPTION’ (840ms)
=======
#> rprojroot (1.3-2 -> 2.0.2) [CRAN]
#> Installing 1 packages: rprojroot
#> 
#>   There is a binary version available but the source version is later:
#>           binary source needs_compilation
#> rprojroot  1.3-2  2.0.2             FALSE
#> installing the source package 'rprojroot'
#>      checking for file ‘/private/var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T/RtmpRHIxqk/remotes15329444a34b1/dsrobertson-onlineFDR-680fc9a/DESCRIPTION’ ...  ✓  checking for file ‘/private/var/folders/8q/bl_mhpkj7m13gjn159ct8rn40000gn/T/RtmpRHIxqk/remotes15329444a34b1/dsrobertson-onlineFDR-680fc9a/DESCRIPTION’
>>>>>>> master
#>   ─  preparing ‘onlineFDR’:
#>      checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  cleaning src
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  building ‘onlineFDR_1.7.1.tar.gz’
#>      
#> 
```

<<<<<<< HEAD
Documentation
-------------

Documentation is hosted at
<a href="https://dsrobertson.github.io/onlineFDR/" class="uri">https://dsrobertson.github.io/onlineFDR/</a>

To view the vignette for the version of this package installed in your
system, start R and enter:

=======
## Documentation

Documentation is hosted at <https://dsrobertson.github.io/onlineFDR/>

To view the vignette for the version of this package installed in your
system, start R and enter:

>>>>>>> master
``` r
browseVignettes("onlineFDR")
#> No vignettes found by browseVignettes("onlineFDR")
```

<<<<<<< HEAD
References
----------

Aharoni, E. and Rosset, S. (2014). Generalized alpha-investing:
definitions, optimality results and applications to public databases.
*Journal of the Royal Statistical Society (Series B)*, 76(4):771–794.

=======
## References

Aharoni, E. and Rosset, S. (2014). Generalized alpha-investing:
definitions, optimality results and applications to public databases.
*Journal of the Royal Statistical Society (Series B)*, 76(4):771–794.

>>>>>>> master
Foster, D. and Stine R. (2008). alpha-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 29(4):429-444.

Javanmard, A., and Montanari, A. (2015). On Online Control of False
<<<<<<< HEAD
Discovery Rate. *arXiv preprint*,
<a href="https://arxiv.org/abs/1502.06197" class="uri">https://arxiv.org/abs/1502.06197</a>.
=======
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.
>>>>>>> master

Javanmard, A., and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON:
an adaptive algorithm for online control of the false discovery rate.
*Proceedings of the 35th International Conference in Machine Learning*,
80:4286-4294.

Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false
discovery rate in biomedical research. *arXiv preprint*,
<<<<<<< HEAD
<a href="https://arxiv.org/abs/1809.07292" class="uri">https://arxiv.org/abs/1809.07292</a>.
=======
<https://arxiv.org/abs/1809.07292>.
>>>>>>> master

Robertson, D.S., Wildenhain, J., Javanmard, A. and Karp, N.A. (2019).
onlineFDR: an R package to control the false discovery rate for growing
data repositories. *Bioinformatics*, 35:4196-4199,
<<<<<<< HEAD
<a href="https://doi.org/10.1093/bioinformatics/btz191" class="uri">https://doi.org/10.1093/bioinformatics/btz191</a>.

Tian, J. and Ramdas, A. (2019a). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *arXiv preprint*,
<a href="https://arxiv.org/abs/1905.11465" class="uri">https://arxiv.org/abs/1905.11465</a>.

Tian, J. and Ramdas, A. (2019b). Online control of the familywise error
rate. *arXiv preprint*,
<a href="https://arxiv.org/abs/1910.04900" class="uri">https://arxiv.org/abs/1910.04900</a>.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous Online
Testing of Multiple Hypotheses. *arXiv preprint*,
<a href="https://arxiv.org/abs/1812.05068" class="uri">https://arxiv.org/abs/1812.05068</a>.
=======
<https://doi.org/10.1093/bioinformatics/btz191>.

Tian, J. and Ramdas, A. (2019a). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *arXiv preprint*,
<https://arxiv.org/abs/1905.11465>.

Tian, J. and Ramdas, A. (2019b). Online control of the familywise error
rate. *arXiv preprint*, <https://arxiv.org/abs/1910.04900>.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous Online
Testing of Multiple Hypotheses. *arXiv preprint*,
<https://arxiv.org/abs/1812.05068>.
>>>>>>> master
