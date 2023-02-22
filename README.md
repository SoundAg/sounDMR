# soundDMR

There are a number of steps in our differential methylation pipeline that run through the use of in-house scripts.  


## Pre-requisites
We need to run Megalodon/Deep-Signal Plant/Bonito to calculate methylation levels at each cytosine for each individual. This package technically works with any-ONt bed file with a format as mentioned in [here](https://www.encodeproject.org/data-standards/wgbs/), in the **Description of bedMethyl file section**. If focussing on a small part fo the genome, these bed files must be subset using bedtools intersect to focus on specific regions of the genome.

## Installation
Make sure you have access to the latest version of the package and go through the installation steps below. This can be done either by downloading the tar file from our [github](https://github.com/SoundAg/DMR_analysis) page or by cloning it in your local environment.
This only needs to be done once.

<div class=".pkgdown-release">

``` r
# Install from CRAN
install.packages("sounDMR")
```

</div>

## Usage
`library(soundDMR)` will load the core packages:

- [tidyverse](https://tidyverse.org)
- [reshape2](https://www.rdocumentation.org/packages/reshape2/versions/1.4.4)
- [changepoint](https://cran.r-project.org/web/packages/changepoint/changepoint.pdf)
- [lme4](https://cran.r-project.org/web/packages/lme4/lme4.pdf)
- [Formula](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula)
- [glmmTMB](https://cran.r-project.org/web/packages/glmmTMB/glmmTMB.pdf)

``` r
library(sounDMR)

#>Loading required package: changepoint
#>Loading required package: zoo
#>
#>Attaching package: ‘zoo’
#>
#>The following objects are masked from ‘package:base’:
#>
#>    as.Date, as.Date.numeric
#>
#>Successfully loaded changepoint package version 2.2.4
#> See NEWS for details of changes.
#>Loading required package: Formula
#>Loading required package: lme4
#>Loading required package: Matrix
#>Loading required package: reshape2
#>Loading required package: tidyverse
#>── Attaching packages ────────────────────────────────────────── tidyverse 1.3.2 ──
#>✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
#>✔ tibble  3.1.8      ✔ dplyr   1.0.10
#>✔ tidyr   1.2.1      ✔ stringr 1.4.1 
#>✔ readr   2.1.3      ✔ forcats 0.5.2 
#>── Conflicts ───────────────────────────────────────────── tidyverse_conflicts() ──
#>✖ tidyr::expand() masks Matrix::expand()
#>✖ dplyr::filter() masks stats::filter()
#>✖ dplyr::lag()    masks stats::lag()
#>✖ tidyr::pack()   masks Matrix::pack()
#>✖ tidyr::unpack() masks Matrix::unpack()
#>Loading required package: glmmTMB

```

## Data availability
Sample data is available on the [github](https://github.com/SoundAg/sounDMR) page which includes methyl bed files for 2 treated and control individuals each along with Experimental-ID file and Gene-corodinates file for the ZoomFrame.


