# soundDMR

There are a number of steps in our differential methylation pipeline that run through the use of in-house scripts.  


## Pre-requisites
We need to run Megalodon/Deep-Signal Plant/Bonito to calculate methylation levels at each cytosine for each individual. This package technically works with any ONT bed file with a format as mentioned in [here](https://www.encodeproject.org/data-standards/wgbs/), in the **Description of bedMethyl file section**. If focussing on a small part fo the genome, these bed files must be subset using bedtools intersect to focus on specific regions of the genome.

``` 
bedtools intersect -a [ONT_methyl.bed] -b [target_regions.bed] -wa > [ONT_methyl_subset.bed]
```
More information on bedtools intersect can be found [here](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).


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
Sample data is available on the [github](https://github.com/SoundAg/sounDMR) page which includes methyl bed files for 2 treated and control individuals each along with Experimental design file and Gene coordinates file for the ZoomFrame.


## DMR Analysis Functionalization

## Authors
Tom Cairns  
Jack Colicchio  
Keerthana N Prabhu


## Purpose
This project intends to create useful, robust functions for DMR analysis. It takes
in a megaframe, cleans with the help of the programmer, then runs the
random effects model for a given formula, and finally identifies the changepoints
in the Z score of interest from the random effects model. 


## Files
*DMR_analysis_workflow*

This file contains the actual workflow for the optimized DMR analysis pipeline. 
It can be used either as a reference for analysis or as the script for analysis
depending on the experiment.


## Package
The sounDMR2 package is contained within the sounDMR2 directory. The current
package can be downloaded from github in the .tar.gz format.


## Guide
This README will not contain a comprehensive guide to how to run the workflow,
but rather an overview. The full guide for our analysis can be found on google
drive.

#### Install and Load Package
The first step is to install the package from this repository. Download the
.tar.gz file and then load into into RStudio by going to Tools -> Install Packages.
Once the package is installed it can be loaded into the working analysis by
using `library(sounDMR2)`.

#### Create Megaframe
The first step is to generate a megaframe that contains methyl bed information
for each individual combined into a single data frame.

#### Create Zoomframe
DMR analysis requires some additional columns apart from just the methyl bed information.
This function creates,
Zeroth_pos : Position adjusted to have position 0 as ATG for the target;
Gene : This can inclkude either Gene Id or Gene names for tracking during DMR analysis;
Zoom_co : This keeps track of whether the given position is within the gene or is an adaptive region.

#### Clean and Rearrange Data
This step cleans the data read into the environment and creates "long" formats
of the count of methylated and unmethylated reads as well as the percent
methylation for each individual.

#### Creating Methylation Summary Data
This section is involved in manipulating the methylation data to summarize
information based on specific variables such as the Percent methylation of each
individual for each cytosine.

#### Creating Differential Methylation Output File
This section creates a data frame where all the output data will be stored. It
creates three columns for each individual: the Z score, the change in methylation
between the treatment and the control, and the read depth for an individual.

#### Group DMR Analysis
This section is where the model is run for our DMR analysis. Here we use a
binomial mixed effects model where the user has the option to change the fixed
and random effects passed into the model depending on the goals of the experiment.

#### Changepoint Analysis
The output from the previous step gets passed into our changepoint analysis
section. Here the user picks a column upon which to run the analysis. We
recommend using one of the Z_score columns created from the fixed effects provided
to the model.
