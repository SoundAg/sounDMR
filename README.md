# soundDMR

There are a number of steps in our differential methylation pipeline that run through the use of in-house scripts.  


## Pre-requisites
We need to run Megalodon/Deep-Signal Plant/Bonito to calculate methylation levels at each cytosine for each individual. This package technically works with any ONT bed file with a format as mentioned in [here](https://www.encodeproject.org/data-standards/wgbs/), in the **Description of bedMethyl file section**. If focusing on a small part fo the genome, these bed files must be subset using `bedtools intersect` to focus on specific regions of the genome.

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


## Workflow
This package follows the general workflow found in the *DMR_analysis_workflow.R* script. It follows the framework outline below:

#### Create Megaframe
The first step is to generate a megaframe that contains methyl bed information for each individual combined into a single data frame.

```
# Import the Geneco file
Geneco <- read.table(file.choose(), header=TRUE, sep=",")

# Find the methyl bed files within the working directory
All_methyl_beds <- list.files(path=".",pattern="*.bed")

# Generate the megaframe
Megaframe <- generate_megaframe(methyl_bed_list=All_methyl_beds, Sample_count = 0, 
                                Methyl_call_type="DSP",  File_prefix="Sample")
```

#### Create Zoomframe
DMR analysis requires some additional columns apart from just the methyl bed information.
This function creates,
Zeroth_pos : Position adjusted to have position 0 as ATG for the target;
Gene : This can inclkude either Gene Id or Gene names for tracking during DMR analysis;
Zoom_co : This keeps track of whether the given position is within the gene or is an adaptive region.

```
Zoomframe <- generate_zoomframe(gene_cord_df = Geneco, MFrame = Megaframe, 
                                Gene_col="Gene_name", filter_NAs=10, 
                                target_info=TRUE, gene_list = Geneco$Gene_name, 
                                File_prefix="Sample")
```

#### Clean and Rearrange Data
This step cleans the data read into the environment and creates "long" formats
of the count of methylated and unmethylated reads as well as the percent
methylation for each individual.

```
experimental_design_df <- read.table(file.choose(), header=TRUE, sep=",")
dmr_obj <- create_dmr_obj(Zoomframe, experimental_design_df)
```

#### Creating Methylation Summary Data
This section is involved in manipulating the methylation data to summarize
information based on specific variables such as the Percent methylation of each
individual for each cytosine. It also creates three columns for each individual:
the Z score, the change in methylation between the treatment and the control, and
the read depth for an individual. There is an optional `additional_summary_cols`
parameter that will allow us to include additional columns given a summary
statistic function such as `mean`, `sd`, or `var` and the name of the column on 
which to run the summary statistics. Multiple tuples can be added into this
list to create additional summary columns. 

```
methyl_summary <- create_methyl_summary(dmr_obj, control = 'C', treated = 'T',
                                        additional_summary_cols = list(c(sd, 'Group')))

# Option to subset methyl_summary
individuals_of_interest = unique(dmr_obj$experimental_design_df$Individual)
methyl_summary <- subset_methyl_summary(methyl_summary, 
                                        individuals_to_keep = individuals_of_interest)
```

#### DMR Analysis
This section is where the model is run for our DMR analysis. We recommend using a
binomial model when comparing groups (e.g. control vs treated) and a beta-binomial
model when comparing individuals vs the control group. The user has the option to change the fixed
and random effects passed into the model depending on the goals of the experiment.

```
# Group Model
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 3, 
                           control = 'C', model = 'binomial', 
                           analysis_type = 'group')
                           
# Individual Model
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 5, 
                           control = 'C', model = 'beta-binomial', 
                           analysis_type = 'individual')
```

#### Changepoint Analysis
The output from the previous step gets passed into our changepoint analysis
section. Here the user picks a column upon which to run the analysis. We
recommend using one of the Z_score columns created from the fixed effects provided
to the model.

```
# The genes of interest for the experiment
target_genes <- unique(dmr_obj$ZoomFrame_filtered$Gene)

methyl_summary <- changepoint_analysis(methyl_summary, CG_penalty = 9, 
                                       CHG_penalty = 4, CHH_penalty = 7, 
                                       target_genes = target_genes,
                                       save_plots = F,
                                       z_col = "Z_GroupT_small")
```


