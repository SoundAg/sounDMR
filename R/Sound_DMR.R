# Copyright 2023 Sound Agriculture Company
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create Gene Percent Data Frames
#'
#' A function to create different subsets of the LongPercent df depending on the
#' column of interest
#'
#' @param LongPercent (df) - containing the methylation data as percents in
#' long format using `dcast()`
#' @param x (str) the column to compare to Gene * Zeroth_pos
#' @param function_name (function) an aggregation function such as mean, sd, or var
#' @return dcast_output (df) - a subset of the LongPercent dataframe
#' @examples
#' # Generate sample data
#' LongPercent <- data.frame(
#'  Gene = c('Gene1', 'Gene1', 'Gene1', 'Gene2', 'Gene2'),
#'  Zeroth_pos = c(1, 1, 1, 2, 2),
#'  Group = c('C', 'C', 'T', 'C', 'T'),
#'  Percent = c(30, 50, 0, 75, 80))
#'
#' # Aggregate based on Group using mean
#' create_gene_percent_x(LongPercent, 'Group', mean)
#'
#' # Aggregate based on Group using standard deviation
#' create_gene_percent_x(LongPercent, 'Group', sd)
#'
#' @import reshape2
#' @import Formula
#' @export

create_gene_percent_x <- function(LongPercent, x = 'Chromosome',
                                  function_name = mean) {
  formula <- Formula::as.Formula(paste('Gene * Zeroth_pos ~', x))
  if (x %in% colnames(LongPercent)) {
    dcast_output <- tryCatch({
      reshape2::dcast(LongPercent, formula, function_name, value.var = 'Percent',
                      na.rm = T)
    }, error = function(e) {
      print('Did not find a recognized function. Try mean or sd')
      return(NA)
    })
  } else {
    print(paste(x, 'is not a column found in the data. Try another column name.'))
    return(NA)
  }
  
  return(dcast_output)
}


#' Create DMR Object
#'
#' A function to clean the input data of methylation data. It filters the data
#' frame to contain only the columns of interest passed in as an argument and
#' removes rows where the mean methylation is 0.
#'
#' @param ZoomFrame (df) containing the input methylation data
#' @param experimental_design_df (df) containing the experimental design
#' @param colnames_of_interest *optional* (list of strings) - the columns to keep
#' in the analysis
#' @return out (list) - is a list of dataframes. This is the dmr object.
#' @import tidyverse
#' @export

create_dmr_obj <- function(ZoomFrame = dataframe,
                           experimental_design_df = dataframe) {
  # Clean the data
  print('Step 1: removing rows that contain 0 methylation')
  Inputpers <- dplyr::select(ZoomFrame, starts_with("Per"))
  Inputpers$MeanMeth <- rowMeans(Inputpers, na.rm=TRUE)
  ZoomFrame_filtered <- ZoomFrame[Inputpers$MeanMeth != 0,]
  
  # Important columns
  colnames_of_interest <- c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                            'Zeroth_pos', 'Plant')
  
  #QC Make sure the correct columns are present
  print('Step 2: checking for missing columns in experimental id')
  for (col in colnames_of_interest) {
    if (!col %in% colnames(ZoomFrame_filtered) & !col %in% c('Individual', 'Plant')){
      print(paste(col, 'not found in the ZoomFrame_filtered'))
    }
  }
  
  # Create columns if they are not present
  if (!'Zeroth_pos' %in% colnames(ZoomFrame_filtered)){
    print('Zeroth_pos column not found, running analysis on Positon column')
    ZoomFrame_filtered$Zeroth_pos <- ZoomFrame_filtered$Position
  }
  if ('Chr' %in% colnames(ZoomFrame_filtered)) {
    print('Changing name of Chr column to Chromosome')
    ZoomFrame_filtered$Chromosome <- ZoomFrame_filtered$Chr
  }
  if (!'Gene' %in% colnames(ZoomFrame_filtered)) {
    print('Gene column not found, running analysis on Chromosome column')
    ZoomFrame_filtered$Gene <- ZoomFrame_filtered$Chromosome
  }
  if (!'Plant' %in% colnames(experimental_design_df) & !'Individual' %in% colnames(experimental_design_df)) {
    stop('No Plant or Individual column name found in Experimental ID. One is necessary to continue')
  }
  if (!'Plant' %in% colnames(experimental_design_df)) {
    print('Plant column not found in the Experimental ID, using Individual')
    experimental_design_df$Plant <- experimental_design_df$Individual
  }
  if (!'Individual' %in% colnames(experimental_design_df)) {
    print('Individual column not found in the Experimental ID, using Plant')
    experimental_design_df$Individual <- experimental_design_df$Plant
  }
  if (!'Individual_Name' %in% colnames(experimental_design_df)) {
    print('Individual_Name column not found in the Experimental_ID, using Plant')
    experimental_design_df$Individual_Name <- experimental_design_df$Plant
  }
  
  # Convert the experimental ID 'Plant' column to character
  experimental_design_df$Plant <- as.character(experimental_design_df$Plant)
  
  # We want to order the input frame in such a way that it will be easy to recreate analysis
  print('Step 3: reordering ZoomFrame_filtered')
  ZoomFrame_filtered <- ZoomFrame_filtered[order(ZoomFrame_filtered[,'Gene'], ZoomFrame_filtered[,'Zeroth_pos']), ]
  
  # Create the Longer dataframes
  #------------------------------
  print('Step 4: creating longer dataframes')
  #Make percent methylation data long, based on each individual and each position
  LongPercent <- pivot_and_subset(ZoomFrame_filtered, 'Per', 'Percent',
                                  colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                           'Zeroth_pos', 'Individual'))
  LongPercent$Individual <- gsub("PerMeth_", "", LongPercent$Individual, perl = T)
  
  #Very similar steps, but for #methylated and unmethylated reads
  LongUnMeth <- pivot_and_subset(ZoomFrame_filtered, 'UnMeth', 'RD',
                                 colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                          'Zeroth_pos', 'Individual'))
  LongMeth <- pivot_and_subset(ZoomFrame_filtered, 'Meth', 'RD',
                               colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                        'Zeroth_pos', 'Individual'))
  
  #Replacing read depth of NA with 0
  LongUnMeth$RD[is.na(LongUnMeth$RD)] = 0
  LongMeth$RD[is.na(LongMeth$RD)] = 0
  
  #Creating column of total read depth
  LongMeth$total_RD <- LongMeth$RD + LongUnMeth$RD
  LongMeth$Individual <- gsub("Meth_","", LongMeth$Individual, perl = T)
  
  #QC to check that LongPercent contains the correct number of rows
  if (nrow(LongPercent) == nrow(ZoomFrame_filtered) * nrow(experimental_design_df)) {
    print('LongPercent contains the expected number of rows')
  } else {
    print(paste('LongPercent contains', nrow(LongPercent), 'rows. Expected:',
                (nrow(ZoomFrame_filtered) * nrow(experimental_design_df)), 'rows'))
  }
  
  print('Step 5: merging long files with experimental_design_df to annotate')
  # Merge Long files with Exp_ID to annotate
  print('Step 6: annotating long files with experimental_design data frame')
  LongPercent <- dplyr::left_join(LongPercent, experimental_design_df, by = c('Individual' = 'ID'))
  LongMeth <- dplyr::left_join(LongMeth, experimental_design_df, by = c('Individual' = 'ID'))
  
  # Aggregate
  print('Step 7: aggregating by plant')
  LongPercent <- LongPercent %>%
    group_by(Gene, Zeroth_pos, Plant, Position, CX, Strand, Group, Chromosome) %>%
    summarize(Percent = mean(Percent, na.rm = T))
  
  LongMeth <- LongMeth %>%
    group_by(Gene, Zeroth_pos, Plant, Position, CX, Strand, Group, Chromosome) %>%
    summarize(total_RD = mean(total_RD, na.rm = T))
  LongMeth <- LongMeth[,c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                          'Zeroth_pos', 'Plant', 'total_RD', 'Group')]
  
  # Save all the outputs
  out <- list()
  out$ZoomFrame_filtered <- ZoomFrame_filtered
  out$LongPercent <- LongPercent
  out$LonUnMeth <- LongUnMeth
  out$LongMeth <- LongMeth
  out$experimental_design_df <- experimental_design_df
  out$Inputpers <- Inputpers
  
  return(out)
}


#' Subset Data Frames to Contain Only Columns of Interest
#'
#' Keep only the necessary columns from a given data frame
#'
#' @param df (df) data frame to subset
#' @inheritParams create_dmr_obj
#' @return subsetted_df (df) - data frame with only the columns of interest
#' @examples
#' # Generate data to subset
#' df <- data.frame(
#'  x = c(1, 2, 3, 4),
#'  y = c(6, 7, 8, 9),
#'  z = c(14, 23, 15, 4))
#'
#' # Subset the data
#' subset_cols(df, c('x', 'z'))
#'
#' @export

subset_cols <- function(df,
                        colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                 'Zeroth_pos', 'Individual')) {
  subsetted_df <- df[,c(colnames_of_interest, colnames(df)[ncol(df)])]
  
  return(subsetted_df)
}



#' Arrange Data into Long Format and Subset
#'
#' This function is used to get data into LongPercent and LongMeth formats to
#' allow for summary analysis. This function utilizes `tidyr::pivot_longer()`
#' and `subset_cols()` from this package.
#'
#' @param data (df) - the `ZoomFrame_filtered`
#' @param starts_with_cols (str) - the string pattern that the columns start with
#' @inheritParams tidyr::pivot_longer
#' @inheritParams create_dmr_obj
#' @return pivoted_and_subsetted_df (df) - the cleaned data frame that has been
#' pivoted to be longer and subsetted to only include the important columns.
#' @import tidyr
#' @export

pivot_and_subset <- function(data,
                             starts_with_cols = 'start',
                             values_to,
                             colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                      'Zeroth_pos', 'Individual')){
  # Pivot data
  pivoted_df <- tidyr::pivot_longer(data, cols = starts_with(starts_with_cols),
                                    names_to = 'Individual', values_to = values_to)
  
  # Subset columns
  pivoted_and_subsetted_df <- subset_cols(pivoted_df, colnames_of_interest)
  
  return(pivoted_and_subsetted_df)
}



#' Building Individual Methylation Columns
#'
#' Set up the output frame to include individual methylation change and Read
#' Depth columns
#'
#' @param Exp_ID_Treated (df) - containing the experimental design
#' @param Output_Frame (df) - the Output_Frame
#' @param GenePercentPlant (df) - summary dataframe built from `create_gene_percent_x()`
#' containing % methylation for each Plant
#' @param GenePercentGroup (df) - summary dataframe built from `create_gene_percent_x()`
#' containing % methylation for groups
#' @param GeneDepthPlant (df) - summary dataframe containing read depth for individuals
#' @param control (str) - the control variable in GenePercentGroup
#' @return Output_Frame (df) - the updated Output_Frame with 3 new columns for each individual in
#' the `Exp_IF_Treated` parameter
#' @export

create_cols_for_individuals <- function(Exp_ID_Treated,
                                        Output_Frame,
                                        GenePercentPlant,
                                        GeneDepthPlant,
                                        GenePercentGroup,
                                        control = 'C') {
  for(id in unique(Exp_ID_Treated$Plant)) {
    Output_Frame$NewZ <- 0
    # Name this column using trickery
    colnames(Output_Frame)[ncol(Output_Frame)] <- paste(id, "Z", sep="_")
    # Create new column that is difference in percent methylation between individual "i" and the control average
    Output_Frame[,ncol(Output_Frame) + 1] <- GenePercentPlant[,id] - GenePercentGroup[[control]]
    # Name this column
    names(Output_Frame)[ncol(Output_Frame)] <- paste(id, "MethChange", sep="_")
    # Create a column that is read depth for that individual
    Output_Frame[,ncol(Output_Frame) + 1] <- GeneDepthPlant[,id]
    # Name this column
    names(Output_Frame)[ncol(Output_Frame)] <- paste(id, "RD", sep="_")
  }
  return(Output_Frame)
}


#' Create Methyl Summary
#'
#' This function creates the methyl_summary that the group and individual DMR
#' functions rely on for their analysis. It includes some quality control checks
#' to make sure the output is in the correct format.
#'
#' @inheritParams create_cols_for_individuals
#' @param dmr_obj the dmr object containing the experimental design and raw data
#' @param treated the string representing the treated group
#' @param additional_summary_cols a nested list of parameters to create additional
#' summary columns. Each nested list will be a tuple where the first value is
#' the summary stats function (e.g. sd, mean, var) and the second value is the
#' string name of the column on which to run the summary statistic function
#'
#' @export
create_methyl_summary <- function(dmr_obj, control = 'C', treated = 'T',
                                  colnames_of_interest = c('Chromosome', 'Gene',
                                                           'Position', 'Strand', 'CX',
                                                           'Zeroth_pos', 'Plant'),
                                  additional_summary_cols = list()) {
  # Create the summary files
  GeneDepthPlant <- dcast(dmr_obj$LongMeth,Gene*Zeroth_pos~Plant,mean,
                          value.var = "total_RD", na.rm=TRUE)
  GenePercentGroup <- create_gene_percent_x(dmr_obj$LongPercent, 'Group', mean)
  GenePercentPlant <- create_gene_percent_x(dmr_obj$LongPercent, 'Plant', mean)
  
  # QC
  # Make sure the GenePercentX dfs are the same length as the Zoomframe_filtered
  # and in the same order
  if (sum(GenePercentPlant$Zeroth_pos == dmr_obj$ZoomFrame_filtered$Zeroth_pos) != nrow(dmr_obj$ZoomFrame_filtered) |
      sum(GenePercentPlant$Gene == dmr_obj$ZoomFrame_filtered$Gene) != nrow(dmr_obj$ZoomFrame_filtered)) {
    print('Output is in a different order. Try running again.')
  }
  
  # Create experimental_design_df_treated
  experimental_design_df_treated <- dmr_obj$experimental_design_df[dmr_obj$experimental_design_df$Group != control,]
  
  # Create the beginnings of the Output statistic data frame
  Output_Frame <- dmr_obj$ZoomFrame_filtered[,colnames_of_interest[-length(colnames_of_interest)]]
  
  # Create 3 new cols for each individual
  Output_Frame <- create_cols_for_individuals(experimental_design_df_treated,
                                              Output_Frame,
                                              GenePercentPlant, GeneDepthPlant,
                                              GenePercentGroup, control)
  
  # QC
  # Checking Output_Frame created the three new cols for each individual
  if (ncol(Output_Frame) == (length(colnames_of_interest) - 1) + 3 * length(unique(experimental_design_df_treated$Plant))) {
    print('Number of columns in methyl_summary is correct')
  } else {
    print('The number of columns in methyl_summary is incorrect. Double check there are no duplicates')
    print(paste('Expected:', (length(colnames_of_interest) - 1) + 3 * length(unique(experimental_design_df_treated$Plant)),
                'Found:', ncol(Output_Frame)))
  }
  
  #Add in any summary statistic columns, such as this one
  Output_Frame <- cbind(Output_Frame,GenePercentPlant[,3:ncol(GenePercentPlant)])
  Output_Frame$Treat_V_Control <- GenePercentGroup[[treated]] - GenePercentGroup[[control]]
  Output_Frame$Control <- GenePercentGroup[[control]]
  Output_Frame$Treated <- GenePercentGroup[[treated]]
  
  # Additional summary columns
  if (length(additional_summary_cols) != 0) {
    for (tuple in additional_summary_cols) {
      GenePercentX = create_gene_percent_x(dmr_obj$LongPercent, x = tuple[[2]],
                                           function_name = get(tuple[[1]]))
      GenePercentX <- GenePercentX %>% select(-c(Gene, Zeroth_pos))
      colnames(GenePercentX) = paste0(colnames(GenePercentX), '_',
                                      tuple[[1]])
      Output_Frame <- cbind(Output_Frame, GenePercentX)
    }
  }
  
  return(Output_Frame)
}


#' Subset Methyl Summary
#'
#' A function to subset the methyl_summary file to only include individuals of
#' interest in the analysis
#'
#' @param methyl_summary (df) - the input data
#' @param individuals_to_keep (list of strings) - the individuals to keep for
#' analysis
#' @return methyl_summary_subset (df) - the data only including the individuals
#' for analysis
#'
#' @export
subset_methyl_summary <- function(methyl_summary, individuals_to_keep) {
  methyl_summary_subset <- methyl_summary %>%
    select(Chromosome, Gene, Position, Strand, CX, Zeroth_pos,
           contains(individuals_to_keep), Treat_V_Control, Control)
  return(methyl_summary_subset)
}

#' Create Fixed Effects
#'
#' A function to create the string combining fixed effects that will be passed
#' into `create_function()`
#'
#' @param fixed (list of strings) - the fixed effects elements
#' @return fixed_effects (str) - the fixed effects properly formatted
#' to be passed into `create_function()`
#' @examples
#' # create multiple independent fixed effects
#' create_fixed_effects(c('Group', 'Individual'))
#' # [1] "Group + Individual"
#'
#' # Create fixed effects with interactions
#' create_fixed_effects(c('Group * Individual'))
#' # [1] "Group * Individual"
#' @export

create_fixed_effects <- function(fixed = c('effect1', 'effect2')) {
  fixed_effects = ''
  if (length(fixed) >= 1) {
    fixed_effects = fixed[1]
  }
  if (length(fixed) > 1) {
    for (i in fixed[2:length(fixed)]) {
      fixed_effects = paste0(fixed_effects, ' + ', i)
    }
  }
  return(fixed_effects)
}


#' Create Random Effects
#'
#' Create a string combining the list of random effects to be passed into the
#' `create_formula()` function
#'
#' @param random (list of strings) - the random effects elements
#' @return random_effects (str) - the random effects properly formatted
#' @examples
#' # Create independent random effects
#' create_random_effects(c("Group", "Individual"))
#' # [1] "(1 | Group) + (1 | Individual)"
#'
#' # Create random effects with an interaction
#' create_random_effects(c("Group * Individual"))
#' # [1] "(1 | Group * Individual)"
#'
#' @export

create_random_effects <- function(random = c('Group', 'ID')) {
  random_effects = ''
  if (length(random) >= 1) {
    random_effects = paste0('(1|', random[1], ')')
  }
  if (length(random) > 1) {
    for (i in random[2:length(random)]) {
      random_effects = paste0(random_effects, ' + ', '(1|', i, ')')
    }
  }
  return(random_effects)
}


#' Create Formula
#'
#' Combine the fixed and random effects into the formula to be used
#' in the model
#'
#' @param fixed (str) - the fixed effects variables. Note: this
#' value **cannot** be 'ID'. ID is used to merge data together downstream for
#' running the model so including the value in the effects will break the model.
#' @param random (str) - the random effects variables. Note: this
#' value **cannot** be 'ID'. ID is used to merge data together downstream for
#' running the model so including the value in the effects will break the model.
#' @return effects_formula (formula) - a formula of the mixed effects
#' @examples
#' # Create formula of independent and random effects
#' create_formula(fixed = c('Group'), random = c('Individual'))
#' # [1] "cbind(Meth, UnMeth) ~ Group + (1 | Individual)"
#'
#' create_formula(fixed = c('Group', 'Individual'), random = c('Plant * Gene'))
#' # [1] "cbind(Meth, UnMeth) ~ Group + Individual + (1 | Plant * Gene)"
#' @import stats
#' @export

create_formula <- function(fixed, random) {
  fixed_effects <- create_fixed_effects(fixed)
  random_effects <- create_random_effects(random)
  
  if (fixed_effects != '' & random_effects != '') {
    effects = paste(random_effects, '+', fixed_effects)
  } else if (fixed_effects != '' & random_effects == '') {
    effects = fixed_effects
  } else if  (fixed_effects == '' & random_effects != '') {
    effects = random_effects
  }
  #  Create the formula
  formula = stats::as.formula(paste0('cbind(Meth, UnMeth) ~ ', effects))
  return(formula)
}


#' Save Model Summary
#'
#' Update the `Output_Frame` dataframe for a given row with the summary statistics
#' provided by the model
#'
#' @param i (int) - the row number of the dataframe
#' @param Output_Frame (df) - all the information
#' @param model_summary (df) - the summary statistics from the model
#' @param ind_name (str) - *optional* containing the name of the z score column for
#' individual DMR analysis
#' @return Output_Frame (df) - the data frame with the proper row updated to include
#' the summary statistics
#' @export

save_model_summary <- function(i, Output_Frame, model_summary, ind_name = '') {
  # Get the row names
  row_names <- rownames(model_summary)[2:length(rownames(model_summary))]
  
  for (name in row_names) {
    if (ind_name == '') {
      # Group Analysis
      Output_Frame[i, paste0('Z_', name, '_small')] = model_summary[name, 'z value']
      Output_Frame[i, paste0('Estimate_', name, '_small')] = model_summary[name, 'Estimate']
      Output_Frame[i, paste0('StdErr_', name, '_small')] = model_summary[name, 'Std. Error']
    } else {
      # Individual Analysis
      Output_Frame[i, ind_name] = model_summary[name, 'z value']
    }
  }
  return(Output_Frame)
}


#' Run Binomial Model
#'
#' This function runs the glmer function using a binomial model for a given
#' optimizer function
#'
#' @param LM (df) - the information to put in the model, usually this is
#' in a "long" format
#' @param i (int) - the row number
#' @param Output_Frame (df) - all the summary information
#' @param formula (formula) - the formula to use in the model
#' @param optimizer_func (str) - the optimizer function to use in the model
#' @return Output_Frame (df) - the data frame containing updated summary information
#' @import lme4
#' @export
#'

run_binomial <- function(LM, i = int, formula,
                         optimizer_func = 'optimizer') {
  binom_model <- glmer(formula, data=LM, family = binomial,
                       glmerControl(check.conv.grad = .makeCC(action = "stop",
                                                              tol = 2e-3, relTol = NULL),
                                    optimizer=optimizer_func, optCtrl = list(maxfun = 500000)))
  # Outputs of this model
  model_summary <- as.data.frame(summary(binom_model)$coefficients)
  
  return(model_summary)
}


#' Run Model
#'
#' Function to run the correct model
#'
#' @param data (df) - the information to put in the model, usually this is
#' in "long" format
#' @inheritParams run_binomial
#' @param individual_name_z (str) - *optional* the name of the z column when running
#' individual DMR analysis
#' @return Output_Frame (df) - updated with the summary data from the model
#' @export
#'

run_model <- function(data, i, Output_Frame, formula, model_type,
                      individual_name_z = ''){
  if (model_type == 'binomial') {
    tryCatch({
      ith_model_summary <- run_binomial(data, i, formula, 'bobyqa')
      
      # If that model didn't converge, it tried again with a different optimizer, allows ~20% more model convergence
    },error = function(e){tryCatch({ print(paste(i, "No bobyqa Converge, trying Nelder"))
      # Run the model with Nelder_Mead optimizer
      ith_model_summary <- run_binomial(data, i, formula, 'Nelder_Mead')
      Output_Frame <- save_model_summary(i, Output_Frame, ith_model_summary,individual_name_z)
      
    }, error=function(e){print(paste(i, "No Converge"))})
    })
  } else if (model_type == 'beta-binomial') {
    tryCatch({
      #The model in the individual version will generally be this, but if complex design, potentially can be multi-factorial
      beta_binomial <- glmmTMB(formula, data=data,
                               family=betabinomial(link = "logit"),
                               control=glmmTMBControl(
                                 optimizer=optim, optArgs=list(method="BFGS")))
      
      # Save the model output
      ith_model_summary <- as.data.frame(summary(beta_binomial)$coefficients$cond)
      Output_Frame <- save_model_summary(i, Output_Frame, ith_model_summary,
                                         individual_name_z)
    }, error=function(e){
      #print(paste(i, "No Converge"))
      paste(i, 'No Converge')
    })
  } else {
    print('Please choose a model type of "binomial" or "beta-binomial".')
  }
  return(Output_Frame)
}

#' Group DMR Analysis
#'
#' To run a binomial model to compare the methylation between groups
#'
#' @param Output_Frame (df) - the read depth and methylation
#' change information
#' @param ZoomFrame_filtered (df) - the percent methylation information
#' @param experimental_design_df (df) - the experimental design
#' @inheritParams create_formula
#' @param reads_threshold (int) - the number of reads that are each
#' methylated and unmethylated. This is important since data containing only one
#' methylated read is unlikely to provide statistical power to our analysis. A
#' value of 3 here means that the read depth is *at least* the depth provided by
#' the threshold. Both the methylated and unmethylated samples need a read depth
#' greater than or equal to this threshold in order to be considered for the model.
#' @inheritParams subset_cols
#' @return Output_Frame (df) - Output_Frame with the summary statistics from the model
#' @export

group_DMR <- function(Output_Frame, ZoomFrame_filtered, experimental_design_df, fixed = c('Group'),
                      random = c('Plant'), reads_threshold = 3, model = 'binomial',
                      colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                               'Zeroth_pos', 'Individual')) {
  #  Create the model formula first
  formula <- create_formula(fixed, random)
  print(formula)
  
  #  Get the number of columns in the Output_Frame dataframe
  original_Output_Frame_col_number <- ncol(Output_Frame)
  
  # The modelling here is the most "delicate" part of the operation.  Options include:
  # (A) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment
  # (B) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment + Generation
  # (C) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment + Generation +Treatment*Generation
  # (D) cbind(Meth, UnMeth) ~ (1|Plant) + Gene_Expression
  # (E) cbind(Meth, UnMeth) ~ (1|Plant) + Phenotype
  
  # Loop to run groupwise analysis for each BP
  for(i in 1:nrow(Output_Frame)){
    
    ZoomFrame_filtered_temp <- ZoomFrame_filtered[i,]  
    # Make long version of input frame for the i'th cytosine row
    LM <- pivot_and_subset(ZoomFrame_filtered_temp, 'Meth', 'Meth',
                           colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                    'Zeroth_pos', 'Individual'))
    LM$Individual <- gsub("Meth_","", LM$Individual, perl = T)
    # For rows with at least X(3) methylated reads across all individuals, move forward
    if(sum(as.numeric(LM$Meth), na.rm=TRUE) >= reads_threshold){
      # Do the same thing for unmethylated reads
      LUM <- pivot_and_subset(ZoomFrame_filtered_temp, 'UnMeth', 'UnMeth',
                              colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                       'Zeroth_pos', 'Individual'))
      LM <- cbind(LM,LUM[,ncol(LUM)])
      # Merge with experimental_design_df to allow for DML testing for that base
      LM <- LM %>% left_join(experimental_design_df, by=c('Individual'='ID'))
      
      # Has to have Y unmethylated reads across all individuials
      if(sum(as.numeric(LM$UnMeth), na.rm=TRUE) >= reads_threshold){
        # Run the model and save the output
        Output_Frame <- run_model(LM, i, Output_Frame, formula, model)
        rm(LM, LUM, ZoomFrame_filtered_temp)
      }
    }
  }
  #  Replace NA values in these columns with 0s
  Output_Frame[,original_Output_Frame_col_number:ncol(Output_Frame)][is.na(Output_Frame[,original_Output_Frame_col_number:ncol(Output_Frame)])] = 0
  
  return(Output_Frame)
}



#' Individual DMR Analysis
#'
#' Purpose: to run individual DMR analysis
#'
#' @inheritParams group_DMR
#' @param control (str) - the control variable
#' @return Output_Frame (df) - a data frame with additional columns
#' @export

individual_DMR <- function(Output_Frame, ZoomFrame_filtered, experimental_design_df,
                           fixed = c('Group'), random = c('Individual'),
                           reads_threshold = 3, control = 'C', model = 'beta-binomial',
                           colnames_of_interest = c('Chromosome', 'Gene', 'Position',
                                                    'Strand', 'CX', 'Zeroth_pos', 'Individual')) {
  
  # Obtain the formula
  formula <- create_formula(fixed, random)
  print(formula)
  
  Exp_ID_Treated <- experimental_design_df[experimental_design_df$Group == 'T',]
  
  # Create a progress bar
  pb = txtProgressBar(min = 0, max = nrow(Output_Frame), initial = 0, style = 3)
  
  #This model compares each of the treated individuals with the whole group of control individuals.
  for(i in 1:nrow(Output_Frame)){
    #First steps up until "For k in ..." are the same as group models.
    LM <- pivot_and_subset(ZoomFrame_filtered[i,], 'Meth', 'Meth',
                           colnames_of_interest)
    LM$Individual <- gsub("Meth_", "", LM$Individual, perl = T)
    
    if(sum(as.numeric(LM$Meth), na.rm=TRUE) >= reads_threshold){
      LUM <- pivot_and_subset(ZoomFrame_filtered[i,], 'UnMeth', 'UnMeth',
                              colnames_of_interest)
      LM <- cbind(LM, LUM[,ncol(LUM)])
      LM <- merge(LM, experimental_design_df, by.x="Individual", by.y="ID", all=FALSE)
      
      #Here we return to the Exp_ID_Treated file, this file has only individuals that are in a treated group. Starting with the first row in this file.
      for(k in 1:nrow(Exp_ID_Treated)){
        #Make subset version of the longmeth file that only includes the control individual and you individual of interest.  Need to confirm Exp_ID_Treated$ID corresponds to LM$Individual
        LMEX <- LM[LM$Group == control | LM$Individual == Exp_ID_Treated$ID[k], ]
        individual_name <- Exp_ID_Treated$Individual_Name[k]
        individual_name_z <- paste0(individual_name, '_Z')
        
        #Now, if has more than 2 meth and unmeth cytosines, runs the statistical model
        if(sum(as.numeric(LMEX$UnMeth), na.rm=TRUE) > reads_threshold & sum(as.numeric(LMEX$Meth), na.rm=TRUE) > reads_threshold){
          # Run the model
          Output_Frame = run_model(LMEX, i, Output_Frame, formula, model,
                                   individual_name_z = individual_name_z)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(Output_Frame)
}


#' DMR Analysis
#'
#' This function runs the DMR analysis based on whether analysis is between
#' groups or for an individual compared to a group.
#'
#' @inheritParams individual_DMR
#' @param analysis_type (str) - either "individual" or "group"
#' @param dmr_obj (list) - the dmr_object containing the experimental design df and zoomFrame_filtered
#' @return Output_Frame (df)
#' @export

find_DMR <- function(Output_Frame, dmr_obj, fixed = c('Group'),
                     random = c('Plant'), reads_threshold = 3,
                     model, control = '', analysis_type) {
  # The required columns
  colnames_of_interest <- c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                            'Zeroth_pos', 'Plant')
  if (tolower(analysis_type) == 'group') {
    Output_Frame = group_DMR(Output_Frame, dmr_obj$ZoomFrame_filtered,
                             dmr_obj$experimental_design_df,
                             fixed = fixed,random = random, colnames_of_interest,
                             reads_threshold = reads_threshold, model = model)
  } else if (tolower(analysis_type) == 'individual') {
    Output_Frame = individual_DMR(Output_Frame, dmr_obj$ZoomFrame_filtered,
                                  dmr_obj$experimental_design_df,
                                  fixed = fixed, random = random,
                                  reads_threshold = reads_threshold,
                                  control = control, model = model)
  } else {
    print(paste(analysis_type, 'analysis type not supported. Please try "individual" or "group"'))
  }
  return(Output_Frame)
}


#'  Find the Columns for Changepoint Analysis
#'
#' This function is used to find the column names on which to run changepoint
#' analysis from the `group_DMR()` or `individual_DMR()` output frames.
#' These columns are usually the Z score column created from the `run_model`
#' model summary. The Z score column is based on the fixed effects that are
#' passed into the model. Note: any column can be passed into the changepoint
#' analysis function, the output of this function are the columns deemed to
#' provide the greatest use to the analysis.
#'
#' @param DMR_output (df) - the DMR Analysis output
#' @return z_cols (list of strings) - a list of strings representing the potential
#' column names that are suggested to be used in the changepoint analysis.
#' @export

find_changepoint_col_options <- function(DMR_output, Output_Frame = Output_Frame) {
  #added_cols <- colnames(DMR_output[c((ncol(Output_Frame) + 1):ncol(DMR_output))])
  z_cols <- colnames(DMR_output)[grepl('Z_', colnames(DMR_output))]
  
  print(z_cols)
}

#' Store Changepoint Information
#'
#' Add the important changepoint summary info to the dataframe from
#' which the changepoints were found
#'
#' @param whole_df (df) - all the information for one gene
#' @param changepoint (df) - the changepoint summary
#' @param col_name (str) - the name of the column upon which the changepoint analysis
#' is being conducted
#' @return whole_df (df) - information for one gene including the
#' changepoint summary
#' @export

add_changepoint_info <- function(whole_df, changepoint, col_name) {
  # The new column names
  MethRegionName <- paste0('MethRegion_', col_name)
  MethRegionLengthName <- paste0('MethRegionLength_', col_name)
  MethGroupName <- paste0('MethGroup_', col_name)
  
  if (typeof(changepoint) != 'NULL') {
    changepoint <<- changepoint
    changepoints <- changepoint@cpts
    
    #Add in the first and last cytosine if needed for the changepoint start and stop sites
    if ((1 %in% changepoints) & (nrow(whole_df) %in% changepoints)) {
      start <- changepoints
      stop <- c(changepoints[-1] - 1, nrow(whole_df))
    } else if ((1 %in% changepoints) & (!nrow(whole_df) %in% changepoints)) {
      start <- changepoints
      stop <- c(changepoints[-1] - 1, nrow(whole_df))
    } else if ((!1 %in% changepoints) & (nrow(whole_df) %in% changepoints)) {
      start <- c(1, changepoints[-length(changepoints)])
      stop <- c(changepoints[-length(changepoints)] - 1, nrow(whole_df))
    } else {
      start <- changepoints
      stop <- changepoints - 1
    }
    
    BP_CG <- data.frame(cbind(start,stop))
    BP_CG$DifMean <- param.est(changepoint)$mean
    
    # Add the start and stop sites at a given position
    for(i in 1:nrow(BP_CG)){
      BP_CG[i, 'BPStart'] <- whole_df[BP_CG[i,'start'], 'Position']
      BP_CG[i, 'BPStop'] <- whole_df[BP_CG[i,'stop'], 'Position']
    }
    
    # Find the MethRegion, MethGroupName, and MethRegionLength
    whole_df[MethRegionName] = 0
    whole_df[MethRegionLengthName] = 0
    for(i in 1:nrow(BP_CG)){
      whole_df[c(BP_CG[i, 'start']:BP_CG[i, 'stop']),MethRegionName] <- BP_CG[i, 'DifMean']
      whole_df[c(BP_CG[i, 'start']:BP_CG[i, 'stop']),MethRegionLengthName] <- abs(as.numeric(BP_CG[i, 'BPStop']) -
                                                                                    as.numeric(BP_CG[i, 'BPStart']))
      whole_df[c(BP_CG[i, 'start']:BP_CG[i, 'stop']),MethGroupName] <- i
    }
    
  } else {
    if (nrow(whole_df) == 1) {
      #This is in case there are no changepoints bc the data isn't big enough
      whole_df[MethRegionName] = whole_df[[col_name]]
      whole_df[MethRegionLengthName] = 1
      whole_df[MethGroupName] = 1
    } else {
      whole_df[nrow(whole_df) + 1,] <- NA
      whole_df[1, MethRegionName] = whole_df[[col_name]]
      whole_df[1, MethRegionLengthName] = 0
      whole_df[1, MethGroupName] = 1
    }
  }
  
  return(whole_df)
}

#' Find CPT Mean
#'
#' Purpose: find the changepoints based on mean values and include error catches
#'
#' @inheritParams changepoint_analysis
#' @param data (df) - the dataframe that has been subsetted for the gene and cytosine
#' context
#' @param penalty (int) - the penalty value to include for the changepoint analysis
#' @return changepoint_object (S4) - containing the changepoint information
#' @import changepoint
#' @export

find_cpt_mean <- function(data, z_col, penalty){
  changepoint_object <- tryCatch(
    {
      # Filter out any rows containing NAs
      data[is.na(data[[z_col]]) == 0, ]
      
      # Run the changepoing analysis
      changepoint::cpt.mean(data[[z_col]], method = 'PELT',
                            penalty = 'Manual', pen.value = penalty)
    }, error = function(e) {
      #print(e)
      print('No changepoints found, data contains 1 or 0 rows')
      return(NULL)
    })
  return(changepoint_object)
}


#' Plot Changepoints
#'
#' Plot the changepoint data
#'
#' @param data (df) - containing the cytosines with their percent methylation
#' @param changepoint_obj (S4) - the S4 object created from finding the changepoints
#' @param gene_name (str) - the name of the gene
#' @param penalty_val (int) - the penalty value being used
#' @param cyt_context (str) - the cytosine context
#' @param z_col (str) - the column upon which the changepoints are calculated
#' @return changepoint_plot (plot) - the plot of the changepoints
#' @import ggplot2
#' @import tidyverse
#' @export

plot_changepoints <- function(data, changepoint_obj, gene_name, penalty_val,
                              cyt_context, z_col) {
  tryCatch({
    # Get the segment data
    segment_data = data.frame(matrix(nrow = 0, ncol = 4))
    colnames(segment_data) = c('x', 'xend', 'y', 'yend')
    for (i in 1:length(changepoint_obj@cpts)) {
      group_data = data[data[[paste0('MethGroup_', z_col)]] == i, 'Zeroth_pos']
      segment_data[i, 'x'] = min(group_data)
      segment_data[i, 'xend'] = max(group_data)
      segment_data[i, 'y'] = changepoint_obj@param.est$mean[i]
      segment_data[i, 'yend'] = changepoint_obj@param.est$mean[i]
    }
    
    # Create basic plot
    plot <- data %>%
      ggplot(aes(x = Zeroth_pos, y = !!sym(z_col))) +
      geom_line(linewidth = 0.5) +
      theme(panel.background = element_blank(), axis.line = element_line()) +
      labs(title = paste('Changepoints for', cyt_context, gene_name, '; penalty',
                         penalty_val)) +
      scale_x_continuous(expand = c(0, 0))
    
    # Add in the segments
    plot <- plot + geom_segment(data = segment_data,
                                aes(x = x, y = y, xend = xend, yend = yend),
                                color = 'red')
    
    return(plot)
  }, error = function(e) {
    return(paste('Figure could not be created for:', gene_name))
  })
}


#' Changepoint Analysis
#'
#' Purpose: run the changepoint analysis for each gene and each cytosine
#' context
#'
#' @param whole_df (df) - the information such as genes, z scores, etc
#' @param CG_penalty,CHG_penalty,CHH_penalty (int) - penalty values for the
#' change point analysis. The higher the value the fewer changepoints that will
#' be created. The lower the value the more changepoints that will be created.
#' @param target_genes (list of strings) - a list of the target genes
#' @param save_plots (boolean) - if TRUE save plots to working directory, default
#' is FALSE
#' @param z_col (str) - the column to run the changepoint analysis on. This
#' can be any column, but for DMR analysis we recommend using the z scores for
#' a fixed effect variable.
#' @return everything (df) - data frame containing the mean changepoint value for the
#' `z_col` column
#' @import tibble
#' @export

changepoint_analysis <- function(whole_df,
                                 CG_penalty = int,
                                 CHG_penalty = int,
                                 CHH_penalty = int,
                                 target_genes = c(),
                                 save_plots = FALSE,
                                 z_col = 'column') {
  everything <- tibble::as_tibble(matrix(ncol = 8))
  colnames(everything) = c('index', 'CX', 'Zeroth_pos', 'Gene', paste0('MeanMeth_',z_col),
                           paste0('MeanMethRegion_',z_col), paste0('MethRegionLength_',z_col),
                           paste('MethylGroup',z_col))
  #  remove the first row
  everything <- everything[-1,]
  genes <- unique(whole_df$Gene)
  large_genes <- c()
  
  # Create a progress bar
  print('Running changepoint analysis:')
  pb = txtProgressBar(min = 0, max = length(genes), initial = 0, style = 3)
  #  Work through the genes and cytosine contexts
  for (i in 1:length(genes)) {
    #print(gene)
    gene <- genes[i]
    # Create the filtered df to contain the gene of interest
    gene_df <- whole_df[whole_df$Gene == gene,]
    
    # Create dfs for cytosine contexts
    gene_df.CG <- gene_df[gene_df$CX == 'CG',]
    gene_df.CHG <- gene_df[gene_df$CX == 'CHG',]
    gene_df.CHH <- gene_df[gene_df$CX == 'CHH',]
    
    # Find the changepoints
    x.CG <- find_cpt_mean(gene_df.CG, z_col, CG_penalty)
    x.CHG <- find_cpt_mean(gene_df.CHG, z_col, CHG_penalty)
    x.CHH <- find_cpt_mean(gene_df.CHH, z_col, CHH_penalty)
    
    # Add in the changepoint info
    gene_df.CG <- add_changepoint_info(gene_df.CG, x.CG, z_col)
    gene_df.CHG <- add_changepoint_info(gene_df.CHG, x.CHG, z_col)
    gene_df.CHH <- add_changepoint_info(gene_df.CHH, x.CHH, z_col)
    
    # Create the plots
    if (gene %in% target_genes & nrow(gene_df) < 100000) {
      plot.CG <- plot_changepoints(gene_df.CG, x.CG, gene, CG_penalty, 'CG', z_col)
      plot.CHG <- plot_changepoints(gene_df.CHG, x.CHG, gene, CHG_penalty, 'CHG', z_col)
      plot.CHH <- plot_changepoints(gene_df.CHH, x.CHH, gene, CHH_penalty, 'CHH', z_col)
      print(plot.CG)
      print(plot.CHG)
      print(plot.CHH)
      if (save_plots == T) {
        png(file=paste0(gene, '_CG_changepoints.png'))
        print(plot.CG)
        dev.off()
        png(file=paste0(gene, '_CHG_changepoints.png'))
        print(plot.CHG)
        dev.off()
        png(file=paste0(gene, '_CHH_changepoints.png'))
        print(plot.CHH)
        dev.off()
      }
    } else if (gene %in% target_genes & nrow(gene_df) >= 100000) {
      large_genes <- c(large_genes, gene)
    }
    
    # Combine everything into the everything dataframe (rename)
    everything <- rbind(everything, gene_df.CG, gene_df.CHG, gene_df.CHH)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  for (gene in large_genes) {
    print0(gene, ' has more than 100,000 cytosines. Plot not created')
  }
  
  return(everything)
}


#' Sound Score
#' @description
#' Takes in an Output_Frame that has been broken in to changepoint regions based on a specific test statistic of interest, and creates an aggregated changepoint region file that includes summary statistics for each region, and the "sound score" a measure of the strength of a DMR within that region.
#'
#' @param changepoint_OF (df) - the OF file that has had changepoint calling ran on it
#' @param Statistic (str) - name of the test statistic that changepoint was run on
#' @param Per_Change (str) - Column that captures the change in methylation between two groups of interest
#' @param other_columns (str) - Names of any other columns that one is interested in calculating region averaqges on.
#' @param CF (str) - True/False if there is a custom function to be calculated in addition to soundscore
#' @param User_Function (str) - The function to use for the custome function
#' @return Aggregated_Changepoint_Object (list) - File that contains information about methylation patterns within each changepoint region
#' @export



sound_score <- function(changepoint_OF = dataframe, Statistic="Z_GroupT_small",
                        Per_Change = "Treat_V_Control", Control="Control",
                        other_columns=c("Estimate_GroupT_small"), CF=FALSE, UserFunction=NA) {
  # Determine proper column names givwen the test statistic of interest
  MethGroup <- paste("MethGroup_",Statistic, sep = "")
  MethRegion_Z <- paste("MethRegion_",Statistic, sep = "")
  MethRegion_Length <- paste("MethRegionLength_",Statistic, sep = "")
  # Create vector of columns that include information that the user wants aggregated across each changepoint region.  Includes Statistic, Per_Change, region length, and any other columns of interest.
  keep_cols <- c(MethRegion_Z, MethRegion_Length,Per_Change, Control,other_columns )
  # Create new dataframe that includes a new column that has a column for every unique changepoint region
  cp_OF <- within(changepoint_OF, cp_group <- paste(Gene,CX,changepoint_OF[[MethGroup]], sep='_'))
  #Calculate statistics and aggregate for every region
  RegionStats <- cp_OF %>%
    group_by(cp_group) %>%
    mutate(Count = n()) %>%
    group_by(cp_group, Gene,CX, Count) %>%
    summarise_at(vars(one_of(keep_cols)), mean)
  Ag_Pos <- cp_OF %>%
    group_by(cp_group) %>%
    summarise(Start=min(Zeroth_pos), Stop=max(Zeroth_pos))
  RegionStats <- cbind(Ag_Pos, RegionStats[,-1])
  #Calculate Sound Statistic
  RegionStats$dmr_score<-(((RegionStats$Count)^(1/3))*(abs(RegionStats[[MethRegion_Z]])*abs(RegionStats[[Per_Change]]))^(1/2))
  RegionStats$dmr_score2<-(((RegionStats$Count)^(1/3))*(abs(RegionStats[[MethRegion_Z]])*abs(asin(sqrt(RegionStats[[Per_Change]]/100+RegionStats[[Control]]/100))-asin(sqrt(RegionStats[[Control]]/100)))^(1/2)))
  # if(CF==TRUE){
  #   RegionStats$CustomFunction<<-UserFunction
  #   RegionStats$CustomFunction_Percentile<-ecdf(RegionStats$CustomFunction)(RegionStats$CustomFunction)
  # }
  for(i in 1:nrow(RegionStats)){
    if (!is.na(RegionStats[[Per_Change]][i])) {
      if(RegionStats[[Per_Change]][i]<0){
        RegionStats$dmr_score[i]<-RegionStats$dmr_score[i]*-1
        RegionStats$dmr_score_Percentile<-ecdf(RegionStats$dmr_score)(RegionStats$dmr_score)
        RegionStats$dmr_score2[i]<-RegionStats$dmr_score2[i]*-1
        RegionStats$dmr_score2_Percentile<-ecdf(RegionStats$dmr_score2)(RegionStats$dmr_score2)
      }
    }
  }
  plot(RegionStats$dmr_score2_Percentile, RegionStats$dmr_score2)
  plot <- RegionStats %>%
    ggplot(aes(x = .data[[MethRegion_Z]], y = Count)) +
    geom_point(aes(color = dmr_score2))
  print(plot)
  SS_Obj <- list(RegionStats, cp_OF)
  names(SS_Obj) <- c("region_summary", "methyl_summary")
  return(SS_Obj)
}

#' split_by_chromosome
#' @description
#' A function to split a bed file into multiple file by chromosome.
#' It will create a new directory per each chromosome.
#' @param input_file (str) - A string with the name of the input file.
#' @return output_filelist(list) - A list with the bedfile names of the output files.
#' @export

split_by_chromosome <- function(input_file) {

  output_filelist <- list()
  # get only dir
  fields <- strsplit(input_file, "/")[[1]]
  input_dir <- paste(fields[1:(length(fields) - 1)], collapse = "/")
  # Get name without extension
  base_name <- tools::file_path_sans_ext(basename(input_file))
  # Open the input file for reading
  con <- file(input_file, "r")
  # Create a list to store the file connections for each chromosome output file
  output_files <- list()
  # Read the input file line by line and process each line
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break # Exit loop if end of file is reached
    # Split the line by tab to get the chromosome
    fields <- strsplit(line, "\t")[[1]]
    chromosome <- fields[1]
    # Create the output directory for the chromosome if not already created
    if (!dir.exists(file.path(input_dir, paste0("chr_", chromosome)))) {
      dir.create(file.path(input_dir, paste0("chr_", chromosome)))
    }
    # Create the output file for the chromosome if not already opened
    if (!(chromosome %in% names(output_files))) {
      output_file <- file.path(input_dir, paste0("chr_", chromosome), paste0(base_name, ".bed"))
      output_filelist <- append(output_filelist, output_file)
      output_files[[chromosome]] <- file(output_file, "w")
    }
    # Write the line to the corresponding output file
    writeLines(line, output_files[[chromosome]])
  }
  # Close all output file connections
  for (chr in names(output_files)) {
    close(output_files[[chr]])
    cat(paste("Chromosome", chr, "data has been written to", output_file, "\n"))
  }
  # Close the input file connection
  close(con)
  
  return(output_filelist)
}

#' split_by_chunk
#' @description
#' A function to split a bed file into multiple chunks.
#'
#' @param input_file (str) - A string witht the name of the input file.
#' @param chunk_size (int) - Maximum amount of nucleotides per chunk.

split_by_chunk <- function(input_file, chunk_size, output_dir = "./chunks/") {
  # Get base name
  base_name <- tools::file_path_sans_ext(basename(input_file))
  # Create the output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # Read the BED file into a data frame
  bed_df <- read.table(input_file, header = FALSE, col.names = c("chromosome", "start", "end", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  nucls <- tail(bed_df$start, n = 1)
  # Get number of chunks
  total_chunks <- ceiling(nucls/chunk_size)
  # Initialize variables for chunk creation
  lower_val <- 0
  upper_val <- chunk_size
  for (chunk_n in 1:total_chunks) {
    chunk <- subset(bed_df, start > lower_val & start <= upper_val)
    output_file <- file.path(output_dir, paste0(base_name, "_chunk_", chunk_n, ".bed"))
    write.table(chunk, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)    
    lower_val <- lower_val + chunk_size
    upper_val <- upper_val + chunk_size
  }
}


#' get_standard_methyl_bed
#' @description
#' A function to create a data frame for every individual in the experimental design without having to re run for every individual separately.
#' This function takes in the methyl_bed file, subsets and then creates Methylated and unmethylated counts for each position to be used in the next steps.
#' @param Methyl_bed (df) - Data frame containing the ONT methylation calls in a bed file for each individual separately
#' @param Sample_ID (str) - A string that is used inplace of sample name to keep it uniform. We are using an enumerator to generate this based on the number of samples/individuals in the experiment.
#' @param Methyl_call_type (str) - A string that includes information about the type of run. Currently this package works on Megalodon , DSP (DeepSignal Plant), Dorado and Bonito. Default call type is Dorado.
#' @param max_read_depth (int) - This parameter can be used to filter out rows with high red depth such as highly repetative regions to not be included in the analysis. Default is 0. 
#' @return Methyl_bed_sub (df) - Standard data frame methyl file for every individual with Meth, Unmeth and Per_Meth columns
#' @import tidyverse
#' @import stringr
#' @export

get_standard_methyl_bed <-function(Methyl_bed="Methyl.bed", Sample_ID = "S1", Methyl_call_type="Dorado", max_read_depth=100) {
  
  # Extract columns of interest based on which process was run
  #Columns of Interest include Chromosome|Position|Strand|Total_reads|Percent_Methylation|Cytosine_context
  
  if (Methyl_call_type %in% c('DSP', 'Bonito', 'Dorado')){
    
    Methyl_bed <- Methyl_bed[,c(1,2,3,4,5,6)]
  }
  else if(Methyl_call_type=="Megalodon"){
    Methyl_bed <- Methyl_bed[,c(1,2,3,4,5,7)]
  }
  #Assign column names
  colnames(Methyl_bed) <- c("Chromosome","Position",paste("Strand", Sample_ID, sep="_"),"Tot_reads", paste("PerMeth", Sample_ID, sep="_") ,paste("CX", Sample_ID, sep="_"))
  #calculate Meth and Unmeth reads from total reads. This is necessary for downstream Differential Methylation Region (DMR) analysis.
  Methyl_bed[[paste("Meth", Sample_ID, sep="_")]] <- round( (Methyl_bed[[paste("PerMeth", Sample_ID, sep="_")]] * Methyl_bed$Tot_reads)/100 )
  Methyl_bed[[paste("UnMeth", Sample_ID, sep="_")]] <- (Methyl_bed$Tot_reads - Methyl_bed[[paste("Meth", Sample_ID, sep="_")]])
  Methyl_bed <- Methyl_bed %>% dplyr::filter(Tot_reads< max_read_depth)
  Methyl_bed <- Methyl_bed %>% select(-Tot_reads)
  
  return(Methyl_bed)
}


#' generate_megaframe
#' @description
#' A function to create a single-combined data frame from individual methyl beds in the experiment
#' This function only works with bedfiles output from only either of the three methylation call algorithms - DeepSignal Plant(DSP), Megalodon and Bonito.
#'
#' @param methyl_bed_list (list) - ONT methyl bed filenames for each individual contained within the directory. This will just be a list of bedfile names.
#' Hint : The input will be the "methyl_bed_list" vector that you create in the previous step.
#' @param Sample_count (int) - This is required to assign proper alphabet codes. If you need to include the samples from a previous round, then enter the total number of samples from the previous round here. Default is 0. By default alphabetizing starts with 'A'.
#' @param Methyl_call_type (str) - A string that includes information about the type of run. Currently this package works on Megalodon , DSP (DeepSignal Plant) and Bonito.
#' @param File_prefix (Flexible str) - This is to add a prefix to all the files that get exported and saved to the working directory while running the function.
#' @return Megaframe(df) - Clean data frame containing combined methyl bed information for every individual in the experiment.
#' @import tidyverse
#' @import stringr
#' @import data.table
#' @examples
#' # Basic usage for methyl_call_type
#' # generate_megaframe(Methyl_call_type="DSP") # OR
#' # get_standard_methyl_bed(Methyl_call_type="Megalodon") # OR
#' # get_standard_methyl_bed(Methyl_call_type="Bonito")
#' @export



generate_megaframe <- function(methyl_bed_list=All_methyl_beds, Sample_count = 0, Methyl_call_type="Dorado", File_prefix="", max_read_depth=100){
  
  #QC
  QC <- missing(methyl_bed_list)
  if(QC==TRUE){
    stop("methyl_bed_list parameter cannot be empty.
         Hint: Use methyl_bed_list vector you create in the previous step")
  }
  
  if(Methyl_call_type==""){
    cat("Methylation call type not given, using default type i.e Dorado \n")
  }
  
  
  if (!(Methyl_call_type %in% c('DSP', 'Megalodon', 'Bonito', 'Dorado'))){
    stop("Methylation call not recognized, use 'DSP' or 'Megalodon' or 'Bonito' or 'Dorado', exiting!")
    
  }
  
  sample_number_list <- c()
  for (i in 1:(Sample_count + length(methyl_bed_list)) ){
    S_enumerator <- paste("S",i,sep="")
    sample_number_list[(length(sample_number_list) + 1)] <- list(S_enumerator)
  }
  
  sample_number_list <- unlist(sample_number_list)
  sample_number <- sample_number_list[(Sample_count+1):(Sample_count+length(methyl_bed_list))]
  
  cat("Creating the Megaframe \n")
  
  mylist <- c()
  experimental_design_df <- data.frame()
  for (i in 1:length(methyl_bed_list)){ #Iterate through methyl beds one by one
    tmpsampleData <- read.csv(methyl_bed_list[i], sep="\t", header=FALSE, nrows = 5)
    classes <- sapply(tmpsampleData, class)
    #replace some columns to null to delete them
    classes[c(3, 4, 5, 7, 8, 9)] <- "NULL"
    #import the bed file
    import_bedfile <- data.frame(purrr::map(methyl_bed_list[i], ~fread(.x, sep="\t", header=FALSE, colClasses = classes)))
    #call get_standard_methyl_bed function to clean up the bed file from each sample
    methyl_data <- get_standard_methyl_bed(Methyl_bed = import_bedfile, Sample_ID = sample_number[i], Methyl_call_type = Methyl_call_type, max_read_depth = max_read_depth )
    #getting the count of nrow for sanity checks
    raw_count <- nrow(methyl_data)
    methyl_data <- unique(methyl_data) #remove duplicates if any
    clean_count <- nrow(methyl_data)
    #QC
    if (raw_count==clean_count){
      cat("QC : No duplicates in ",methyl_bed_list[i] ,", proceeding \n")
    }
    else {
      cat("QC : Duplicates found in ",methyl_bed_list[i] ,", cleaning data before proceeding \n")
    }
    mylist[(length(mylist) + 1)] <- list(methyl_data) #append it to a list
    #get a list of alphabet codes and bed files - this will be saved in the experimental design starter
    Bedfile_comb <- data.frame(sample_number[i],methyl_bed_list[i])
    experimental_design_df <- rbind(experimental_design_df,Bedfile_comb)
    
  }
  write.table(experimental_design_df, paste(File_prefix, "Experimental_design_starter.csv",sep="_"), row.names=F, col.names = c("ID","Bedfile"), sep=",")
  
  cat("The experimental design file is now available in current directory! \n")
  
  #merge the methyl beds from diff samples into a signle large data frame
  combined_methyl_beds <-Reduce(function(x, y) full_join(x, y, by=c("Chromosome", "Position")), c(mylist) )
  
  #get Strand and CX columns to coalesce.
  Strands <- combined_methyl_beds %>% select(starts_with("Strand_")) %>% colnames()
  combined_methyl_beds$Strand <- do.call(dplyr::coalesce, combined_methyl_beds[Strands])
  
  Cxs <- combined_methyl_beds %>% select(starts_with("CX_")) %>% colnames()
  combined_methyl_beds$CX <- do.call(dplyr::coalesce, combined_methyl_beds[Cxs])
  
  #Remove unwanted columns
  combined_methyl_beds <- combined_methyl_beds  %>% select(-(starts_with("Strand_")), -starts_with(("CX_")))
  #Rearrange
  Megaframe <- combined_methyl_beds[,c(1:2,(ncol(combined_methyl_beds)-1),(ncol(combined_methyl_beds)),3:(ncol(combined_methyl_beds)-2) )]
  
  #sanity check - to ensure no NAs in Strand and CX columns after coalesce
  if ( (sum(is.na(Megaframe$Strand))==0 ) &
       sum(is.na(Megaframe$CX))==0 ) {
    cat('QC : Megaframe looks good \n')
  } else {
    cat('QC: Strand and CX should not have NAs, re-run the megaframe function \n')
  }
  
  write.table(Megaframe, paste(File_prefix, "MegaFrame.csv",sep="_"), row.names=F, sep=",")
  
  cat("Megaframe is now available in current directory and in the R-env! \n")
  
  
  
  #QC : Filter rows/sample with missing data
  Megaframe$NAs <- rowSums(is.na(Megaframe))
  #make a histogram
  
  cat("QC: The plot provides information about missing data that can be filtered out in the next step by using the filter_NAs parameter \n")
  QCplot <- suppressMessages(ggplot(Megaframe, aes(x=NAs/3))+geom_histogram(bins=30) +
                               labs(title = "Missing data per cytosine") +
                               xlab("Individual") +
                               ylab("Count of rows with missing data")) + theme_bw()
  
  print(QCplot)
  rm(QCplot)
  
  return (Megaframe)
}



#' add_zoom_coords
#'
#' A function to add Zoom codes based on the gene positions.
#' Currently the codes for this is as below
#' 1- Anything that is only between gene start and gene stop
#' 2- Anything that is between Adaptive start and Adaptive Stop
#' 0- Anything that doesn't fall within in the above - to ensure we don't include these in the DMR analysis.
#'
#' @param target (df) - Subset of ONT-methyl bed to positions pertaining to a single gene at a time
#' @param geneco_index (int) - This is index that refers to the row number of the geneco df, required to assign proper Alphabet codes.
#' @param gcoord_exist (Boolean) - This is to use the function only the gene_cord_df file has the location of the gene.
#' @param Gene_col (str) - Use this column to specifify wether to add Gene Names or Ids in the Zoomframe.
#' @return Zoomframe (df) - Similar to Megaframe except this includes more information on targets, positions zero'ed to ATG for each target and a few other information with an additional column that included zoom codes
#' @import tidyverse
#' @import stringr
#' @export

add_zoom_coords <- function(target, gene_cord_df, geneco_index, gcoord_exist=TRUE, Gene_col="Gene_name") {
  
  if(gcoord_exist==TRUE){
    for (i in 1:nrow(target)) {
      if (gene_cord_df[[Gene_col]][geneco_index]==target$Gene[i]){
        if((target$Position[i]>=Geneco$Low[geneco_index]) & (target$Position[i]<=Geneco$High[geneco_index]) ) {
          target$Zoom_co[i] <- 1 #Gene body region
          
        } else if ((target$Position[i]>=Geneco$Adapt_Low[geneco_index]) & (target$Position[i]<=Geneco$Adapt_High[geneco_index])) {
          target$Zoom_co[i] <- 2 #Adaptive sequence region
          
        } else {
          target$Zoom_co[i] <- 0 #Region beyond adaptve sequence
        }
      } else {
        print ("Gene Names don't match, please check the gene_coordinate file")
      }
    }
  }
  else {
    target$Zoom_co <- "NA"
  }
  return(target)
  
}



#' Create Zoomframe
#'
#' A function to create a single-merged data frame from individual methyl beds in the experiment
#'
#' @param gene_cord_df (df) - Data frame containing gene-coordinate info
#' @param MFrame (df) - Megaframe data from the previous function
#' @param File_prefix (Flexible str) - This is to add a prefix to all the files that get exported while running the function.
#' @param filter_NAs (int) - Select this parameter based on the histogram plot. This will filter out NAs based on per sample
#' @param target_info (Boolean) - This takes in TRUE or FALSE. Enter TRUE only if the megaframe contains target genes that need to be differentiated from non-targets.
#' @param gene_list (list) - Provide a list of target genes to distiguish from non-target genes within the Zoomframe. By deafult it will take in the All the genes from the gene coordinates file.
#' @inheritParams add_zoom_coords
#' @return Zoomframe (df) - Similar to Megaframe except this includes more information on targets, positions zero'ed to ATG for each target and a few other information.
#' @import tidyverse
#' @import stringr
#' @export

generate_zoomframe <- function(gene_cord_df, MFrame, Gene_col, target_info=TRUE, gene_list=gene_cord_df[[Gene_col]], File_prefix="") {
  
  
  cat("Creating the ZoomFrame! \n")
  
  #create an empty df()
  Final_gene_set <- data.frame()
  for (i in 1:nrow(gene_cord_df)){
    if( any(gene_cord_df$Chromosome[i]==MFrame$Chromosome) ){ #make sure the chromosomes match between Mframe and gene_cord_df
      if ( max(MFrame$Position)>=gene_cord_df$Adapt_Low[i] && (min(MFrame$Position)<=gene_cord_df$Adapt_High[i]) ) {
        Gene_subset <- MFrame[MFrame$Chromosome %in% gene_cord_df$Chromosome[i], ] #subset based on the gene
        Gene_subset <- Gene_subset %>% dplyr::filter(Position >= (gene_cord_df$Adapt_Low[i]) & Position <=(gene_cord_df$Adapt_High[i]) )
        Gene_subset$Gene <- gene_cord_df[[Gene_col]][i] #Add-in the gene/geneID names
        if (gene_cord_df$Strand[i]=="+") {
          Gene_subset$Zeroth_pos <- (Gene_subset$Position - gene_cord_df$Low[i]) #Computing the Zeroth position to center everything around ATG.
        }
        else if (gene_cord_df$Strand[i] == "-") {
          Gene_subset$Zeroth_pos <-  (gene_cord_df$High[i]- Gene_subset$Position) #reorienting the anti-sense genes
        }
        #Call the add_zoom_coords function to add Zoom_co-ordinates
        Target_df <- add_zoom_coords(target=Gene_subset, gene_cord_df=gene_cord_df, geneco_index=i, gcoord_exist=TRUE, Gene_col=Gene_col)
        Final_gene_set <- rbind(Final_gene_set,Target_df) #append it to a Final dataframe
      }
    }
    else {
      print(c(i,"Chromosomes don't match, check the gene_cord_df file "))
    }
  }
  Final_gene_set <- Final_gene_set[,c(1,(ncol(Final_gene_set)-2),2:4,5:(ncol(Final_gene_set)-3),ncol(Final_gene_set)-1,ncol(Final_gene_set))]
  
  cat("Zoomframe generated, Adding in target info column; Almost done!\n")
  
  #Clean up columns with NAs
  Meth_Unmeth <- Final_gene_set %>% select(starts_with("Meth"), starts_with("UnMeth")) %>% colnames()
  cat("Columns to change NAs -> 0s\n" , Meth_Unmeth)
  #convert NAs to 0s. Here we are not changing the Percent methylation column
  Final_gene_set[Meth_Unmeth][is.na(Final_gene_set[Meth_Unmeth])] <- 0
  
  if (target_info==TRUE){
    for (i in 1:nrow(Final_gene_set)){
      if (Final_gene_set$Gene[i] %in% gene_list){
        Final_gene_set$Target_info[i] <- "T"
      }
      else {
        Final_gene_set$Target_info[i] <- "NT"
      }
    }
  }
  
  
  write.table(Final_gene_set, paste(File_prefix, "ZoomFrame.csv",sep="_"), row.names=F, sep=",")
  cat("\nZoomframe is available in your current directory!")
  
  return(Final_gene_set)
}

#' Create Methylframe
#'
#' A function to create either Megaframe (A single data frame containing all samples bedMethyl data) or Zoomframe (A version of the Megaframe which includes the Gene Name, Position Zero'ed in on ATG of the respective gene and other columns to help with downstream analysis)
#' NOTE:
#' 1.  In both cases, Megaframe data will be exported in the current working directory and Zoomframe will only be exported if the Gene coordinates file is provided
#' 2. If this function is run using 'gene_info= TRUE' then it needs additional parameters like gene_cordinate_file,  Gene_column,target_info, gene_list etc to precoeed. If no such file is present, enter gene_info as false and
#'
#' @param methyl_bed_list (list) - ONT methyl bed filenames for each individual contained within the directory. This will just be a list of bedfile names. The input will be the "methyl_bed_list" vector that you create in the previous step.
#' @param Sample_count (int) - This is required to assign proper alphabet codes. If you need to include the samples from a previous round, then enter the total number of samples from the previous round here. Default is 0. By default alphabetizing starts with 'A'.
#' @param Methyl_call_type (str) - A string that includes information about the type of run. Currently this package works on Megalodon , DSP (DeepSignal Plant) and Bonito.
#' @param filter_NAs (int) - Select this parameter based on the histogram plot. This will filter out NAs based on per sample, by default this is 0.
#' @param gene_info (str) - This takes in a boolean variable. If the gene info- coordinates, gene name etc are present then make sure to have the gene-cordinates.csv file in the right format (as shown in the sample data on github).
#' @param gene_cordinate_file (str) - File containing gene-coordinate info. It is important to be in a specific format and should have but not limited to the following columns : Chromosome | Gene_Name | Low | High | Adapt_Low | Adapt_High
#' Low and High : These are Start and Stop cordinates of the gene. Low is always the lower coordinate which could be start for the positive stranded gene and stop for the negative stranded gene and vice versa.
#' Adapt_Low and Adapt_High : Cordinates for adaptive regions around the Lower and Higher co-ordinate of the gene respectively
#' @param File_prefix (Flexible str) - This is to add a prefix to all the files that get exported and saved to the working directory while running the function.
#' @inheritParams generate_megaframe
#' @inheritParams generate_zoomframe
#' @return Megaframe(df) or Zoomframe(df) - Clean data frame containing combined methyl bed information for every individual in the experiment.
#' @import tidyverse
#' @examples
#' # Basic usage for methyl_call_type
#' # 1. With gene coordinate file
#' # generate_methylframe(methyl_bed_list= <list_of_BedMethyl files>, gene_info = TRUE, gene_cordinate_file = <File with gene info>, Gene_col=<Gene Name/Gene ID>, target_info=TRUE, gene_list = <list of target genes>)
#' # 2. Without gene cooridnate file
#' # generate_methylframe(methyl_bed_list= <list_of_BedMethyl files>, gene_info = FALSE )
#' @export

generate_methylframe <-function(methyl_bed_list=All_methyl_beds, Sample_count = 0,
                                Methyl_call_type="Dorado", filter_NAs=0, max_read_depth=100,
                                gene_info = FALSE, gene_coordinate_file = NULL, Gene_column='',
                                target_info=FALSE, gene_list = gene_coordinate_file[[Gene_column]],
                                File_prefix="Sample")
{
  
  if (typeof(gene_info) != 'logical'){
    stop("gene_info variable needs to be a boolen variable. Default is False")
  }
  
  QC <- is.null(gene_coordinate_file)
  if (gene_info==TRUE & (QC==TRUE | Gene_column=='') ) {
    stop("gene_info is TRUE. Please provide gene-coordinates file, additional values such as Gene_column and re-run the function. Look into documentation for additional information \n")
  }
  
  if (gene_info==TRUE & target_info==FALSE ) {
    stop("gene_info is TRUE, this means the target_info also needs to be TRUE. \n")
  }
  
  Megaframe <- generate_megaframe(methyl_bed_list=methyl_bed_list, Sample_count = 0,
                                  Methyl_call_type=Methyl_call_type, max_read_depth=max_read_depth, File_prefix=File_prefix)
  
  cat('\n NOTE: Filtering NAs default is set to 0, See documentation for ideas on how to use the filter \n')
  
  Megaframe <- Megaframe[Megaframe$NAs<=(filter_NAs*3),]
  
  if (gene_info==TRUE) {
    
    Zoomframe <- generate_zoomframe(gene_cord_df=gene_coordinate_file, MFrame = Megaframe,
                                    Gene_col=Gene_column,
                                    target_info=FALSE, gene_list=gene_list ,
                                    File_prefix=File_prefix)
    
    return(Zoomframe)
  }
  else {
    
    #Duplicating the column for downstream analysis since the functions look for a Zeroth_pos column
    Megaframe$Zeroth_pos <- Megaframe$Position
    
    return(Megaframe)
    
  }
  
}




#' Boot_score
#' @description
#' Bootstrap analysis comparing target region with other regions of genome
#'
#' @param sound_score_obj (list) - Aggregated Changepoint Object created from sound_score
#' @param target_gene (str) - name of the target gene that was hit with oligos
#' @param target_start (str) - Start position relative to gene start that was hit with oligos
#' @param target_end (str) - Stop position relative to gene start that was hit with oligos
#' @param nboots (str) - Number of bootstrap replicates
#' @param scoring_col_name (str) - Name of the score column to run bootstrapping on (one of dmr_score or dmr_score2)
#' @return Boot_Obj (list) - File that contains information updated on target region info and bootstrap values
#' @export


boot_score<-function(sound_score_obj = NA, target_gene= NA, target_start=-1000, target_end=0, nboots=1000, scoring_col_name ="dmr_score", direction_DMR="positive"){
  #Make necessary objects
  rs<-sound_score_obj$region_summary
  ms<-sound_score_obj$methyl_summary
  target_rs<-rs[rs$Gene==target_gene,]
  
  #create data frame to be filled with results of bootstrapping
  boot_out <- data.frame(matrix(ncol = 5, nrow = nboots+1))
  
  #provide column names
  colnames(boot_out) <- c('Gene', 'CG_Score', 'CHG_Score', "CHH_Score", 'Target')
  
  #Calculate precision adjusted score for each change point region
  target_rs <- target_rs %>%
    group_by(cp_group) %>%
    mutate(distance_from_target= min(c(abs(Start-target_start), abs(Start-target_end),abs(Stop-target_start), abs(Stop-target_end))),
           distance_from_target=ifelse(Start<=target_start & Stop >= target_end | Start>=target_start & Stop <= target_end , 0, distance_from_target),
           adjusted_soundscore= ifelse(distance_from_target>10000,0,(1-(0.0001*distance_from_target))*!!as.name(scoring_col_name)))
  
  if(direction_DMR=="negative"){
    target_rs$adjusted_soundscore<-target_rs$adjusted_soundscore*(-1)
  }
  
  if(direction_DMR=="absolute"){
    target_rs$adjusted_soundscore<-abs(target_rs$adjusted_soundscore)
  }
  
  #Find strongest DMR around target gene
  boot_out$CG_Score[1]<-max(target_rs$adjusted_soundscore[target_rs$CX=="CG"])
  boot_out$CHG_Score[1]<-max(target_rs$adjusted_soundscore[target_rs$CX=="CHG"])
  boot_out$CHH_Score[1]<-max(target_rs$adjusted_soundscore[target_rs$CX=="CHH"])
  
  # Write info to Boot_Object
  boot_out$Target[1]<- 1
  boot_out$Gene[1]<-target_gene
  boot_out$Start_CG[1]<-target_rs$Start[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CG"]) & target_rs$CX=="CG"]
  boot_out$Stop_CG[1]<-target_rs$Stop[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CG"]) & target_rs$CX=="CG"]
  boot_out$Start_CHG[1]<-target_rs$Start[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CHG"]) & target_rs$CX=="CHG"]
  boot_out$Stop_CHG[1]<-target_rs$Stop[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CHG"]) & target_rs$CX=="CHG"]
  boot_out$Start_CHH[1]<-target_rs$Start[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CHH"]) & target_rs$CX=="CHH"]
  boot_out$Stop_CHH[1]<-target_rs$Stop[target_rs$adjusted_soundscore==max(target_rs$adjusted_soundscore[target_rs$CX=="CHH"]) & target_rs$CX=="CHH"]
  
  # Make non-target frame
  nontarget_ms<-ms[ms$Gene!=target_gene,]
  
  # Create lookup table with info for each gene
  nontarget_region_lookup <- nontarget_ms %>%
    group_by(Gene) %>%
    summarise(Start=min(Zeroth_pos), Stop=max(Zeroth_pos))
  
  # Determine if gene is far enough from end of contig to be used for bootstrapping
  nontarget_ms$far_enough<-FALSE
  for(i in 1:nrow(nontarget_region_lookup)){
    nontarget_ms$far_enough[nontarget_ms$Gene == nontarget_region_lookup$Gene[i]] <- pmin(abs(nontarget_ms$Zeroth_pos[nontarget_ms$Gene == nontarget_region_lookup$Gene[i]]-nontarget_region_lookup$Start[i]), abs(nontarget_ms$Zeroth_pos[nontarget_ms$Gene == nontarget_region_lookup$Gene[i]]-nontarget_region_lookup$Stop[i])) > 10000
  }
  #subset to only include these rows
  nontarget_ms_good_distance<-nontarget_ms[nontarget_ms$far_enough==TRUE,]
  
  # sample from rows
  boot_positions<-sample_n(nontarget_ms_good_distance, nboots, replace=TRUE)
  
  # run bootstrapping
  for(i in 1:nrow(boot_positions)){
    boot_rs <- rs[rs$Gene==boot_positions$Gene[i],]
    target_start <- boot_positions$Zeroth_pos[i]
    target_end <- boot_positions$Zeroth_pos[i]
    boot_rs <- boot_rs %>%
      group_by(cp_group) %>%
      mutate(distance_from_target= min(c(abs(Start-target_start), abs(Start-target_end),abs(Stop-target_start), abs(Stop-target_end))),
             distance_from_target=ifelse(Start<=target_start & Stop >= target_end | Start>=target_start & Stop <= target_end , 0, distance_from_target),
             adjusted_soundscore= ifelse(distance_from_target>10000,0,(1-(0.0001*distance_from_target))*!!as.name(scoring_col_name )))
    
    # Write info to Boot_Object
    
    if(direction_DMR=="negative"){
      boot_rs$adjusted_soundscore<-boot_rs$adjusted_soundscore*(-1)
    }
    
    if(direction_DMR=="absolute"){
      boot_rs$adjusted_soundscore<-abs(boot_rs$adjusted_soundscore)
    }
    
    
    boot_out$CG_Score[i+1]<-max(boot_rs$adjusted_soundscore[boot_rs$CX=="CG"])
    boot_out$CHG_Score[i+1]<-max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHG"])
    boot_out$CHH_Score[i+1]<-max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHH"])
    boot_out$Target[i+1]<-0
    boot_out$Gene[i+1]<-boot_positions$Gene[i]
    boot_out$Start_CG[i+1]<-boot_rs$Start[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CG"]) & boot_rs$CX=="CG"][1]
    boot_out$Stop_CG[i+1]<-boot_rs$Stop[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CG"]) & boot_rs$CX=="CG"][1]
    boot_out$Start_CHG[i+1]<-boot_rs$Start[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHG"]) & boot_rs$CX=="CHG"][1]
    boot_out$Stop_CHG[i+1]<-boot_rs$Stop[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHG"]) & boot_rs$CX=="CHG"][1]
    boot_out$Start_CHH[i+1]<-boot_rs$Start[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHH"]) & boot_rs$CX=="CHH"][1]
    boot_out$Stop_CHH[i+1]<-boot_rs$Stop[boot_rs$adjusted_soundscore==max(boot_rs$adjusted_soundscore[boot_rs$CX=="CHH"]) & boot_rs$CX=="CHH"][1]
    
  }
  
  bo_CG <- boot_out[!duplicated(boot_out[c(1,2)]),]
  bo_CHG <- boot_out[!duplicated(boot_out[c(1,3)]),]
  bo_CHH <- boot_out[!duplicated(boot_out[c(1,4)]),]
  
  bo_CG_order <- bo_CG[order(-bo_CG$CG_Score),]
  bo_CHG_order <- bo_CHG[order(-bo_CHG$CHG_Score),]
  bo_CHH_order <- bo_CHH[order(-bo_CHH$CHH_Score),]
  
  bo_CG_order$rank<-c(seq(1, nrow(bo_CG_order), by=1))
  bo_CHG_order$rank<-c(seq(1, nrow(bo_CHG_order), by=1))
  bo_CHH_order$rank<-c(seq(1, nrow(bo_CHH_order), by=1))
  
  print(paste("Precision Adjusted CG DMR score of:", round(boot_out$CG_Score,3)[1], " For a CG bootstrap p-value of: ", (bo_CG_order[bo_CG_order$Target==1,]$rank)/nrow(bo_CG_order)))
  print(paste("Precision Adjusted CHG DMR score of:", round(boot_out$CHG_Score,3)[1], " For a CHG bootstrap p-value of: ", (bo_CHG_order[bo_CHG_order$Target==1,]$rank)/nrow(bo_CHG_order)))
  print(paste("Precision Adjusted CHH DMR score of:", round(boot_out$CHH_Score,3)[1], " For a CHH bootstrap p-value of: ", (bo_CHH_order[bo_CHH_order$Target==1,]$rank)/nrow(bo_CHH_order)))
  
  print(rbind(target_rs[target_rs$adjusted_soundscore==boot_out$CG_Score[1] & target_rs$adjusted_soundscore!=0,], target_rs[target_rs$adjusted_soundscore==boot_out$CHG_Score[1] & target_rs$adjusted_soundscore!=0,], target_rs[target_rs$adjusted_soundscore==boot_out$CHH_Score[1]& target_rs$adjusted_soundscore!=0,]))
  
  print(paste("Final Bootstrap Adjusted CG DMR Score:", round(bo_CG[1,2]*((bo_CG[1,2]-mean(bo_CG[-1,2], na.rm = TRUE))/sd(bo_CG[,2], na.rm=TRUE)),2)))
  print(paste("Final Bootstrap Adjusted CHG DMR Score:", round(bo_CHG[1,3]*((bo_CHG[1,3]-mean(bo_CHG[-1,3], na.rm = TRUE))/sd(bo_CHG[,3], na.rm = TRUE)),2)))
  print(paste("Final Bootstrap Adjusted CHH DMR Score:", round(bo_CHH[1,4]*((bo_CHH[1,4]-mean(bo_CHH[-1,4], na.rm = TRUE))/sd(bo_CHH[,4], na.rm = TRUE)),2)))
  print(paste("Less than 1: Nothing there"))
  print(paste("1-2: Subtle shifts in methylation"))
  print(paste("2-3: Moderate methylation shifts near oligo treatment"))
  print(paste("3-5: Significant methylation shifts near oligo treatment"))
  print(paste("5+: Very Strong evidencce of DMR associated with oligo treatment"))
  
  
  Boot_Obj <- list(boot_out, target_rs)
  names(Boot_Obj) <- c("bootstrap_scores", "target_rs")
  return(Boot_Obj)
  
}
