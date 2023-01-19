#' Create Gene Percent Data Frames
#'
#' A function to create different subsets of the LongPercent df depending on the
#' column of interest
#'
#' @param LongPercent data frame containing the methylation data as percents in
#' long format using `dcast()`
#' @param x A string representing the column to compare to Gene * Zeroth_pos
#' @param function_name an aggregation function such as mean, sd, or var
#' @return dcast_output a subset of the LongPercent dataframe
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


#' Clean Input Data
#'
#' A function to clean the input data of methylation data. It filters the data
#' frame to contain only the columns of interest passed in as an argument and
#' removes rows where the mean methylation is 0.
#'
#' @param ZoomFrame data frame containing the input methylation data
#' @param Exp_ID data frame containing the experimental design
#' @param colnames_of_interest *optional* list of strings of the columns to keep
#' in the analysis
#' @return ZoomFrame_filtered a data frame subset of `ZoomFrame` containing only
#' the columns passed in to `colnames_of_interest`
#' @import tidyverse
#' @export

clean_data <- function(ZoomFrame = dataframe,
                       Exp_ID = dataframe,
                       colnames_of_interest = c('Chromosome', 'Gene', 'Position',
                                               'Strand', 'CX', 'Zeroth_pos',
                                               'Individual')) {
  # Clean the data
  print('Step 1: removing rows that contain 0 methylation')
  Inputpers <- dplyr::select(ZoomFrame, starts_with("Per"))
  Inputpers$MeanMeth <- rowMeans(Inputpers, na.rm=TRUE)
  ZoomFrame_filtered <- ZoomFrame[Inputpers$MeanMeth != 0,]
  print('Step 1: complete')

  #QC Make sure the correct columns are present
  print('Step 2: checking for missing columns in experimental id')
  for (col in colnames_of_interest) {
    if (!col %in% colnames(ZoomFrame_filtered)){
      if (col != 'Individual' & col != 'Plant') {
        print(paste(col, 'not found in the ZoomFrame_filtered'))
      }
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
  if (!'Plant' %in% colnames(Exp_ID) & !'Individual' %in% colnames(Exp_ID)) {
    stop('No Plant or Individual column name found in Experimental ID. One is necessary to continue')
  }
  if (!'Plant' %in% colnames(Exp_ID)) {
    print('Plant column not found in the Experimental ID, using Individual')
    Exp_ID$Plant <- Exp_ID$Individual
  }
  if (!'Individual' %in% colnames(Exp_ID)) {
    print('Individual column not found in the Experimental ID, using Plant')
    Exp_ID$Individual <- Exp_ID$Plant
  }
  if (!'Individual_Name' %in% colnames(Exp_ID)) {
    print('Individual_Name column not found in the Experimental_ID, using Plant')
    Exp_ID$Individual_Name <- Exp_ID$Plant
  }
  print('Step 2: complete')

  # Convert the experimental ID 'Plant' column to character
  Exp_ID$Plant <- as.character(Exp_ID$Plant)

  #We want to order the input frame in such a way that it will be easy to recreate analysis
  print('Step 3: reordering ZoomFrame_filtered')
  ZoomFrame_filtered <- ZoomFrame_filtered[order(ZoomFrame_filtered[,'Gene'], ZoomFrame_filtered[,'Zeroth_pos']), ]
  print('Step 3: complete')

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

  print('Step 3: complete')
  #Replacing read depth of NA with 0
  LongUnMeth$RD[is.na(LongUnMeth$RD)] = 0
  LongMeth$RD[is.na(LongMeth$RD)] = 0

  #Creating column of total read depth
  LongMeth$total_RD <- LongMeth$RD + LongUnMeth$RD
  LongMeth$Individual <- gsub("Meth_","", LongMeth$Individual, perl = T)

  #QC to check that LongPercent contains the correct number of rows
  if (nrow(LongPercent) == nrow(ZoomFrame_filtered) * nrow(Exp_ID)) {
    print('LongPercent contains the expected number of rows')
  } else {
    print(paste('LongPercent contains', nrow(LongPercent), 'rows. Expected:',
                (nrow(Zoom_Frame_filtered) * nrow(Exp_ID)), 'rows'))
  }
  print('Step 4: complete')

  print('Step 4: merging long files with Exp_ID to annotate')
  # Merge Long files with Exp_ID to annotate
  print('Step 5: annotating long files with Experimental ID')
  LongPercent <- dplyr::left_join(LongPercent, Exp_ID, by = c('Individual' = 'ID'))
  LongMeth <- dplyr::left_join(LongMeth, Exp_ID, by = c('Individual' = 'ID'))
  print('Step 5: complete')

  # Aggregate
  print('Step 6: aggregating by plant')
  LongPercent <- LongPercent %>%
    group_by(Gene, Zeroth_pos, Plant, Position, CX, Strand, Group, Chromosome) %>%
    summarize(Percent = mean(Percent, na.rm = T))

  LongMeth <- LongMeth %>%
    group_by(Gene, Zeroth_pos, Plant, Position, CX, Strand, Group, Chromosome) %>%
    summarize(total_RD = mean(total_RD, na.rm = T))
  LongMeth <- LongMeth[,c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                          'Zeroth_pos', 'Plant', 'total_RD', 'Group')]
  print('Step 6: complete')

  print('Step 5: complete')

  # Save all the outputs
  out <- list()
  out$ZoomFrame_filtered <- ZoomFrame_filtered
  out$LongPercent <- LongPercent
  out$LonUnMeth <- LongUnMeth
  out$LongMeth <- LongMeth
  out$Exp_ID <- Exp_ID
  out$Inputpers <- Inputpers

  return(out)
}


#' Subset Data Frames to Contain Only Columns of Interest
#'
#' Keep only the necessary columns from a given data frame
#'
#' @param df the data frame to subset
#' @inheritParams clean_data
#' @return subsetted_df data frame with only the columns of interest
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
#' @param data the `ZoomFrame_filtered`
#' @param starts_with_cols the string pattern that the columns start with
#' @inheritParams tidyr::pivot_longer
#' @inheritParams clean_data
#' @return pivoted_and_subsetted_df is the cleaned data frame that has been
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
#' @param Exp_ID_Treated a data frame containing the experimental design
#' @param OF the Output_Frame
#' @param GenePercentPlant summary dataframe built from `create_gene_percent_x()`
#' containing % methylation for each Plant
#' @param GenePercentGroup summary dataframe built from `create_gene_percent_x()`
#' containing % methylation for groups
#' @param GeneDepthPlant summary dataframe containing read depth for individuals
#' @param control a string representing the control variable in GenePercentGroup
#' @return OF the updated Output_Frame with 3 new columns for each individual in
#' the `Exp_IF_Treated` parameter
#' @export

create_cols_for_individuals <- function(Exp_ID_Treated,
                                        OF,
                                        GenePercentPlant,
                                        GeneDepthPlant,
                                        GenePercentGroup,
                                        control = 'C') {
  for(id in unique(Exp_ID_Treated$Plant)) {
    OF$NewZ <- 0
    # Name this column using trickery
    colnames(OF)[ncol(OF)] <- paste(id, "Z", sep="_")
    # Create new column that is difference in percent methylation between individual "i" and the control average
    OF[,ncol(OF) + 1] <- GenePercentPlant[,id] - GenePercentGroup[[control]]
    # Name this column
    names(OF)[ncol(OF)] <- paste(id, "MethChange", sep="_")
    # Create a column that is read depth for that individual
    OF[,ncol(OF) + 1] <- GeneDepthPlant[,id]
    # Name this column
    names(OF)[ncol(OF)] <- paste(id, "RD", sep="_")
  }
  return(OF)
}


#' Find the Difference Between Two Columns
#'
#' This function is intended to find the difference between treatment and
#' control groups, but can be used to find the difference between any two columns
#' in the given data
#'
#' @param data a data frame, usually the GenePercentGroup with treatment and
#' control columns
#' @param treatment a string of the name of the treatment column
#' @param control a string of the name of the control column
#' @return treat_v_control a **vector** containing the difference between
#' treatment and control
#' @export

find_col_diff <- function(data,
                          treatment = "treatment",
                          control = "control") {
  treat_v_control <- data[[treatment]] - data[[control]]
  return(treat_v_control)
}


#' Create Output Frame
#'
#' This function creates the Output_Frame that the group and individual DMR
#' functions rely on for their analysis. It includes some quality control checks
#' to make sure the output is in the correct format.
#'
#' @inheritParams create_cols_for_individuals
#' @param data the `ZoomFrame_filtered` produced in the clean_data function
#'
#' @export
create_output_frame <- function(Exp_ID_Treated, data, GenePercentPlant,
                                GeneDepthPlant, GenePercentGroup,
                                control = 'C') {
  #Create the beginnings of the Output statistic data frame
  Output_Frame <- data[,colnames_of_interest[-length(colnames_of_interest)]]

  # create 3 new cols for each individual
  Output_Frame <- create_cols_for_individuals(Exp_ID_Treated, Output_Frame,
                                              GenePercentPlant, GeneDepthPlant,
                                              GenePercentGroup, control)

  # QC
  # Checking OF created the three new cols for each individual
  if (ncol(Output_Frame) == (length(colnames_of_interest) - 1) + 3 * length(unique(Exp_ID_Treated$Plant))) {
    print('Number of columns in Output_Frame is correct')
  } else {
    print('The number of columns in Output_Frame is incorrect. Double check there are no duplicates')
    print(paste('Expected:', (length(colnames_of_interest) - 1) + 3 * length(unique(Exp_ID_Treated$Plant)),
                'Found:', ncol(Output_Frame)))
  }

  #Add in any summary statistic columns, such as this one
  Output_Frame <- cbind(Output_Frame,GenePercentPlant[,3:ncol(GenePercentPlant)])
  Output_Frame$Treat_V_Control <- find_col_diff(GenePercentGroup, 'T', 'C')
  Output_Frame$Control <- GenePercentGroup$C

  return(Output_Frame)
}

 #' Filter Experimental ID
 #'
 #' This function is used to filter the experimental id based on a condition.
 #' Currently it is only useful for filtering categorical data, but in time
 #' will be used for all types of data.
 #'
 #' @param experimental_id a data frame containing the experimental design information
 #' @param col_name the name of the column by which to filter
 #' @param list_to_save a list of strings representing the values to save
 #' @return sub_ID the filtered experimental_id
 #' @examples
 #' # Generate data
 #' experimental_id <- data.frame(ID = c('A', 'B', 'C', 'D'),
 #'                               val = c(1, 2, 3, 4))
 #'
 #' # Filter the data
 #' create_sub_ID(experimental_id, 'ID', c('A', 'C'))
 #'
 #' @export

 create_sub_ID <- function(experimental_id, col_name = 'ID',
                           list_to_save = c('A', 'B', 'C')) {
  sub_ID <- experimental_id[experimental_id[[col_name]] %in% list_to_save,]
  return(sub_ID)
}


#' Create Fixed Effects
#'
#' A function to create the string combining fixed effects that will be passed
#' into `create_function()`
#'
#' @param fixed a list of strings containing the fixed effects elements
#' @return fixed_effects a string containing the fixed effects properly formatted
#' to be passed into `create_function()`
#' @examples
#' # create multiple independent fixed effects
#' create_fixed_effects(c('Group', 'Individual'))
#'
#' # Create fixed effects with interactions
#' create_fixed_effects(c('Group * Individual'))
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
#' @param random a list of strings containing the random effects elements
#' @return random_effects a string of the random effects properly formatted
#' @examples
#' # Create independent random effects
#' create_random_effects(c("Group", "Individual"))
#'
#' # Create random effects with an interaction
#' create_random_effects(c("Group * Individual"))
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
#' @param fixed a string containing the fixed effects variables. Note: this
#' value **cannot** be 'ID'. ID is used to merge data together downstream for
#' running the model so including the value in the effects will break the model.
#' @param random a string containing the random effects variables. Note: this
#' value **cannot** be 'ID'. ID is used to merge data together downstream for
#' running the model so including the value in the effects will break the model.
#' @return effects_formula a formula of the mixed effects
#' @examples
#' # Create formula of independent and random effects
#' create_formula(fixed = c('Group'), random = c('Individual'))
#'
#' create_formula(fixed = c('Group', 'Individual'), random = c('Plant' * 'Gene'))
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
#' @param i an integer for the row number of the dataframe
#' @param Output_Frame the data frame containing all the information
#' @param model_summary the data frame containing the summary statistics from
#' the model
#' @param ind_name (optional) containing the name of the z score column for
#' individual DMR analysis
#' @return Output_Frame the data frame with the proper row updated to include
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
#' @param LM containing the information to put in the model, usually this is
#' in a "long" format
#' @param i an integer representing the row number
#' @param Output_Frame the data frame containing all the summary information
#' @param formula the formula to use in the model
#' @param optimizer_func a string representing the optimizer function to use in
#' the model
#' @return Output_Frame the data frame containing updated summary information
#' @import lme4
#' @export
#'

run_binomial <- function(LM, i = int, Output_Frame, formula,
                      optimizer_func = 'optimizer') {
  binom_model <- glmer(formula, data=LM, family = binomial,
                       glmerControl(check.conv.grad = .makeCC(action = "stop",
                                                        tol = 2e-3, relTol = NULL),
                              optimizer=optimizer_func, optCtrl = list(maxfun = 500000)))
  # Outputs of this model
  model_summary <- as.data.frame(summary(binom_model)$coefficients)
  Output_Frame <- save_model_summary(i, Output_Frame, model_summary)

  return(Output_Frame)
}


#' Run Model
#'
#' Function to run the correct model
#'
#' @param data containing the information to put in the model, usually this is
#' in "long" format
#' @inheritParams run_binomial
#' @param individual_name_z (optional) the name of the z column when running
#' individual DMR analysis
#' @return Output_Frame updated with the summary data from the model
#' @export
#'

run_model <- function(data, i, Output_Frame, formula, model_type,
                      individual_name_z = ''){
  if (model_type == 'binomial') {
    tryCatch({
      Output_Frame <- run_binomial(data, i, Output_Frame, formula, 'bobyqa')

      # If that model didn't converge, it tried again with a different optimizer, allows ~20% more model convergence
    },error = function(e){tryCatch({ print(paste(i, "No bobyqa Converge, trying Nelder"))
      # Run the model with Nelder_Mead optimizer
      Output_Frame <- run_binomial(data, i, Output_Frame, formula, 'Nelder_Mead')

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
      model_summary <- as.data.frame(summary(beta_binomial)$coefficients$cond)
      Output_Frame <- save_model_summary(i, Output_Frame, model_summary,
                                          individual_name_z)
    }, error=function(e){
      #print(paste(i, "No Converge"))
      paste(i, 'No Converge')
      })
  } else {
    print('Please choose a model type of "binomial" or "beta-binomial".')
    Output_Frame <- Output_Frame
  }
  return(Output_Frame)
}

#' Group DMR Analysis
#'
#' To run a binomial model to compare the methylation between groups
#'
#' @param Output_Frame data frame containing the read depth and methylation
#' change information
#' @param ZoomFrame_filtered data frame containing the percent methylation information
#' @param Exp_ID data frame of the experimental design
#' @inheritParams create_formula
#' @param reads_threshold integer representing the number of reads that are each
#' methylated and unmethylated. This is important since data containing only one
#' methylated read is unlikely to provide statistical power to our analysis. A
#' value of 3 here means that the read depth is *at least* the depth provided by
#' the threshold. Both the methylated and unmethylated samples need a read depth
#' greater than or equal to this threshold in order to be considered for the model.
#' @inheritParams subset_cols
#' @return Output_Frame data frame containing the summary statistics from the model
#' @export

group_DMR <- function(Output_Frame, ZoomFrame_filtered, Exp_ID, fixed = c('Group'),
                      random = c('Plant'), reads_threshold = 3, model = 'binomial',
                      colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                               'Zeroth_pos', 'Individual')) {
  #  Create the model formula first
  formula <- create_formula(fixed, random)
  print(formula)

  #  Get the number of columns in the Output_Frame dataframe
  original_OF_col_number <- ncol(Output_Frame)

  # The modelling here is the most "delicate" part of the operation.  Options include:
  # (A) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment
  # (B) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment + Generation
  # (C) cbind(Meth, UnMeth) ~ (1|Plant) + Treatment + Generation +Treatment*Generation
  # (D) cbind(Meth, UnMeth) ~ (1|Plant) + Gene_Expression
  # (E) cbind(Meth, UnMeth) ~ (1|Plant) + Phenotype

  # Loop to run groupwise analysis for each BP
  for(i in 1:nrow(Output_Frame)){
    # Make long version of input frame for the i'th cytosine row
    LM <- pivot_and_subset(ZoomFrame_filtered[i,], 'Meth', 'Meth',
                           colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                      'Zeroth_pos', 'Individual'))
    LM$Individual <- gsub("Meth_","", LM$Individual, perl = T)
    # For rows with at least X(3) methylated reads across all individuals, move forward
    if(sum(as.numeric(LM$Meth), na.rm=TRUE) >= reads_threshold){
      # Do the same thing for unmethylated reads
      LUM <- pivot_and_subset(ZoomFrame_filtered[i,], 'UnMeth', 'UnMeth',
                              colnames_of_interest = c('Chromosome', 'Gene', 'Position', 'Strand', 'CX',
                                                       'Zeroth_pos', 'Individual'))
      LM <- cbind(LM,LUM[,ncol(LUM)])
      # Merge with Exp_ID to allow for DML testing for that base
      LM <- merge(LM, Exp_ID, by.x="Individual",by.y="ID", all=FALSE)

      # Has to have Y unmethylated reads across all individuials
      if(sum(as.numeric(LM$UnMeth), na.rm=TRUE) >= reads_threshold){
        # Run the model and save the output
        Output_Frame <- run_model(LM, i, Output_Frame, formula, model)
      }
    }
  }
  #  Replace NA values in these columns with 0s
  Output_Frame[,original_OF_col_number:ncol(Output_Frame)][is.na(Output_Frame[,original_OF_col_number:ncol(Output_Frame)])] = 0

  return(Output_Frame)
}



#' Individual DMR Analysis
#'
#' Purpose: to run individual DMR analysis
#'
#' @inheritParams group_DMR
#' @param control a string to represent the control variable
#' @return Output_Frame a data frame with additional columns
#' @export

individual_DMR <- function(Output_Frame, ZoomFrame_filtered, Exp_ID_Treated,
                           fixed = c('Group'), random = c('Individual'),
                           reads_threshold = 3, control = 'C', model = 'beta-binomial',
                           colnames_of_interest = c('Chromosome', 'Gene', 'Position',
                                                    'Strand', 'CX', 'Zeroth_pos', 'Individual')) {

  # Obtain the formula
  formula <- create_formula(fixed, random)
  print(formula)

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
      LM <- merge(LM, Exp_ID, by.x="Individual", by.y="ID", all=FALSE)

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
#' @param analysis_type either "individual" or "group"
#' @return Output_Frame
#' @export

DMR <- function(Output_Frame, ZoomFrame_filtered, Exp_ID, fixed = c('Group'),
                random = c('Plant'), colnames_of_interest, reads_threshold = 3,
                model, control = '', analysis_type) {
  if (tolower(analysis_type) == 'group') {
    Output_Frame = group_DMR(Output_Frame, ZoomFrame_filtered, Exp_ID,
                             fixed = fixed,random = random, colnames_of_interest,
                             reads_threshold = reads_threshold, model = model)
  } else if (tolower(analysis_type) == 'individual') {
    Output_Frame = individual_DMR(Output_Frame, ZoomFrame_filtered, Exp_ID,
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
#' @param DMR_output the DMR Analysis output dataframe
#' @return z_cols a list of strings representing the potential column names that
#' are suggested to be used in the changepoint analysis.
#' @export

find_changepoint_col_options <- function(DMR_output, OF = Output_Frame) {
  #added_cols <- colnames(DMR_output[c((ncol(OF) + 1):ncol(DMR_output))])
  z_cols <- colnames(DMR_output)[grepl('Z_', colnames(DMR_output))]

  print(z_cols)
}

#' Store Changepoint Information
#'
#' Add the important changepoint summary info to the dataframe from
#' which the changepoints were found
#'
#' @param whole_df data frame containing all the information for one gene
#' @param changepoint the changepoint summary
#' @param col_name the name of the column upon which the changepoint analysis
#' is being conducted
#' @return whole_df data frame containing information for one gene including the
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
#' @param data the dataframe that has been subsetted for the gene and cytosine
#' context
#' @param penalty the penalty value to include for the changepoint analysis
#' @return an S4 object containing the changepoint information
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
#' @param data containing the cytosines with their percent methylation
#' @param changepoint_obj the S4 object created from finding the changepoints
#' @param gene_name the name of the gene
#' @param penalty_val the penalty value being used
#' @param cyt_context the cytosine context
#' @param z_col the column upon which the changepoints are calculated
#' @return changepoint_plot the plot of the changepoints
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
      group_data = data[data[[paste0('MethGroup_', z_col)]] == i, 'Position']
      segment_data[i, 'x'] = min(group_data)
      segment_data[i, 'xend'] = max(group_data)
      segment_data[i, 'y'] = changepoint_obj@param.est$mean[i]
      segment_data[i, 'yend'] = changepoint_obj@param.est$mean[i]
    }

    # Create basic plot
    plot <- data %>%
      ggplot(aes(x = Position, y = !!sym(z_col))) +
      geom_line(size = 0.5) +
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
#' @param whole_df containing the information such as genes, z scores, etc
#' @param CG_penalty,CHG_penalty,CHH_penalty penalty values for the
#' change point analysis. The higher the value the fewer changepoints that will
#' be created. The lower the value the more changepoints that will be created.
#' @param z_col the name of the column to run the changepoint analysis on. This
#' can be any column, but for DMR analysis we recommend using the z scores for
#' a fixed effect variable.
#' @return everything data frame containing the mean changepoint value for the
#' `z_col` column
#' @import tibble
#' @export

changepoint_analysis <- function(whole_df,
                                 CG_penalty = int,
                                 CHG_penalty = int,
                                 CHH_penalty = int,
                                 z_col = 'column') {
  everything <- tibble::as_tibble(matrix(ncol = 8))
  colnames(everything) = c('index', 'CX', 'Position', 'Gene', paste0('MeanMeth_',z_col),
                           paste0('MeanMethRegion_',z_col), paste0('MethRegionLength_',z_col),
                           paste('MethylGroup',z_col))
  #  remove the first row
  everything <- everything[-1,]
  genes <- unique(whole_df$Gene)

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
    plot.CG <- plot_changepoints(gene_df.CG, x.CG, gene, CG_penalty, 'CG', z_col)
    plot.CHG <- plot_changepoints(gene_df.CHG, x.CHG, gene, CHG_penalty, 'CHG', z_col)
    plot.CHH <- plot_changepoints(gene_df.CHH, x.CHH, gene, CHH_penalty, 'CHH', z_col)
    print(plot.CG)
    print(plot.CHG)
    print(plot.CHH)

    # Combine everything into the everything dataframe (rename)
    everything <- rbind(everything, gene_df.CG, gene_df.CHG, gene_df.CHH)

    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(everything)
}

#' tidyMeg - Tidy data coming from ONT methylation calls for each individual
#'
#' A function to create a clean data frame for every individual in the experimental design without having to re run for every individual separately.
#'
#' @param Megfile Data frame containing the ONT methylation calls in a bed file for each individual
#' @param sample_ID A string that takes in Alphabet code to assign to every individual sample in the experiment
#' @return Megfile_sub Clean data frame methyl file for every individual with Meth, Unmeth and Per_Meth columns
#' @import tidyverse
#' @import stringr
#' @export

tidy <-function(Megfile = ONT_bedfile, sample_ID = "E", call_type="DSP") {

  # Extract columns of interest based on which process was run
  if (call_type=='DSP' | call_type=='Bonito'){
    Megfile_sub <- Megfile[,c(1,2,6,10,11,12)]
  }
  else if(call_type=="Megalodon"){
    Megfile_sub <- Megfile[,c(1,2,6,10,11,13)]
  }
  colnames(Megfile_sub) <- c("Chromosome","Position",paste("Strand", sample_ID, sep="_"),"Tot_reads", paste("PerMeth", sample_ID, sep="_") , paste("CX", sample_ID, sep="_")) #Assign column names.
  #calculate Meth and Unmeth reads from total reads. This is necessary for downstream DMR analysis.
  Megfile_sub[[paste("Meth", sample_ID, sep="_")]] <- round( (Megfile_sub[[paste("PerMeth", sample_ID, sep="_")]] * Megfile_sub$Tot_reads)/100 )
  Megfile_sub[[paste("UnMeth", sample_ID, sep="_")]] <- (Megfile_sub$Tot_reads - Megfile_sub[[paste("Meth", sample_ID, sep="_")]])
  Megfile_sub <- Megfile_sub %>% select(-Tot_reads)

  return(Megfile_sub)
}


#' GenerateMegF
#'
#' A function to create a single-merged data frame from individual methyl beds in the experiment
#'
#' @param Bedfiles ONT methyl bed files for each individual contained within the directory
#' @param N_prev_sample This is required to assign proper Alphabet codes. If you need to include te samples from a previous round, then mention the total number of samples here. Default is 0. In this case Alphabetizing starts with 'A'
#' @param  project_info This is to add a prefix to all the files that get exported while running the function.
#' @return Megaframe Clean data frame containing methyl bed information for every individual in the experiment
#' @import tidyverse
#' @import stringr
#' @export



GenerateMegF <- function(Bedfiles=All_beds, call_type="DSP", project_info=""){


  Alphabet_code <- c()
  for(i in 1:10) {
    code <- paste(LETTERS,LETTERS[i],sep="")
    Alphabet_code[(length(Alphabet_code) + 1)] <- list(code)
  }

  Alphabet_code <- unlist(Alphabet_code)
  IDcat <- c( LETTERS[1:length(All_beds)] )

  
  cat("Creating the Megaframe; sit tight!\n")
  mylist <- c()
  Exp_Id <- data.frame()
  for (i in 1:length(Bedfiles)){ 
    for (j in 1:length(IDcat)){
      if(i==j){
        #import the bed file
        Methyl_data <- data.frame(purrr::map(Bedfiles[i], ~read.csv(.x, sep="\t", header=FALSE)))
        #call tidy cats for each sample
        clean_data <- tidy(Megfile = Methyl_data, sample_ID = IDcat[j], call_type="DSP" )
        raw <- nrow(tidy)
        tidy <- unique(tidy)
        clean <- nrow(tidy)
        if (raw==clean){
          cat("No duplicates in ",Bedfiles[i] ,", proceeding \n")
        }
          else {
            cat("Duplicates found in ",Bedfiles[i] ,", cleaning data before proceeding \n")
          }
        mylist[(length(mylist) + 1)] <- list(tidy) #append it to a list
        Bedfile_comb <- data.frame(IDcat[j],Bedfiles[i])
        Exp_Id <- rbind(Exp_Id,Bedfile_comb)
        #print(c(All_beds[i], IDcat[j]) )
      }
      else {
        
      }
    }
  }
  write.table(Exp_Id, paste(project_info, "Experimental_design_starter.csv",sep="_"), row.names=F, col.names = c("ID","Library"), sep=",")

  cat("The matrix file is now available in current directory!\n")

  #merge the megalodon output from diff samples #make sure to add the combined runs if applicable
  bedfiles_Merged <-Reduce(function(x, y) merge(x, y, by=c("Chromosome", "Position"), all=TRUE), c(mylist) )

  #get Strand and CX columns to coalesce.
  Strands <- bedfiles_Merged %>% select(starts_with("Strand_")) %>% colnames()
  bedfiles_Merged$Strand <- do.call(dplyr::coalesce, bedfiles_Merged[Strands])

  Cxs <- bedfiles_Merged %>% select(starts_with("CX_")) %>% colnames()
  bedfiles_Merged$CX <- do.call(dplyr::coalesce, bedfiles_Merged[Cxs])
  
  #Remove unwanted columns
  bedfiles_Merged <- bedfiles_Merged  %>% select(-(starts_with("Strand_")), -starts_with(("CX_")))
  #Rearrange
  Megaframe <- bedfiles_Merged[,c(1:2,(ncol(bedfiles_Merged)-1),(ncol(bedfiles_Merged)),3:(ncol(bedfiles_Merged)-2) )]
  
  #sanity check - to ensure no NAs in Strand and CX columns after coalesce
  if ( (sum(is.na(Megaframe$Strand))==0 ) & 
       sum(is.na(Megaframe$CX))==0 ) {
    cat('QC : Megaframe looks good, Proceed to Zoomframe \n')
  } else {
    cat('QC: Strand and CX should not have Nas, re-run the megaframe function \n')
  }
  
  write.table(Megaframe, paste(project_info, "MegaFrame.csv",sep="_"), row.names=F, sep=",")

  cat("Megaframe is now available in current directory and in the Renv!")
  
  colnames(Exp_Id) = c("ID","Library")
  
  
  #Check the number of rows/sample with missing data to exclude
  rowSums(is.na(Megaframe))->Megaframe$NAs
  #make a histogram
  QCplot <- suppressMessages(ggplot(Megaframe, aes(x=NAs/3))+geom_histogram(bins=30) + 
    labs(title = "Missing data per sample") + 
    xlab("Sample") + 
    ylab("Number of rows with missing data"))
  
  print(QCplot)
  
  Megaframe_list <- list(Megaframe,Exp_Id)
  
  return (Megaframe_list)
}


#' Zoom - Adding in codes based on gene co-ordinates
#'
#' A function to add Zoom codes based on the gene positions.
#' Currently the codes for this is as below
#' 1- Region that is only between gene start and gene stop
#' 2- Region that is between Adaptive start and Adaptive Stop excluding the gene body
#' This helps subsetting easier for downstream analysis
#'
#' @param target ONT methyl bed files for each individual contained within the directory
#' @param j This is required to assign proper Alphabet codes. If you need to include the samples from a previous round, then mention the total number of samples here. Default is 0. In this case Alphabetizing starts with 'A'.
#' @param gcoord_exist This is to add a prefix to all the files that get exported while running the function.
#' @param Gene_col Use this column to specificity whether to add Gene Names or Ids in the Zoom frame.
#' @return Zoomframe Similar to Megaframe except this includes more information on targets, positions Zero'd to ATG for each target and a few other information with an additional column that included zoom codes
#' @import tidyverse
#' @import stringr
#' @export

Zoom <- function(target=CPL3, j=2, gcoord_exist=TRUE, Gene_col="Gene.Name") {


  if(gcoord_exist==TRUE){
    for (i in 1:nrow(target)) {
      if (Geneco[[Gene_col]][j]==target$Gene[i]){
        if((target$Position[i]>=Geneco$Low[j]) & (target$Position[i]<=Geneco$High[j]) ) {
          #print("entering 1")
          target$Zoom_co[i] <- 1
        } else {
          #print("entering 0")
          target$Zoom_co[i] <- 2
        }
      }   else {
        print ("Chromosomes/Column names don't match, check your Gene Coordinates file")
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
#' @param Geneco ONT methyl bed files for each individual contained within the directory
#' @param MFrame This is required to assign proper Alphabet codes. If you need to include te samples from a previous round, then mention the total number of samples here. Default is 0. In this case Alphabetizing starts with 'A'.
#' @param project_info This is to add a prefix to all the files that get exported while running the function.
#' @param filter_NAs Select this parameter based on the histogram. Ideally It should be within 1-5
#' @param gene_list Provide a list of genes that were targeted
#' @inheritParams Zoom
#' @return Zoomframe Similar to Megaframe except this includes more information on targets, positions zero'ed to ATG for each target and a few other information.
#' @import tidyverse
#' @import stringr
#' @export

getZoomF <- function(Geneco = Geneco, MFrame = MegaFrame, Gene_col="Gene.Name", filter_NAs=1, project_info="") {

  #set the filter based on how stringent it needs to be based on the plot
  MFrame[MFrame$NAs<(filter_NAs*3),]->MFrame

  cat("Creating the ZoomFrame, sit tight!\n")

  #create an empty df()
  Final_gene_set <- data.frame()
  for (i in 1:nrow(Geneco)){
    if( any(Geneco$Chromosome[i]==MFrame$Chromosome) ){ #make sure the chromosomes match between Mframe and geneco
      if ( max(MFrame$Position)>=Geneco$Adapt_Low[i] && (min(MFrame$Position)<=Geneco$Adapt_High[i]) ) {
        Gene_subset <- MFrame[MFrame$Chromosome %in% Geneco$Chromosome[i], ] #subset based on the gene
        Gene_subset$Gene <- Geneco[[Gene_col]][i] #Add-in the gene/geneID names
        if (Geneco$Strand[i]=="+") {
          Gene_subset$Zeroth_pos <- (Gene_subset$Position - Geneco$Low[i]) #Computing the Zeroth position to center everything around ATG.
        }
        else if (Geneco$Strand[i] == "-") {
          Gene_subset$Zeroth_pos <-  (Geneco$High[i]- Gene_subset$Position) #reorienting the anti-sense genes
        }
        #Call the Zoom function to add Zoom_co-ordinates
        Target_df <- Zoom(Gene_subset,i, TRUE)
        Final_gene_set <- rbind(Final_gene_set,Target_df) #append it to a Final dataframe
        #print(i)
      }
    }
    else {
      print(c(i,"Genes don't match"))
    }
  }
  Final_gene_set <- Final_gene_set[,c(1,(ncol(Final_gene_set)-2),2:4,5:(ncol(Final_gene_set)-4),ncol(Final_gene_set)-1,ncol(Final_gene_set))]

  cat("Zoomframe generated, Adding in target info column and changing NAs; Almost done!\n")
  
  #Clean up Meth and Unmeth columns with NAs
  Meth_Unmeth <- Final_gene_set %>% select(starts_with("Meth"), starts_with("UnMeth")) %>% colnames()
  cat("Columns to change NAs -> 0s\n" , Meth_Unmeth)
  #convert NAs to 0s. Here we are not changing the Percent methylation column
  Final_gene_set[Meth_Unmeth][is.na(Final_gene_set[Meth_Unmeth])] <- 0

  
  write.table(Final_gene_set, paste(project_info, "ZoomFrame.csv",sep="_"), row.names=F, sep=",")
  cat("\nZoomframe is available in your current directory!")

  return(Final_gene_set)
}
