library(sounDMR2)
library(tidyverse)
library(reshape2)

#-------------------------Set working directory ---------------------------------------------------------------------------------------
setwd("/Volumes/ONT/Working_directory/Pepsi/Atlantic_beds/") 
#OR
#use the Session in the menu bar and click on change working directory

#---------------------------------------------------MegaFrame section--------------------------------------------------------------------------------------------

#Read in gene target file 
#IT IS SUPER IMPORTANT TO HAVE IT IN A STANDARD FORMAT TO ENSURE THAT THE CODE DOESN'T BREAK. THE CURRENT FORMAT IS GIVEN BELOW. 
#Chromosome | Low	| High | Adapt_Start	| Adapt_End	| Gene Name	| Strand	| Plots	| Flanking | Gene_length | etc. 
Geneco <- read.csv("GeneCo.csv", header=TRUE)


#get a list of all the bed files in your path 
#It is preferred to put all the bed files that you want to work with in a single directory and change the path to that directory as mentioned. 
All_beds <- list.files(path=".",pattern="*re_aggregated.subset.bed")

#LETTERS and letters are inbuilt classes that generate a list of Alphabets for upper and lower case respectively.
#Subset to have as many number of Letters as there are individuals 
#NOTE: If the number of samples exceeds 26 and it is less than 52, make sure to subset just as much. 
#IDcat<-c(LETTERS[1:(length(All_beds))]) #change All_beds to remaining beds when necessary



#call the function and create a Megaframe
Megaframe <- GenerateMegF(Bedfiles = All_beds, N_prev_sample = 16, project_info = "Pepsi_R2_redo" )

#sanity check - to ensure no NAs in Strand and CX columns after coalesce
if ( (sum(is.na(Megaframe$Strand))==0 ) & 
      sum(is.na(Megaframe$CX))==0 ) {
  print('Megaframe looks good, Proceed to Zoomframe')
} else {
  print('Strand and CX should not have Nas, re-run the megaframe function')
}


#------------------------------------------------ZoomFrame section--------------------------------------------------------------------------------------------

#QCs
#---Check for NAs-----------------------------------
rowSums(is.na(Megaframe))->Megaframe$NAs
#make a histogram
ggplot(Megaframe, aes(x=NAs/3))+geom_histogram() + 
  labs(title = "Missing data per sample") + 
  xlab("Sample") + 
  ylab("Number of rows with missing data")

#---------------------------------------------------

targets <- c("PPO1")
#If all of the genes in the Genco are targets then use targets <- Genco$Gene . 

Zoomframe <- getZoomF(Geneco = Geneco, MFrame = Megaframe, Gene_col="Gene.Name", filter_NAs = 1, gene_list = targets, project_info = "Pepsi_R2_redo")

#It is always important to check the data. 
head(Zoomframe)
View(Zoomframe)

#QC checks
#table with 
with(Zoomframe,table(Gene,Zoom_co))

#Adaptive region unless there are overlapping gene coordinates
#IF3 <- IF3 %>% filter(Zoom_co!=0)
#double check taking an arbitrary Zeroth_pos value with the Position to see if it lines up. 
#For example if Zeroth_pos=1 then the Position should be Low(from Geneco)+1
Zoomframe %>% filter(Zeroth_pos==1)








