library(sounDMR2)


#Read in gene target file 
#IT IS SUPER IMPORTANT TO HAVE IT IN A STANDARD FORMAT TO ENSURE THE CODE DOESN'T BREAK. THE CURRENT FORMAT IS GIVEN BELOW. 
#Chromosome | Low	| High | Gene.Name	| Strand	| Gene_length | Adapt_Low | Adapt_High
Geneco <- read.table("Gene_coordinates.tsv", header=TRUE, sep="\t")


#Get a list of all the bed files in your path 
#The below command works only if all the bed files that you want to work with are in the working directory.
#Make sure to change the pattern based on your methyl bed file
All_beds <- list.files(path=".",pattern="*_methyl.bed")


#Create Megaframe
Megalist <- GenerateMegF(Bedfiles = All_beds, call_type="DSP", project_info = "Test" )

MegaFrame <- Megalist[[1]]
Exp_ID <- Megalist[[2]]

#Create Zoomframe
Zoomframe <- getZoomF(Geneco = Geneco, MFrame = MegaFrame, Gene_col="Gene.Name", filter_NAs = 1, project_info = "Test")







