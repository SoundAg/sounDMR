# sounDMR 2.3.1
        - In create_dmr_obj():
              - Changed calculation of Percent column in LongPercent to weighted mean
              - Fixed bug in LongMeth total_RD calculation to use sum of read depths instead of mean

# sounDMR 2.3.0
	- Removed split_by_chromosome function from workflow due to bugs 
	- Improved notes for running whole genome analysis
	- Improved memory issues with group_DMR()
	- Updated figures produced by changepoint analysis

# sounDMR 2.2.2
	- Updated generate_methylframe() documentation

# sounDMR 2.2.1
	- Updated split_by_chromosome() to ignore scaffolds in the bed file
	- Bugfix for megaframe - was only saving unfiltered megaframe. Will now write out Megaframe after the NAs filter.

# sounDMR 2.2.0
	- find_DMR() local variable bugfix
	- Parameter fix for generate_methylframe() 
	- Minor updates to upstream functions
	- Added split_by_chunk that enables creating manageable chunks of the bedfile
	- updated split_by_chrom0some

# sounDMR 2.1.3
	- Bugfix for exporting max_read_depth parameter
	- minor enhancements to run generate_megaframe faster

# sounDMR 2.1.2

	- Minor readme updates
	- Minor updates to the functions
	- Bugfix in exporting split_by_chromosome function 

# sounDMR 2.1.1

* Package updated 
	-Includes bugfix for saving files with the prefix provided by the user
	-Includes bugfix for filtering the bedfile containing > max_read_depth paarmeter in the generate_methylframe the value for which can be provided by the user. 
	
# sounDMR 2.1.0

* Package updated to handle large bed files incase of whole genome sequencing
	-Includes fix for memory leaks
	-A built-in method to split large bed files by chromosome

# sounDMR 2.0.0

* Initial public release!
