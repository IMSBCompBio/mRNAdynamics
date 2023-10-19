# mRNAdynamics

This repository is associated with our work on modelling of mRNA metabolism. A detailed read can be found in our latest preprint on bioRxiv, called "Nuclear export is a limiting factor in eukaryotic mRNA metabolism" (doi: https://doi.org/10.1101/2023.05.04.539375). Please use [SlamDunk](https://t-neumann.github.io/slamdunk/) as a mapping tool. Raw and processed sequencing data used in this work has been deposited in NCBI's Gene Expression Omnibus and are accessible through GEO Series accession number [GSE233546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233546).

## Step1: Mapping
Method: 	SlamDunk 

Input: 		raw reads

Output:  	aligned reads (BAM file)

## Step 2: Count nucleotide pairings for each aligned read (‚read summary table‘)

Method: 	Script `pipeline.py`, function `snp_corrected_rsummary` command line call `python3 pipeline.py denovo_rsummary [ARGS]`

Input: 		filtered aligned reads (filtered BAM file), reference sequence (FASTA), SNP file (gzipped vcf)

Output: ‚read summary table file‘: a file which stores the following information on each read:
* reference sequence name mapped to
* read length
* positions of insertions (position before insertion event, w.r.t. reference, 1-based) 
* positions of deletions (position before deletion event, w.r.t. reference, 1-based)
* table listing the number of nucleotide pairings as well as insertions and deletions within the read

## Step 3: Sort aligned reads to 3’UTR regions they overlap with & divide between forward and reverse 3’UTR regions

Method:	Script `pipeline.py`, function `antisense_sense_rsummary` command line call `python3 pipeline.py filter_rsummary [ARGS]`

Input:		filtered aligned reads (filtered BAM file), read summary table file, BED file defining 3’UTR regions

Output:		two read summary table files, one for all reverse and one for all forward 3’UTR regions 

## Step 4: Count the number of labeling conversions (‚conversion counts table‘)

Method:	Script `pipeline.py`, function `new_convcount_table`

command line call `python3 pipeline.py denovo_cctable [ARGS]`

Input:		mode (a string that defines what kind of input table(s) is used, in our case this is ‚summary‘), a list with the reverse and forward read summary table file

Output:		‚conversion counts table file‘: a file that stores a reads‘ total number of positions where it could potentially be labeled (all mRNA transcript ‚U‘ positions), and the reads‘ number of labeling conversions actually observed


## Step 5: Estimate the labeling efficiency (for each time point and compartment and biological sample separately!) (‚conversion efficiency table file‘)

Method: 	Script `pipeline.py`, function `conversion_efficiency_table` command line call `python3 pipeline.py denovo_convefftable [ARGS]`

Input:		conversion counts table file, base error (probability of observing a labeling conversion due to a sequencing error), sequencing error (probability of loosing a true labeling conversion due to a sequencing error)

Output:		‚conversion efficiency table file‘: a file that assigns each reference sequence name in the conversion counts table (in our case, each 3’UTR identifier) a labeling efficiency estimate (and a new/total ratio estimate that we do not need for the moment)

Notes: 
It is very important to estimate labeling efficiency for each measurement separately (each time point for each compartment in each biological sample; technical replicates like sequencing lanes can be thrown together) since labeling efficiency increases over time to a different extend in nucleus and cytosol!

## Step 6: Estimate the ‚true‘ new/total ratios

Method: 	Script `pipeline.py`, function `conversion_efficiency_table` command line call `python3 pipeline.py denovo_convefftable [ARGS]`

Input:		conversion counts table file, base error (probability of observing a labeling conversion due to a sequencing error), sequencing error (probability of loosing a true labeling conversion due to a sequencing error), labeling efficiency estimate

Output:		‚conversion efficiency table file‘: a file that assigns each reference sequence name in the conversion counts table (in our case, each 3’UTR identifier) a new/total ratio estimate (and a labeling efficiency which doesn’t make sense for this analysis step since we fix labeling efficiency with our final estimate from Step 4)

Notes: We use the same EM algorithm as in Step 5 which is designed to estimate labeling efficiency alongside with new/total ratios, but we fix labeling efficiency to our final estimate from Step 4 such that the algorithm estimates the new/total ratios for a fixed labeling efficiency determined in the step before.

## Step 7: Create Summary Tables

Method: 	Script `pipeline.py`, function `new_summary_table` command line call `python3 pipeline.py denovo_summarytable [ARGS]`

Input:		conversion effiency table file, base error (probability of observing a labeling conversion due to a sequencing error), sequencing error (probability of loosing a true labeling conversion due to a sequencing error), labeling efficiency estimate

Output:		‚summary table file‘: a summary table for a given a conversion counts table and conversion efficiency table file
    and writes the table to an output file. 

If summary tables for all time points have been generated, they can be merged by calling `merge_summarytable [ARGS]` using the tables as input.

## Step 8: Compute expression level and prepare rawdata

Method:    Script `helper_functions.R`, functions `expr_level_regression` and `summarytable_to_rawdata`

Input:		A merged summary table file, containing the samples for one time-series. The `expr_level_regression` function is used to calculate the expression level of each gene along the time-series (needed for assignment of reliability criteria later). The `summarytable_to_rawdata` converts the summarytable to a data frame that is used for the model fit.

Output:		An RData-object which stores the rawdata that is used for the estimation of metabolic rates.

## Step 9: Estimation of RNA metabolic rates


Method:    Script `modelfit.R`, functions `fitGene`

Input:	The function fitGene takes each entry from the rawdata generated in Step 8. 

Output:    An RData-object which stores the estimation of RNA metabolic rates.

## Notes
A full list of RNA metabolic rates we obtained from the time-series data used in our preprint "Nuclear export is a limiting factor in eukaryotic mRNA metabolism" (doi: https://doi.org/10.1101/2023.05.04.539375) will also be available with publication of this manuscript.
