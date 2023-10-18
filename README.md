# mRNAdynamics

## Step1: Mapping
Method: 	slamdunk 

Input: 		raw reads

Output:  	aligned reads (BAM file)

## Step 2: Count nucleotide pairings for each aligned read (‚read summary table‘)

Method: 	Script confusion_table_v2, function snp_corrected_rsummary command line call python3 confusion_table_v2 denovo_rsummary [ARGS]

Input: 		filtered aligned reads (filtered BAM file), reference sequence (FASTA), SNP file (gzipped vcf)

Output: ‚read summary table file‘: a file which stores the following information on each read:
* reference sequence name mapped to
* read length
* positions of insertions (position before insertion event, w.r.t. reference, 1-based) 
* positions of deletions (position before deletion event, w.r.t. reference, 1-based)
* table listing the number of nucleotide pairings as well as insertions and deletions within the read

## Step 3: Sort aligned reads to 3’UTR regions they overlap with & divide between forward and reverse 3’UTR regions

Method:	Script confusion_table_v2, function antisense_sense_rsummary command line call python3 confusion_table_v2 filter_rsummary [ARGS]

Input:		filtered aligned reads (filtered BAM file), read summary table file, BED file defining 3’UTR regions

Output:		two read summary table files, one for all reverse and one for all forward 3’UTR regions 

## Step 4: Count the number of labeling conversions (‚conversion counts table‘)

Method:	Script confusion_table_v2, function new_convcount_table

command line call python3 confusion_table_v2 denovo_cctable [ARGS]

Input:		mode (a string that defines what kind of input table(s) is used, in our case this is ‚summary‘), a list with the reverse and forward read summary table file

Output:		‚conversion counts table file‘: a file that stores a reads‘ total number of positions where it could potentially be labeled (all mRNA transcript ‚U‘ positions), and the reads‘ number of labeling conversions actually observed


## Step 5: Estimate the labeling efficiency (for each time point and compartment and biological sample separately!) (‚conversion efficiency table file‘)

Method: 	Script confusion_table_v2, function conversion_efficiency_table command line call python3 confusion_table_v2 denovo_convefftable [ARGS]

Input:		conversion counts table file, base error (probability of observing a labeling conversion due to a sequencing error), sequencing error (probability of loosing a true labeling conversion due to a sequencing error)

Output:		‚conversion efficiency table file‘: a file that assigns each reference sequence name in the conversion counts table (in our case, each 3’UTR identifier) a labeling efficiency estimate (and a new/total ratio estimate that we do not need for the moment)

Notes: 
It is very important to estimate labeling efficiency for each measurement separately (each time point for each compartment in each biological sample; technical replicates like sequencing lanes can be thrown together) since labeling efficiency increases over time to a different extend in nucleus and cytosol!

## Step 6: Estimate the ‚true‘ new/total ratios

Method: 	Script confusion_table_v2, function conversion_efficiency_table command line call python3 confusion_table_v2 denovo_convefftable [ARGS]

Input:		conversion counts table file, base error (probability of observing a labeling conversion due to a sequencing error), sequencing error (probability of loosing a true labeling conversion due to a sequencing error), labeling efficiency estimate

Output:		‚conversion efficiency table file‘: a file that assigns each reference sequence name in the conversion counts table (in our case, each 3’UTR identifier) a new/total ratio estimate (and a labeling efficiency which doesn’t make sense for this analysis step since we fix labeling efficiency with our final estimate from Step 4)

Notes: We use the same EM algorithm as in Step 5 which is designed to estimate labeling efficiency alongside with new/total ratios, but we fix labeling efficiency to our final estimate from Step 4 such that the algorithm estimates the new/total ratios for a fixed labeling efficiency determined in the step before.
