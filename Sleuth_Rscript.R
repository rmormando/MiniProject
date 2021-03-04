# R script to run Sleuth
# mini project

# load in packages
library(sleuth)
library(data.table)
library(dplyr)


#upload the sample covariates file with the path, sample names, and conditions needed
# creates a table with the file values
cov_file <- read.table("sample_covariates.txt", header=TRUE, stringsAsFactors=FALSE, sep = ',')

# sleuth = group of kallistos
# This method takes a list of samples with kallisto results
# returns a sleuth object with the defined normalization of the data across samples
so <- sleuth_prep(cov_file)

# differential expression analysis
# fit a model comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')

# fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')

# perform the likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


#filtering significant results
#filter most significant results (FDR/qval < 0.05)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
sig_sleuth <- sleuth_significant %>% select(target_id, test_stat, pval, qval)

#write target id, test stat, pval and qval for significant transcript
#include header, tab-delimit
write.table(sig_sleuth, file="sleuth_out.txt",quote= FALSE,row.names= FALSE)
