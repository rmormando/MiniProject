# Computational Biology 483 - Mini Project 

This repository contains Python code from an LUC Computational Biology project.


## Software Requirements
- Linux/Unix
- Python3
- Biopython
- Kallisto
- Bowtie2
- SPAdes


## Packages Imported for Python
- os
- csv
- BioPython (Entrez and SeqIO)


## Packages Imported for R
- sleuth
- dplyr


## Run the Pipeline:

To run this script, first copy the following line of code to clone this repository to your local machine:

```
git clone https://github.com/rmormando/MiniProject_RitaMormando.git
```

Then use this line to run through the python script:

```
python3 final_minicode.py
```


## Installation Location:
```
https://github.com/rmormando/MiniProject_RitaMormando.git
```


## Using the Test Data:

You can run through the script using the given SRR files included in the python code. However, if you want to save time and run through this program quickly use the test samples that are included on this repository.

The test samples contain the first 10000 lines of the paired-end fastq files.


## Files included in Repo:

- final_minicode.py
> the whole pipeline which is composed of multiple functions that give each outputs

- Sleuth_Rscript.R
> R script which outputs a text file which identifies the diffrence between 2 timepoints in the expressed genes

- MiniProject.log
> a log file that has all of the information from the data ran with the script


#### Output Files included in Repo:

- HCMV_index.idx
> index created to be used by kallisto

- EF999921_CDS.fasta
> FASTA file with CDS of EF999921

- myBlastResults.csv
> CSV file that has contains the results after run through BLAST

- longest_contig.fasta
> FASTA file that contains the longest contig found while running SPAdes

- out_sleuth.txt
> TXT file used by R and then written to the log file


## Assignment Questions:

1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following
transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by
constructing the path based on the SRR numbers for each of these samples). 
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360 
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375 

2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI
accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with
kallisto. You will need to extract the CDS features from the GenBank format. Write the following to your log file (replace # with the number of coding sequences in the HCMV genome):
The HCMV genome (EF999921) has # CDS.

3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially
expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth. Write the following details for
each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item:
target_id test_stat pval qval

4. It has been proposed that HCMV disease and pathogenesis may be related to the genetic diversity of the virus
(Renzette et al. https://www.ncbi.nlm.nih.gov/pubmed/25154343/). Which publicly available strains are most similar to
these patient samples? To compare to other strains, we will assemble these transcriptome reads. We don’t expect
assembly to produce the entire genome, but enough to be useful in BLAST. Virus sequencing experiments often include
host DNAs. It is difficult to isolate the RNA of just the virus (as it only transcribes during infection of the host cell). Before
assembly, let’s make sure our reads map to HCMV. Using Bowtie2, create an index for HCMV (NCBI accession EF999921).
Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each
transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would
write to the log (numbers here are arbitrary):
Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.

5. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write
the SPAdes command used to the log file.

6. Write Python code to calculate the number of contigs with a length > 1000 and write the # to the log file as follows
(replace # with the correct integer):
There are # contigs > 1000 bp in the assembly.

7. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in
length) and write this # to the log file as follows (replace # with the correct integer):
There are # bp in the assembly.

8. Write Python code to retrieve the longest contig from your SPAdes assembly. Use the longest contig as blast+ input to
query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily. You will need to make a local
database of just sequences from the Betaherpesvirinae subfamily. Identify the top 10 hits. For the top 10 hits, write the
following to your log file: Subject accession, Percent identity, Alignment length, Start of alignment in query, End of
alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
Include the following header row in the log file, followed by the top 10 hits, and tab-delimit each item:
sacc pident length qstart qend sstart send bitscore evalue stitle 


## Author
Rita Mormando

Bioinformatics Program, Loyola University Chicago
