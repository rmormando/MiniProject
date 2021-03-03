# import all of the packages
import os
#import system
#import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq


'''
Problem 1

First: retrieve the transcriptomes of given SRR numbers and convert them to paired-end fastq files.

You can use wget (by constructing the path based on the SRR numbers for each of these samples).

'''
# initialize SRR list with all of the accession numbers listed
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

# completed this
'''
def inputFiles(SRR):
    # fetches the links through SRA and downloads the files
    # wget = uses the specified path (link) to download the files onto working directory
    getFiles = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    
    # split the SRA files into paired reads
    # fastq-dump = uncompresses the data
    splitFiles = 'fastq-dump -I --split-files '+ SRR + '.1'
    
    # os.system() = executes the command in a subshell
    os.system(getFiles)
    os.system(splitFiles)
'''


'''
Problem 2

First: Build a transcriptome index for HCMV (NCBI accession EF999921)

Second: Extract the CDS features from the GenBank format.

Third: Write the output to your log file
                                                                                                                                                                                                                                                                                                                                                    
'''    

# initialize the log file (will be updated along the way)
# 'a' = to append to the file and not delete everything along the way
log_file = open('MiniProject.log', 'a')

# completed this
'''
# extract the CDS and make a transcriptome index    
def Transcriptome_Index():
    
    out_Fasta = open("EF999921.fasta", "w") # FASTA outfile for HCMV sequence
    out_CDS = open("EF999921_CDS.fasta","w") # FASTA outfile for CDS sequence
    
    Entrez.email = 'rmormando@luc.edu' # use school email to access Entrez
    
    # get nucleotide sequence in FASTA format from NCBI's [Nucleotide] database
    HCMV_file = Entrez.efetch(db = 'nucleotide', id= ['EF999921'], rettype= 'fasta') 
    
    # set the gathered info as a list called records
    records = list(SeqIO.parse(HCMV_file, "fasta")) 
    
    # use .seq object to write out FASTA file of whole sequence to the file
    # include > and separate by new line for FASTA format
    out_Fasta.write('>' + str(records[0].description)+ '\n' + str(records[0].seq))
    out_Fasta.close() # all the data has been added so close the file
    
    # get GenBank sequence from NCBI's [GenBank] database
    GB_file = Entrez.efetch(db = 'nucleotide', id = ['EF999921'], rettype= 'gb', retmode='text')
    
    # initialize tracker variable to count the number of CDS features genome has
    track = 0
    
    # this writes out the CDS sequences with a >identifier to a new file 
    # also adds a count for recording number of CDS sequences
    for record in SeqIO.parse(GB_file, 'genbank'):
        
            for feature in record.features:
                
                # find the CDS regions and then add them to the transcriptome index file = outFile
                if feature.type == "CDS": 
                    
                    track += 1 # add to the counter
                    
                    out_CDS.write('>' + str(feature.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'","") + '\n' + str(feature.location.extract(record).seq) +'\n')                
                    
    out_CDS.close() # done writing to the file so close
    
    log_file.write('The HCMV genome (EF999921) has ' + str(track) + ' CDS.' + '\n') # output how many CDS regions there are in the EF999921 genbank file

    log_file.close() # done writing to the file so close
'''


# run kallisto with the SRR numbers and files
def kallisto(SRR):
    
    # command to get kallisto
    kallisto_index = 'time kallisto index -i HCMV_index.idx EF999921_CDS.fasta'
    os.system(kallisto_index)
    
    # command to run kallisto
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o ./' + SRR + ' -b 30 -t 4 ' + SRR + '_1.fastq ' + SRR + '_2.fastq'
    os.system(run_kallisto)

# create file for input into sleuth
def input_Sleuth(SRR):
    covFile = open('covFile.txt', 'w') # make the input file for sleuth
    timePoint1 = '2dpi' # initialize time point 1
    timePoint2 = '6dpi' # initialize time point 2
    
    # first line in file (labels)
    covFile.write('Sample' + '\t' + 'Condition' + '\t' + 'Path' + '\n')
    
    # based on SRR number, write condiditon and path to output file
    for i in SRR:
        path = '/data/rmormando/MiniProject_RitaMormando/' + i
        
        if int(i[3:]) % 2 == 0: # if its even its at the first time point (2dpi)
            # write to the file using the format    
            covFile.write(str(i) + '\t' + timePoint1 + '\t' + str(path) + '\n')
        
        else: # if its odd then its at the second time point
            # write to the file using the format
            covFile.write(str(i) + '\t' + timePoint2 + '\t' + str(path) + '\n')
            
    covFile.close() # done writing to the file, close

# run sleuth in R
# read the output and add them to the log file
def Sleuth():
    runSleuth = 'Rscript sleuth.R'
    os.system(runSleuth)
    output = open('sleuth_out.txt','r')
    listoflines= output.readlines()
    
   # write to the log file the top ten sleuth hits
    with open('MiniProject.log' ,'a') as output:
        for i in range(len(file)):
                output.write(str(file[i]) + '\n') #print out the top sleuth hits
    output.close()


# ---- ALL FUNCTION CALLS ---- #
    
for i in SRR:
 #   inputFiles(i)
    kallisto(i)
    input_Sleuth(i)
#Transcriptome_Index()
Sleuth()
