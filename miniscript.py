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

path = os.getcwd()

# completed this
'''
def inputFiles(SRR):
    # fetches the links through SRA and downloads the files
    # wget = uses the specified path (link) to download the files onto working directory
    getFiles = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
 
    # split the SRA files into paired reads
    splitFiles = 'fastq-dump -I --split-files ' + SRR + '.1'
    
    # os.system() = executes the command in a subshell
    os.system(getFiles)
    os.system(splitFiles)
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
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_'  + ' -b 30 -t 4 ' + SRR + '.1_1.fastq ' + SRR + '.1_2.fastq'
    os.system(run_kallisto)
'''

# run Sleuth using an R script
def Sleuth(SRR):
    # make the input file for sleuth
    run_R = 'Rscript sample_covariates.R ' + SRR[0] + ' ' + SRR[1]+ ' ' + SRR[2] + ' ' + SRR[3]
    os.system(run_R)

    # grab and run through R script
    run_R_again = 'Rscript Sleuth_Rscript.R' 
    os.system(run_R_again)                                                                                                                                                                      
    file = open('sleuth_out.txt').readlines()

    # write out the top sleuth hits to the log file
    with open('MiniProject.log' ,'a') as output:
        for i in range(len(file)):
                output.write(str(file[i]) + '\n') 
    output.close()
'''


# still wont work
# try new sleuth function

def SleuthInput(SRR):
    #input file for sleuth
    covFile = open('sample_covariates.txt','w+')
    #initial line in file
    covFile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to outnput file
    for i in SRR:
        path = i
        if int(i[3:]) % 2 == 0: # this line keeps messing me up
            covFile.write(str(i)+ '\t' + '2dpi' + '\t'+ path + '\n')
        else:
            covFile.write(str(i)+ '\t' + '6dpi' + '\t'+ path + '\n')
    covFile.close()


def Sleuth():
    runSleuth = 'Rscript Sleuth_Rscript.R'
    os.system(runSleuth)
    output = open('sleuth_output.txt','r').readlines()
    #listoflines= output.readlines()
    #for line in listoflines:
     #   logging.info(line)
    for i in output:
        log_file.write(i + 'n')


# wont work

# builds a bow tie index for HCMV
def build_Bowtie(SRR):
    #builds initial index using reference genome-- use transcriptome_index CDS fasta file?
    bowtie_command = 'bowtie2-build ./EF999921.fasta EF999921_1'
    os.system(bowtie_command)
    
    #maps transcriptome reads to the index we just created, generating sam file
    #need al-conc to make fastq files in output 
    bowtie_command2 = 'bowtie2 --quiet --no-unal --al-conc EF999921_' + SRR + '.fastq -x EF999921_1 -1 '+ SRR+ '.1_1.fastq -2 ' + SRR+ '.1_2.fastq -S EF999921_' + SRR+ '.sam'
    #bowtie_cmd = 'bowtie2 --quiet --no-unal --al-conc BOW_'+SRR+'.fastq -x EF99992_1 -1 '+SRR+ '_1.fastq -2'+SRR+'_2.fastq -S '+SRR+ '.sam'
    os.system(bowtie_command2)

# save the mapped reads to the index
# count the number of reads in each transcriptome before and after the bowtie mapping
def bowtie_original(SRR, number):

    bowtie_SRR1 = open('EF999921_' + SRR + '.1.fastq').readlines() # count the number of read pairs (mapped)
    bowtie_SRR2 = open('EF999921_' + SRR + '.2.fastq').readlines() # count the number of read pairs (mapped)
    donor = ''
    if number == 1:
        donor += 'Donor 1 (2dpi)' 
    elif number == 2:
        donor += 'Donor 1 (6dpi)'
    elif number == 3:
        donor += 'Donor 3 (2dpi)'
    elif number == 4:
        donor += 'Donor 3 (6dpi)' 
    len_bowtie = ((len(bowtie_SRR1)+len(bowtie_SRR2))/8) # get the length of the bowtie pair

    original1 = open(SRR + '.1_1.fastq').readlines() # count the number of read pairs (og)
    original2 = open(SRR + '.1_2.fastq').readlines() # count the number of read pairs (og)
    original = (len(original1) + len(original2))/8 # get the length of the original pair

    # write out to the log file using given format
    with open('MiniProject.log', 'a') as output:
        output.write(donor + " had " + str(original) + ' read pairs before Bowtie2 filtering and ' + str(len_bowtie) + ' read pairs after \n')
        output.close()



# ---- ALL FUNCTION CALLS ---- #

number = 1
for i in SRR:
   # inputFiles(i)
    kallisto(i)
   # Sleuth(i)
    build_Bowtie(i)
#Transcriptome_Index()
#Sleuth()
    bowtie_original(i, number)
   # SleuthInput(i)
Sleuth()
SleuthInput(SRR)
