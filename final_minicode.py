# import all of the packages
import os
import csv

from Bio import SeqIO
from Bio import Entrez


'''
Problem 1

First: retrieve the transcriptomes of given SRR numbers and convert them to paired-end fastq files.

You can use wget (by constructing the path based on the SRR numbers for each of these samples).

'''

# initialize the log file (will be updated along the way)
# 'a' = to append to the file and not delete everything along the way
log_file = open('MiniProject.log', 'a')


# initialize SRR list with all of the accession numbers listed
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

# initialize path variable to the current working directory
path = os.getcwd()


# get the input files and download as paired-read FASTQ files
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
Problem 2

First: Build a transcriptome index for HCMV (NCBI accession EF999921)

Second: Extract the CDS features from the GenBank format.

Third: Write the output to your log file
                                                                                                                                                                                                                                                                                                                                                    
''' 

# extract the CDS and make a transcriptome index for HCMV
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
                    
                    out_CDS.write('>' + str(feature.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'","") + '\n' + str(feature.location.extract(record).seq) +'\n') # write to the CDS fasta file
                    
    out_CDS.close() # done writing to the file so close
    
    log_file.write('The HCMV genome (EF999921) has ' + str(track) + ' CDS.' + '\n') # output how many CDS regions there are in the EF999921 genbank file

'''
Problem 3

First: Quantify the TPM of each CDS in each transcriptome using kallisto 

Second: Use kallisto results as input to find differentially expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth

Third: Write out the information to the log file using a header row and tab-delimit for each item:
    
    target_id test_stat pval qval
    
'''

# run kallisto with the SRR numbers and files
def kallisto(SRR):
    
    # command to get kallisto using the CDS file created before
    kallisto_index = 'time kallisto index -i HCMV_index.idx EF999921_CDS.fasta'
    os.system(kallisto_index)
    
    # command to run kallisto
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + SRR + ' -b 30 -t 4 ' + SRR + '.1_1.fastq ' + SRR + '.1_2.fastq'
    os.system(run_kallisto)

# code to make sleuth input file
def SleuthInput(SRR):
    #input file for sleuth
    covFile = open('sleuth_input.txt','w')
    #initial line in file
    covFile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to outnput file
    for i in SRR:
        paths = path + '/' + 'results_' + i
        #print(paths)
        if int(i[3:]) % 2 == 0: # if its even, 2dpi
            covFile.write(str(i)+ '\t' + '2dpi' + '\t'+ paths + '\n')
        else: # if its odd, 6dpi
            covFile.write(str(i)+ '\t' + '6dpi' + '\t'+ paths + '\n')
    covFile.close()

# code to run throgh the sleuth script in R
def Sleuth():
    runSleuth = 'Rscript Sleuth_Rscript.R'
    os.system(runSleuth)
    output = open('out_sleuth.txt','r').readlines()
    for i in output:
        log_file.write(i + '\n')


'''
Problem 4

First: Use Bowtie2 to create an index for HCMV (NCBI accession EF999921).

Second: Save the reads that map to the HCMV index for use in assembly. 

Third: Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping.

'''

# builds a bow tie index for HCMV
def build_Bowtie(SRR):
    #builds initial index using reference genome-- use transcriptome_index CDS fasta file?
    bowtie_command = 'bowtie2-build ./EF999921.fasta EF999921_1'
    os.system(bowtie_command)
    
    #maps transcriptome reads to the index we just created, generating sam file
    #need al-conc to make fastq files in output 
    bowtie_command2 = 'bowtie2 --quiet --no-unal --al-conc EF999921_' + SRR + '.fastq -x EF999921_1 -1 ' + SRR + '.1_1.fastq -2 ' + SRR + '.1_2.fastq -S EF999921_' + SRR + '.sam'
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

'''
Problem 5

First: Use the Bowtie2 output reads to assemble all four transcriptomes together to produce 1 assembly via SPAdes. 

Second: Write the SPAdes command used to the log file.

'''

# code to run through SPAdes using the command line
def Spades(SRR1, SRR2, SRR3, SRR4):
    spades_command = 'spades -k 55,77,99,127 --only-assembler -t 2 --pe1-1 EF999921_'+ SRR1 + '.1.fastq --pe1-2 EF999921_'+ SRR1 + '.2.fastq --pe2-1 EF999921_'+ SRR2 + '.1.fastq --pe2-2 EF999921_' + SRR2 + '.2.fastq --pe3-1 EF999921_' + SRR3 + '.1.fastq --pe3-2 EF999921_' + SRR3 +'.2.fastq --pe4-1 EF999921_' + SRR4 + '.1.fastq --pe4-2 EF999921_' + SRR4 + '.2.fastq -o Spades/'    
    with open('MiniProject.log','a') as output:
        output.write(spades_command + '\n') # write out the collected information to the log file
        output.close()
    os.system(spades_command) # execute the command using the command line

'''
Problem 6

First: Write Python code to calculate the number of contigs with a length > 1000

Second: Write the # to the log file as follows (replace # with the correct integer):
There are # contigs > 1000 bp in the assembly.

''' 

# count the number of contigs that are greater than 1000 bp
def count_Contigs():
    newFile = open('ContigsNum.fasta', 'w') # open up new file
    track = 0 # initialize the tracking variable to 0
    
    inFile = SeqIO.parse('./Spades/contigs.fasta', 'fasta') # parse out SPAdes output as a FASTA file
    #if the sequence len is greater than 1000 add to count
    #Add sequences greater than 1000 to outfile
    for record in inFile:
        len_seq = len(record.seq) # get the length of each sequence
        
        if len_seq > 1000: # if the sequence lenfth is greater than 1000...
            
            track += 1 # track that seqeunce add to the counter
            
            newFile.write('> '+ str(record.id) + '\n' + str(record.seq) + '\n') # add the record ID and sequence of that element to the file
    
    newFile.close() # only one condition to write to the file, so close after the end of the loop
    
    log_file.write('There are ' + str(track) + ' contigs > 1000 bp in the assembly.' + '\n') # write out to the log file

'''
Problem 7

First: Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length) and write this # to the log file as follows (replace # with the correct integer):
    
    There are # bp in the assembly.
'''

# gets the sum of the total amount of base pairs in the assembly
def length_Contigs():
    contig_file = SeqIO.parse('ContigsNum.fasta', 'fasta') # parse through the file that contains all of the sequences greater than 1000, write as FASTA
    total_Len = [] # initialize empty list for the total length
    #add each sequence len to a list
    
    for record in contig_file: # for all in the file...
        
        seq_len = len(record.seq) # collect all of their individual lengths
        
        total_Len.append(int(seq_len)) # add those lengths, as ints, to the list
    #sum the list
    total = sum(total_Len) # get the sum of all of those lengths
   
    log_file.write('There are ' + str(total) + ' bp in the assembly.' + '\n') # write out the total length (in bp) to the log file 

'''
Problem 8

First: Write Python code to retrieve the longest contig from your SPAdes assembly. 

Second: Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily. You will need to make a local database of just sequences from the Betaherpesvirinae subfamily. 

Third: Identify the top 10 hits. For the top 10 hits, write the following to your log file: Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.

Include the following header row in the log file, followed by the top 10 hits, and tab-delimit each item:
    sacc pident length qstart qend sstart send bitscore evalue stitle
'''

# function to get the longest contig from the SPAdes assembly
# SPAdes orders the sequences in order from longest to shortest when assembled
def longest_Contig():
    long_count = next(SeqIO.parse('ContigsNum.fasta', 'fasta')) # next only pulls out the first sequence in a FASTA file
    longest = open('longest_contig.fasta', 'w') # make a new file to store the longest contig
    SeqIO.write(long_count, longest, 'fasta') # write to the file

# split into three separate functions (database call, blast call, and blast parse)
# function to make the blast local database
def blast_db():
    os.system("makeblastdb -in Beta_seqs.fasta -out Betaherpesvirinae_db -title Betaherpesvirinae_db -dbtype nucl")

# function to call blast command with tab delimiters (the 10)
def blast():
    os.system('blastn -query longest_contig.fasta -db Betaherpesvirinae_db -out myBlastResults.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"')

# parse through the csv file and pull out the top 10 hits and write them to the log file
def blast_parse():
    headers = ["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"] # top line in file
    log_file.write('sacc' + '\t' + 'pident' + '\t' + 'length' + '\t' + 'qstart' + '\t' + 'qend' + '\t' + 'sstart' + '\t' + 'send' + '\t' + 'bitscore' + '\t' + 'evalue' + '\t' + 'stitle') # write out the header row to the log file
    blast_results = open("BLAST_results.csv", "r") # open the csv file
    csv_rows = csv.DictReader(blast_results, headers, delimiter = ",") # read in each row of the csv file as dictionary
    track = 0 # tracker for number of hits
    with open("MiniProject.log", "a") as out: # add top 10 hits to the log file
        for row in csv_rows:
            if track >= 9: #count used to get only top 10 hits
                break
            else:
                # pull out values for each key 
                out1 = str(row["sacc"]) # subject accession
                out2 = str(row["pident"]) # percent identity
                out3 = str(row["length"]) # alignment length
                out4 = str(row["qstart"]) # start of alignment in query
                out5 = str(row["qend"]) # end of alignment in query
                out6 = str(row["sstart"]) # start of alignment in subject
                out7 = str(row["send"]) # end of alignment in subject
                out8 = str(row["bitscore"]) # bitscore
                out9 = str(row["evalue"]) # evalue
                out10 = str(row["stitle"]) # subject title
                out.write(out1 + "\t" + out2 + "\t" +  out3 + "\t" +  out4 + "\t" +  out5 + "\t" +  out6 + "\t" +  out7 + "\t" +  out8 + "\t" +  out9 + "\t" +  out10 + "\n") # then write out all of the data found in tab delimited format to the out file
            track += 1 # track and make sure only getting top ten hits

            
            
# ---- ALL FUNCTION CALLS ---- #

number = 1

for i in SRR: # for all of the functions that are using the given SRR values
   inputFiles(i)
   kallisto(i)
   build_Bowtie(i)
   bowtie_original(i, number)
   number += 1
   
Transcriptome_Index()
Sleuth()
SleuthInput(SRR)
Spades(SRR[0], SRR[1], SRR[2], SRR[3])
count_Contigs()
length_Contigs()
longest_Contig()
blast_db()
blast()
blast_parse()
