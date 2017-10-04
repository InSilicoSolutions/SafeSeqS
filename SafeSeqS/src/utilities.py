from collections import namedtuple

# Utilities  Reverse Compliment
#            Load Primers
#            Primer file Format must be tab delimited and contain: AmpMatchName, Gene, 
#            Read1PrimerSeq, Read2PrimerSeq, AmpliconSeq, TargetSeqLength, Chrom, 
#            ReadStrand, Hg19Start, Hg19End
#

#define the primer data record
PrimerRecord = namedtuple('PrimerRecord', ['ampMatchName', 'gene', 'read1', 'read2', 'ampSeq', 'target_len', 'chrom', 'readStrand', 'testSeq_start', 'testSeq_end'])

#load primers from file into dictionary
def load_primers(filename):
    primer_dict = {}
        
    primer_fh = open(filename,'r')
       
    for line in primer_fh:
        #For each line, create a data record with the value of the line (tab separated)
        p = PrimerRecord(*line.strip().split('\t'))
        
        #we will need the reverse compliment of the read2 sequence
        temp_read2 = reverse_compliment(p.read2)

        #Primer file has the Hg19Start. This is position where the full reference sequence begins in the chromosome. 
        #We will need to know where the target sequence starts. Adjust the HG19start position for Read1 or Read2 length.
        if p.readStrand == '+':
        #forward primers' target sequence starts right after Read1
            temp_seq_start = str(int(p.testSeq_start) + len(p.read1))
            temp_ampSeq = p.ampSeq[len(p.read1):(len(p.read1) + int(p.target_len))].upper()
        else:
        #if this is a reverse primer we are reading backwards from the end of the full reference sequence. 
        #Adjust seq_start position with Read2 length and reverse compliment the target sequence
            temp_seq_start = str(int(p.testSeq_start) + len(p.read2))
            temp_ampSeq = reverse_compliment(p.ampSeq[len(p.read1):(len(p.read1) + int(p.target_len))].upper())

        #replace the originals with the revised values
        p = p._replace(read2 = temp_read2, testSeq_start = temp_seq_start, ampSeq = temp_ampSeq)

        #Put the record in the dictionary - In this case the primer name is in the record and is the dictionary key. 
        primer_dict[p.ampMatchName] = p
            
    primer_fh.close()
    return primer_dict


#reverse the order of the string and substitute the opposite nucleotide
def reverse_compliment(string):

    revCompSeq = ""
    for c in string[::-1]:
        if c == "A" or c == "a":
            revCompSeq = revCompSeq + "T"
        elif  c == "C" or c == "c":
            revCompSeq = revCompSeq + "G"
        elif  c == "G" or c == "g":
            revCompSeq = revCompSeq + "C"
        elif  c == "T" or c == "t":
            revCompSeq = revCompSeq + "A"
        else:
            revCompSeq = revCompSeq + "N"

    return revCompSeq

