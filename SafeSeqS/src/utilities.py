from collections import namedtuple
import os

# Utilities  Reverse Compliment
#            Load Primers
#            Primer file Format must be tab delimited and contain: AmpMatchName, Gene, 
#            Read1PrimerSeq, Read2PrimerSeq, AmpliconSeq, TargetSeqLength, Chrom, 
#            ReadStrand, Hg19Start, Hg19End
#

#define common data records for easier reference
PrimerRecord = namedtuple('PrimerRecord', ['ampMatchName', 'gene', 'read1', 'read2', 'ampSeq', 'target_len', 'chrom', 'readStrand', 'hg19_start', 'hg19_end'])
ChromRefRecord = namedtuple('ChromRefRecord', ['chrom', 'hg19_start', 'hg19_end'])
ReferenceRecord = namedtuple('ReferenceRecord', ['chrom', 'position', 'baseFrom', 'baseTo', 'value'])
ReadRecord = namedtuple('ReadRecord', ['read_hdr', 'read_seq', 'read_qual', 'barcode', 'bc_quality', 'uid', 'uid_qual'])
UniqueSeqRecord = namedtuple('UniqueSeqRecord', ['seqUID', 'read_seq', 'read_cnt', 'primerMatch', 'read1_match', 'read1_pos', 'read2_match', 'read2_pos'])
WellFamilyRecord = namedtuple('WellFamilyRecord', ['seqUID', 'read_seq', 'barcode', 'uid', 'read_cnt', 'primer'])
UidStatsRecord = namedtuple('UidStatsRecord', ['barcode', 'uid', 'family_cnt', 'family_diversity','family_good_cnt', 'amplicon_diversity', 'primer', 'usable'])
AlignRecord = namedtuple('AlignRecord', ['seqUID', 'read_seq', 'read_cnt','primer', 'test_seq', 'read1', 'read1_pos', 'read2', 'read2_pos',
                                         'indel_cnt', 'mismatch_cnt', 'ins_bases', 'del_bases', 'snp_cnt', 'cosmic_cnt', 'corr_mismatch_cnt'])
ChangeRecord = namedtuple('ChangeRecord', ['seqUID', 'chrom', 'position', 'type', 'baseFrom', 'baseTo', 'cycle', 'snp', 'cosmic'])


def get_log_dir(output):
    #strip the output file from the output string to get the output directory
    results_directory = os.path.dirname(output)
    #add \log to make the log directory
    log_directory = os.path.join(results_directory, os.path.pardir, "log")
    #if the log directory does not exist, create it
    if not os.path.isdir(log_directory):
        os.makedirs(log_directory)
    return log_directory


#load primers from file into dictionary for alignment
def load_primers(filename):
    primer_dict = {}
        
    primer_fh = open(filename,'r')
       
    for line in primer_fh:
        #For each line, create a data record with the value of the line (tab separated)
        p = PrimerRecord(*line.strip().split('\t'))
        
        #ensure read1 sequence is UPPER case
        temp_read1 = p.read1.upper()
        
        #we will need the reverse compliment of the read2 sequence
        temp_read2 = reverse_compliment(p.read2)

        #Primer file has the Hg19Start. This is position where the full reference sequence begins in the chromosome. 
        #We will need to know where the target sequence starts. Adjust the hg19start position for Read1 or Read2 length.
        if p.readStrand == '+':
        #forward primers' target sequence starts right after Read1
            temp_seq_start = str(int(p.hg19_start) + len(p.read1))
            temp_ampSeq = p.ampSeq[len(p.read1):(len(p.read1) + int(p.target_len))].upper()
        else:
        #if this is a reverse primer we are reading backwards from the end of the full reference sequence. 
        #Adjust seq_start position with Read2 length and reverse compliment the target sequence
            temp_seq_start = str(int(p.hg19_start) + len(p.read2))
            temp_ampSeq = reverse_compliment(p.ampSeq[len(p.read1):(len(p.read1) + int(p.target_len))].upper())

        #replace the originals with the revised values
        p = p._replace(read1 = temp_read1, read2 = temp_read2, hg19_start = temp_seq_start, ampSeq = temp_ampSeq)

        #Put the record in the dictionary - In this case the primer name is in the record and is the dictionary key. 
        primer_dict[p.ampMatchName] = p
            
    primer_fh.close()
    return primer_dict


def condense_ref_data(primer_file, db_file, subset_file):
    #load a dictionary of chromosomes along with a list of start and end positions used by the primers for the study
    chromRefs = load_chrom_pos(primer_file)
    
    #write subset of reference data
    refDB_fh = open(db_file,'r')
    refSubset_fh = open(subset_file,'w')
        
    #write all ref items (COSMICs or SNPs) for chroms being used in the run whose position is within primer start and end position. 
    for line in refDB_fh:
        l = line.strip().split('\t')
        #There must be a value in the SNP or Somatic Count column
        if len(l) == 5:
            chrom = l[0].lower()
            #if this reference record is for a chromosome in the study, check the position
            if chrom in chromRefs:
                #check all position pairs being used by primers in the study
                for positions in chromRefs[chrom]:
                    if l[1] >= positions[0] and l[1] <= positions[1]:
                        #be sure to use lowercase chrom and uppercase Base From and To
                        refSubset_fh.write('\t'.join([chrom, l[1], l[2].upper(), l[3].upper(), l[4]]) + '\n')

    refDB_fh.close()
    refSubset_fh.close()  

    return


#load chromosome positions from primer file into dictionary for use in selection of SNPs and COSMIC
def load_chrom_pos(filename):
    chromRefs = {}
        
    primer_fh = open(filename,'r')
       
    for line in primer_fh:
        #For each line, create a data record with the value of the line (tab separated)
        p = PrimerRecord(*line.strip().split('\t'))
        chrom = p.chrom.lower()
        #We will need to know the chrom and the hg19start & hg19end positions that are being targeted by the primers in our study.
        if chrom in chromRefs:
            chromRefs[chrom].append((p.hg19_start, p.hg19_end))
        else:
            chromRefs[chrom] = [(p.hg19_start, p.hg19_end)]           
                        
    primer_fh.close()
    return chromRefs


def load_references(filename, ref_type):
    references = {}

    if os.path.isfile(filename):
        
        #extract the study directory from the filename, look for the COSMIC subset file there
        input_fh = open(os.path.join(os.path.dirname(filename), os.path.pardir, ref_type + ".txt"),'r')
        
        for line in input_fh:
            #For each line, create a data record with the value of the line (tab separated)
            r = ReferenceRecord(*line.strip().split('\t'))
            #Store the reference record (SNP or COSMIC) in the dictionary.
            if r.chrom in references:
                references[r.chrom].append((r))
            else:
                references[r.chrom] = [r]
                         
        input_fh.close()
    else:
        references['empty'] = True
        
    return references


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

