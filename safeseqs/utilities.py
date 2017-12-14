from collections import namedtuple
import os
import logging


#define common data records for easier reference
BarcodeMapRecord = namedtuple('BarcodeMapRecord', ['barcodeNumber', 'barcode', 'wbcPlateNumber', 'template', 'purpose', 'gEsWellOrTotalULUsed', 'mutOrTotalGEsWell', 'ampMatchName', 'row', 'col'])

PrimerRecord = namedtuple('PrimerRecord', ['ampMatchName', 'gene', 'read1', 'read2', 'ampSeq', 'target_len', 'chrom', 'readStrand', 'hg19_start', 'hg19_end'])
ChromRefRecord = namedtuple('ChromRefRecord', ['chrom', 'hg19_start', 'hg19_end'])
ReferenceRecord = namedtuple('ReferenceRecord', ['chrom', 'position', 'baseFrom', 'baseTo', 'value'])
ReadRecord = namedtuple('ReadRecord', ['read_hdr', 'read_seq', 'read_qual', 'barcode', 'bc_quality', 'uid', 'uid_qual'])
UniqueSeqRecord = namedtuple('UniqueSeqRecord', ['seqUID', 'read_seq', 'read_cnt', 'primerMatch', 'read1_match', 'read1_pos', 'read2_match', 'read2_pos'])
WellFamilyRecord = namedtuple('WellFamilyRecord', ['seqUID', 'read_seq', 'barcode', 'uid', 'read_cnt', 'primer'])
UidStatsRecord = namedtuple('UidStatsRecord', ['barcode', 'uid', 'family_cnt', 'family_diversity','family_good_cnt', 'amplicon_diversity', 'primer', 'usable', 'orig_cnt', 'orig_fgr_cnt'])
AlignRecord = namedtuple('AlignRecord', ['seqUID', 'read_seq', 'read_cnt','primer', 'test_seq', 'read1', 'read1_pos', 'read2', 'read2_pos',
                                         'indel_cnt', 'mismatch_cnt', 'ins_bases', 'del_bases', 'snp_cnt', 'cosmic_cnt', 'corr_mismatch_cnt'])
ChangeRecord = namedtuple('ChangeRecord', ['seqUID', 'chrom', 'position', 'type', 'baseFrom', 'baseTo', 'cycle', 'snp', 'cosmic'])
SuperMutantTabsRecord = namedtuple('SuperMutantTabsRecord', ['primer','barcode', 'uid', 'chrom', 'position', 'mutType', 'baseFrom', 'baseTo', 'mut_cnt', 'fam_good_reads', \
                                                             'minMis', 'avgMis', 'maxMis','minCM', 'avgCM', 'maxCM','minInd', 'avgInd', 'maxInd', \
                                                             'minInsB', 'avgInsB', 'maxInsB','minDelB', 'avgDelB', 'maxDelB','minQual', 'avgQual', 'maxQual'])
WellAmpTabRecord = namedtuple('WellAmpTabRecord', ['barcode', 'primer', 'total_UIDs', 'sumGoodReads'])


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
        
        #skip header line
        if not p.hg19_start.isdigit():
            continue
        
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
        
    #write all ref items (COSMICs or dbSNPs) for chroms being used in the run whose position is within primer start and end position. 
    for line in refDB_fh:
        l = line.strip().split('\t')
        #There must be a value in the dbSNP or Somatic Count column
        if len(l) >= 5 and l[4] != '':
            #be sure to use lowercase chrom and uppercase Base From and To
            l[0] = l[0].lower()
            l[2] = l[2].upper()
            l[3] = l[3].upper()
            #if this reference record is for a chromosome in the study, check the position
            if l[0] in chromRefs:
                #check all position pairs being used by primers in the study
                for positions in chromRefs[l[0]]:
                    if l[1] >= positions[0] and l[1] <= positions[1]:
                        refSubset_fh.write('\t'.join(l) + '\n')

    refDB_fh.close()
    refSubset_fh.close()  

    return


#load chromosome positions from primer file into dictionary for use in selection of dbSNPs and COSMIC
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

    if os.path.isfile(os.path.join(os.path.dirname(filename), os.path.pardir, ref_type + ".txt")):
        
        #extract the runname directory from the filename, look for the COSMIC subset file there
        input_fh = open(os.path.join(os.path.dirname(filename), os.path.pardir, ref_type + ".txt"),'r')
        
        for line in input_fh:
            #For each line, create a data record with the value of the line (tab separated)
            r = ReferenceRecord(*line.strip().split('\t'))
                
            #Store the reference record (dbSNP or COSMIC) in the dictionary.
            if r.chrom in references:
                references[r.chrom].append((r))
            else:
                references[r.chrom] = [r]
                         
        input_fh.close()
    else:
        references['empty'] = True
        
    return references


def get_barcode_details(filename, barcode):
    num = 'NULL'
    mapBarcode = 'NULL'
    template = 'NULL'
    purpose = 'NULL'
    GEs = 'NULL'
    
    if os.path.isfile(filename):

        input_fh = open(filename,'r')
        
        for line in input_fh:
            #For each line, create a list of the fields
            b = line.strip().split('\t')
                
            #Store the barcode info in the dictionary.
            if b[1] == barcode:
                num = b[0]
                mapBarcode = b[1]
                template = b[3]
                purpose = b[4]
                GEs = b[5]
                break
                         
        input_fh.close()
        
    return (num, mapBarcode, template, purpose, GEs)


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


def load_good_reads(filename):
    good_reads = {}
    gr_fh = open(filename,'r')
    #loop through the list of good reads for the barcode and store them in a list
    for line in gr_fh:
        l = line.strip().split('\t')
        good_reads[l[0]] = True
    gr_fh.close()    
    logging.info('Loaded %s good reads', str(len(good_reads)))


    return(good_reads)


