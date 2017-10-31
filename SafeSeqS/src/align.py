import argparse
import os
import sys, traceback
import logging
from collections import namedtuple
import utilities

#import cProfile
# Align      This program takes the output after the Unique Reads have been identified .
#            The output is a tab delimited file with one line per unique read sequence with both a Read1 match and Read2 match.
#
#            It creates a tab delimited Alignment Sequence file with: 
#                ReadSequence
#                ReadCount
#                PrimerName
#                TestSequence
#                Read1PrimerMatch
#                Read1MatchPosition
#                Read2PrimerMatch
#                Read2MacthPosition
#                IndelCount
#                MismatchCount
#                InsertedBases
#                DeletedBases
#                later more columns may be added
#                SNPcount
#                COSMICcount
#                CorrectedMismatchCount
#
#            It creates a tab delimited Changes file with: 
#                seqUID
#                Chrom
#                Position
#                MutType
#                BaseFrom
#                BaseTo
#                AlignmentCycle
#                SNP
#                COSMIC
#

#define the Unique Read Sequence input data record for easier reference
UniqueSeqRecord = namedtuple('UniqueSeqRecord', ['seqUID', 'read_seq', 'read_cnt', 'primerMatch', 'read1_match', 'read1_pos', 'read2_match', 'read2_pos'])

#define the Well Family input data record for easier reference
WellFamilyRecord = namedtuple('WellFamilyRecord', ['seqUID', 'read_seq', 'barcode', 'uid', 'read_cnt'])

#create a global dictionary to store primers form the primer file
primer_dict = {}

#create a global dictionary to store seqUIDs and primers for reads that match Good Read criteria
good_reads = {}
 
#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='UniqueReadSequences file.', required=True)
    parser.add_argument('-o', '--output', help='SequenceAlignment file.', required=True)
    parser.add_argument('-c', '--changes', help='Changes file.', required=True)
    parser.add_argument('-w','--wellfamilies', help='Well Families file.', required=True)
    parser.add_argument('-u', '--UIDstats', help='UIDstats file.', required=True)
    parser.add_argument('-p', '--primerset', help='The primer set file name.', required=True)
    parser.add_argument('-mm', '--max_mismatches_allowed', help='The maximum number of mismatches that a read can contain and still be a good read.', required=True)
    parser.add_argument('-mi', '--max_indels_allowed', help='The maximum number of indels that a read can contain and still be a good read.', required=True)
    parser.add_argument('-mg', '--min_good_reads', help='The minimum good reads a UIDfamily must contain to be usable.', required=True)
    parser.add_argument('-mf', '--min_frac_good_reads', help='The minimum percent of reads that must be good reads for a UIDfamily to be usable.', required=True)
    parser.add_argument('-un', '--use_UIDs_with_Ns', help='flag indicating whether uids with Ns are usable.', required=True)

    args = parser.parse_args()
    return args    


def get_log_dir(output):
    #strip the output file from the output string to get the output directory
    results_directory = os.path.dirname(output)
    #add \log to make the log directory
    log_directory = os.path.join(results_directory, os.path.pardir, "log")
    #if the log directory does not exist, create it
    if not os.path.isdir(log_directory):
        os.makedirs(log_directory)
    return log_directory


def perform_align(args):
    global primer_dict
    global good_reads
    cosmics = {}
    snps = {}
    
    args.max_mismatches_allowed = int(args.max_mismatches_allowed)
    args.max_indels_allowed = int(args.max_indels_allowed)
    args.min_good_reads = int(args.min_good_reads)
    args.min_frac_good_reads = float(args.min_frac_good_reads)/100

    
    logfile = os.path.join(get_log_dir(args.output),"Align"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('ALIGN PROCESSING STARTED for %s', args.input)

    try :
        input_fh = open(args.input,'r')
        output_fh = open(args.output,'w')
        changes_fh = open(args.changes,'w')
        
        #load primers from file into dictionary
        primer_dict = utilities.load_primers(args.primerset)
        
        for line in input_fh:
            #For each line, create a data record with the value of the line (tab separated)
            u = UniqueSeqRecord(*line.strip().split('\t'))
             
            test_seq = 'NULL'
            mismatch_cnt = 0
            snp_cnt = 0
            cosmic_cnt = 0
            indel_cnt = 0
            ins_bases = 0
            del_bases = 0            

            if u.read1_match == 'Perfect Match' and u.read2_match == 'Perfect Match':
                #retrieve the specific primer info from primer_dict
                p = primer_dict[u.primerMatch]
                
                #test sequence is found in the read sequence between the read1 and read2 primer sequences
                if p.readStrand == '+':
                    test_seq = u.read_seq[len(p.read1):int(u.read2_pos)-1].upper()
                else: 
                    #if this is a reverse primer, perform reverse compliment on the test sequence
                    test_seq = utilities.reverse_compliment(u.read_seq[len(p.read1):int(u.read2_pos) - 1].upper())

                #if the test sequence matches the target amplicon sequence, we have a perfect match. Simply write the alignment record with no changes.  
                if test_seq != p.ampSeq:
                #if there is a mismatch, insertions or deletions (indels) and/or Single Base Substitutions (SBS) exist, finda and record the changes      
                    pos = -1 #initialize the indel loop counter to -1 so that if there is no indel, SBS mismatch check won't think it was position 0
                    indel_len = len(test_seq) - len(p.ampSeq)
                    
                    # A difference in lengths indicates that an InDel exists. Find it and modify the test sequence to add/remove it.
                    if indel_len != 0:
                        indel_cnt += 1

                        pos, BaseFrom, BaseTo, test_seq = get_indel(p.ampSeq, test_seq)
                        if indel_len > 0: #insertion
                            ins_bases = indel_len
                            chrom_pos = int(p.hg19_start) + pos - 1

                            if p.readStrand == '+': 
                                cycle = int(u.read1_pos) + len(p.read1) + pos
                            else:
                                cycle = int(u.read2_pos) - pos - 1 #if the primer indicates that the read is reverse, we count backwards from the beginning of read2
                            
                        else:
                            del_bases = indel_len #deletion
                            chrom_pos = int(p.hg19_start) + pos

                            if p.readStrand == '+':
                                cycle = int(u.read1_pos) + len(p.read1) + pos - 1
                            else:
                                cycle = int(u.read2_pos) - pos #if the primer indicates that the read is reverse, we count backwards from the beginning of read2

                        changes_fh.write('\t'.join([u.seqUID, p.chrom, str(chrom_pos), "Indel", BaseFrom, BaseTo, str(cycle)]) + '\n')

                    #look for any SBS mismatches
                    #the indel logic will have adjusted for any length mismatch    
                    for i in range(0, len(p.ampSeq)):
                        if p.ampSeq[i] != test_seq[i] and test_seq[i] != "N":
                            chrom_pos = int(p.hg19_start) + i
                            cosmic_fnd = ''
                            snp_fnd = ''

                            #cycle is exact position of the SBS on original read sequence, we will need it to evaluate the read quality score later
                            #if there was an indel, adjust for its length
                            if p.readStrand == '+':
                                cycle = int(u.read1_pos) + len(p.read1) + i
                                if i > pos:
                                    cycle = cycle + indel_len #adjust by +/- length of indel
                            else:
                                cycle = int(u.read2_pos) - i - 1 
                                if i > pos:
                                    cycle = cycle - indel_len #adjust by +/- length of indel
                        
                            #check for COSMIC, then SNP
                            if not cosmics:
                                #load cosmics into dictionary
                                cosmics = utilities.load_references(args.output, 'COSMIC')

                            if p.chrom in cosmics:
                                for r in cosmics[p.chrom]: #go through each entry with COSMICS on this chrom
                                    if str(chrom_pos) == r.position and p.ampSeq[i] == r.baseFrom and test_seq[i] == r.baseTo:
                                        cosmic_fnd = r.value
                                        cosmic_cnt += 1
                                        break
                            #if no COSMIC, then check for SNP
                            if cosmic_fnd == '':
                                if not snps:
                                    #load snps into dictionary
                                    snps = utilities.load_references(args.output, 'SNP')
    
                                if p.chrom in snps:
                                    for r in snps[p.chrom]: #go through each entry with SNPs on this chrom
                                        if str(chrom_pos) == r.position and p.ampSeq[i] == r.baseFrom and test_seq[i] == r.baseTo:
                                            snp_fnd = r.value
                                            snp_cnt += 1
                                            break
                                    
                            changes_fh.write('\t'.join([u.seqUID, p.chrom, str(chrom_pos), 'SBS', p.ampSeq[i], test_seq[i], str(cycle), snp_fnd, cosmic_fnd]) + '\n')
                            #increment number of mismatches found in this sequence alignment record
                            mismatch_cnt +=1
                
                corrected_mismatch_cnt = mismatch_cnt - snp_cnt - cosmic_cnt
                output_fh.write('\t'.join([u.seqUID, u.read_seq, u.read_cnt, p.ampMatchName, test_seq, u.read1_match, u.read1_pos, u.read2_match, u.read2_pos, str(indel_cnt), str(mismatch_cnt), str(ins_bases), str(del_bases), str(snp_cnt), str(cosmic_cnt), str(corrected_mismatch_cnt)]) + '\n')

                #determine if this unique read sequence meets the settings for a Good Read
                if corrected_mismatch_cnt <= args.max_mismatches_allowed and indel_cnt <= args.max_indels_allowed:
                    #save this ID and it's primer for building UIDstats later
                    good_reads[u.seqUID] = p.ampMatchName
                        
        input_fh.close()
        output_fh.close()
        changes_fh.close()
        
        calculate_uid_stats(args)

        logging.info('ALIGN PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)


# Find the indel location and substitution values using the reference sequence and test sequence.
def get_indel(ref, test):
    max_score = -1
    pos = 0
    BaseFrom = 'NULL'
    BaseTo = 'NULL'
    best = ''
    
    indel_len = len(test) - len(ref)
    if indel_len < 0: #test sequence is too short, a deletion occurred
        filler = "N" * abs(indel_len)
        #step thru the test sequence, adding a filler of chars. check if resulting substring is best match to reference string
        for i in range(0, (len(test)+1)):
            if i == len(test):  #for last time thru, just add the filler to the end of test string
                temp = test + filler
            else:               #put together a temp test string by inserting the filler at loop counter position
                temp = test[0:i] + filler + test[i:len(test)]
            score = compare_sequences(ref, temp)
            if score >= max_score: #check if this temp test string outscores previous iterations
                max_score = score
                pos = i
                best = temp
                BaseFrom = ref[i:(i+abs(indel_len))]

    elif indel_len > 0: #test sequence is too long, an insertion occurred with the length of "indel_len"
        #step thru the test sequence, removing a chunk of chars. check if resulting substring is best match to reference string
        for i in range(0, (len(ref)+1)): 
            if i == len(ref): #for last time thru, just lop the number of inserted chars off the end of test string
                temp = test[0:i]
            else:                 #put together a temp test string by removing the number of inserted chars at loop counter position
                temp = test[0:i] + test[(i + indel_len):len(test)]
            score = compare_sequences(ref, temp)
            if score >= max_score: #check if this temp test string outscores previous iterations
                max_score = score
                pos = i
                best = temp
                BaseTo = test[pos:(pos+indel_len)]
                    
    return pos, BaseFrom, BaseTo, best


# compare two strings. They must be the same length to get a score. Increment score for each matching character.
def compare_sequences(reference, test):
    score = 0
    
    if len(reference) == len(test):
        if reference == test:
            score = len(reference)
        else:
            for i in range(0, len(reference)):
                if reference[i] == test[i]:
                    score += 1
                    
    return score

def calculate_uid_stats(args):
    #first create a dictionary of UIDs from Well Family File with each seqUID and count.
    #this will organize and sort them from the Well Family file
    UIDs = {}

    input_fh = open(args.wellfamilies,'r')
    #read through well families, creating a dictionary of UIDs with their seqUIds(ids for unique read sequences) and read counts
    for line in input_fh:
        #For each line, create a data record with the value of the line (tab separated)
        w = WellFamilyRecord(*line.strip().split('\t'))
        barcode = w.barcode
        if w.uid not in UIDs:
            #if this is the first we have seen of this UID, create a dictionary entry to hold all of its seqUIds
            UIDs[w.uid] = {}

        #store the read count for the read sequence
        if w.seqUID not in UIDs[w.uid]:
            UIDs[w.uid][w.seqUID] = int(w.read_cnt)

    input_fh.close()  

    #now, step through the UIDs dictionary and assemble it's UIDstats
    us_fh = open(args.UIDstats,'w')

    for uid in UIDs:
        family_cnt = 0
        family_good_cnt = 0
        family_diversity = 0
        primers = []
        amplicon = ''
        usable = 0
        
        for seqUID in UIDs[uid]:
            family_diversity += 1
            family_cnt = family_cnt + UIDs[uid][seqUID]
            
            if seqUID in good_reads:
                family_good_cnt = family_good_cnt + UIDs[uid][seqUID]
                if good_reads[seqUID] not in primers:
                    primers.append(good_reads[seqUID])

        amplicon_diversity = len(primers)
        if amplicon_diversity == 1:
            amplicon = primers[0]


        if uid.find("N") > -1 and args.use_UIDs_with_Ns == 'N': #if the UID had an N AND the flag says skip Ns, don't use this one
            continue
        elif family_good_cnt >= args.min_good_reads and (family_good_cnt/family_cnt) > args.min_frac_good_reads:
            usable = 1
        
        us_fh.write('\t'.join([barcode, uid, str(family_cnt), str(family_diversity), str(family_good_cnt), str(amplicon_diversity), amplicon, str(usable)]) + '\n')
 
    us_fh.close()
    return

    
def main():
     
    args = get_args()
    
    try:
        perform_align(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()
#if __name__ == "__main__": cProfile.run('main()')