import argparse
import os
import sys
import traceback
import logging
from safeseqs import utilities


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

#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='UniqueReadSequences file.', required=True)
    parser.add_argument('-o', '--output', help='SequenceAlignment file.', required=True)
    parser.add_argument('-c', '--changes', help='Changes file.', required=True)
    parser.add_argument('-g', '--goodReads', help='File containing seqUIds and read seqs from good reads. To be used by Optical Dups and UIDstats processes.', required=True)
    parser.add_argument('-p', '--primerset', help='The primer set file name.', required=True)
    parser.add_argument('-mm', '--max_mismatches_allowed', help='The maximum number of mismatches that a read can contain and still be a good read.', required=True)
    parser.add_argument('-mi', '--max_indels_allowed', help='The maximum number of indels that a read can contain and still be a good read.', required=True)

    args = parser.parse_args()
    return args    


def perform_align(args):
    cosmics = {}
    snps = {}
    
    args.max_mismatches_allowed = int(args.max_mismatches_allowed)
    args.max_indels_allowed = int(args.max_indels_allowed)
    
    logfile = os.path.join(utilities.get_log_dir(args.output),"Align"+ str(os.getpid()) + ".log")
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
        gr_fh = open(args.goodReads,'w')

        
        #load primers from file into dictionary
        primer_dict = utilities.load_primers(args.primerset)
        
        for line in input_fh:
            #For each line, create a data record with the value of the line (tab separated)
            u = utilities.UniqueSeqRecord(*line.strip().split('\t'))
             
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
                #if there is a mismatch, insertions or deletions (indels) and/or Single Base Substitutions (SBS) exist, find a and record the changes      
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

                        changes_fh.write('\t'.join([u.seqUID, p.chrom, str(chrom_pos), "Indel", BaseFrom, BaseTo, str(cycle), 'NULL', '0']) + '\n')

                    #look for any SBS mismatches
                    #the indel logic will have adjusted for any length mismatch    
                    for i in range(0, len(p.ampSeq)):
                        if p.ampSeq[i] != test_seq[i] and test_seq[i] != "N":
                            chrom_pos = int(p.hg19_start) + i
                            cosmic_fnd = '0'
                            snp_fnd = 'NULL'

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
                        
                            #check for COSMIC, then dbSNP
                            if not cosmics:
                                #load cosmics into dictionary
                                cosmics = utilities.load_references(args.output, 'COSMIC')

                            if p.chrom in cosmics:
                                for r in cosmics[p.chrom]: #go through each entry with COSMICS on this chrom
                                    if str(chrom_pos) == r.position and p.ampSeq[i] == r.baseFrom and test_seq[i] == r.baseTo:
                                        cosmic_fnd = r.value
                                        cosmic_cnt += 1
                                        break
                            #if no COSMIC, then check for dbSNP
                            if cosmic_fnd == '0':
                                if not snps:
                                    #load dbsnps into dictionary
                                    snps = utilities.load_references(args.output, 'dbSNP')
    
                                if p.chrom in snps:
                                    for r in snps[p.chrom]: #go through each entry with dbSNPs on this chrom
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
                    #save this ID for building UIDstats later
                    gr_fh.write('\t'.join([u.seqUID, u.read_seq]) + '\n')


        input_fh.close()
        output_fh.close()
        changes_fh.close()
        gr_fh.close()
        
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


def main():
     
    args = get_args()
    
    try:
        perform_align(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()