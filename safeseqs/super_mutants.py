import argparse
import os
import sys
import traceback
import logging
import zlib
from safeseqs import utilities

#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--aligns', help='The align file name', required=True)
    parser.add_argument('-c','--changes', help='The changes file name', required=True)
    parser.add_argument('-u','--uidStats', help='The UID stats file name', required=True)
    parser.add_argument('-f','--familyReads', help='The well family file name', required=True)
    parser.add_argument('-r','--reads', help='The reads file name', required=True)
    parser.add_argument('-o', '--output', help='The Super Mutant Tabulation file name', required=True)
    parser.add_argument('-mm', '--max_mismatches_allowed', help='The maximum number of mismatches that a read can contain and still be a good read.', required=True)
    parser.add_argument('-mi', '--max_indels_allowed', help='The maximum number of indels that a read can contain and still be a good read.', required=True)
    parser.add_argument('-aa', '--ascii_adj', help='The constant that must be subtracted from the ascii quality score to get the integer score.', required=True)

    args = parser.parse_args()
    return args    


def perform_sup_mut_tab(args):
    usableReads = {}
    chgPositions = {}
    qualScores = {}
    supMutTab = {}
    
    args.max_indels_allowed = int(args.max_indels_allowed)
    args.max_mismatches_allowed = int(args.max_mismatches_allowed)
    args.ascii_adj = int(args.ascii_adj)
      
    logfile = os.path.join(utilities.get_log_dir(args.output),"SuperMut"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('SUPER MUTANT TAB PROCESSING STARTED for %s', args.uidStats)

    try :
        #create a collection of read sequences that belong to a usable UID
        usableReads, barcode = load_usable_reads(args)

        #create a collection of changes that exist for good reads
        chgPositions = load_good_chg_positions(args)
        
        #aggregate quality scores for positions with changes on good reads that belong to a usable UID
        qualScores = load_read_qual_at_chg_pos(args, chgPositions, usableReads)
        
        supMutTab = aggregate_scores(qualScores, chgPositions, usableReads)
                    
        #write the dictionary to output file
        create_sup_mut_tab_file(args, barcode, supMutTab, usableReads)

        logging.info('SUPER MUTANT TAB PROCESSING COMPLETED for %s', args.uidStats)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)

def load_usable_reads(args):
    usableReads = {}
    barcode = ''
    #load a dictionary of usable UIDs from the UIDstats table
    uidStats_fh = open(args.uidStats,'r')
    for line in uidStats_fh:
        u = utilities.UidStatsRecord(*line.strip().split('\t'))
        barcode = u.barcode
        if u.usable == '1': # if the UID is usable, we store its family good read count 
            usableReads[u.uid] = (u.family_good_cnt, {}) # empty dictionary will store unique read records for UID

    uidStats_fh.close()        
    logging.info('Total usable UIDs %s', str(len(usableReads)))
    
    #use well family reads to collect the unique read sequences for the usable UIDs. store seqUID and family count
    familyReads_fh = open(args.familyReads,'r')
    cnt = 0
    for line in familyReads_fh:
        w = utilities.WellFamilyRecord(*line.strip().split('\t'))
        if w.uid in usableReads: #if read is for a usable UID, we will need to collect its read sequence, seqUID, and family read count
            cnt +=1
            reads = usableReads[w.uid][1] #reads is a pointer to the dictionary inside the tuple for this UID

            #compress the read to use it as key for dictionary
            compRead = zlib.compress (w.read_seq.encode('utf-8'),9)
            reads[compRead] = (w.seqUID, w.read_cnt) #store the seqUID and family read count, keyed on the read_sequence
    familyReads_fh.close()
    logging.info('Reads added to usable UIDs %s', str(cnt))

    return(usableReads,barcode)

    
def load_good_chg_positions(args):
    chgPositions = {}
    #load a dictionary of seqUIDs whose align record has changes within mismatch and indel limits
    aligns_fh = open(args.aligns,'r')
    for line in aligns_fh:
        a = utilities.AlignRecord(*line.strip().split('\t'))
        #there should be at least one change (mismatch or indel) AND indel and mismatches w/o snp or cosmics must be w/i settings
        if (int(a.indel_cnt) > 0 or int(a.mismatch_cnt) > 0)  and \
            int(a.indel_cnt) <= args.max_indels_allowed and int(a.corr_mismatch_cnt) <= args.max_mismatches_allowed:
            
            stats = [a.primer, int(a.indel_cnt), int(a.mismatch_cnt), int(a.ins_bases), int(a.del_bases), int(a.corr_mismatch_cnt)]
            chgPositions[a.seqUID] = (stats, {}) # empty dictionary will store positions where changes occurred in this read sequence

    aligns_fh.close()
    logging.info('Total good aligns with changes %s', str(len(chgPositions)))
        
    #use change records to collect the position and identifying details for each change on the good aligns
    changes_fh = open(args.changes,'r')
    cnt = 0
    for line in changes_fh:
        c = utilities.ChangeRecord(*line.strip().split('\t'))
        if c.seqUID in chgPositions: #if change is for a good align, we will need to collect its positionc.position (cycle), and identifying details
            cnt +=1
            stats, changes = chgPositions[c.seqUID] #changes is a pointer to the dictionary inside the tuple for this seqUID
            #position on Chrom is unique, cycle (position on readSeq) can be duplicate if Indel and SBS occur next to each other
#            changes[int(c.position)] = [(c.chrom, c.position, c.type, c.baseFrom, c.baseTo), c.cycle] #store the change details, keyed on the position in the chrom
            #combine cycle (position on readSeq) with type of change to create a unique key 
            changes[c.cycle+c.type] = [(c.chrom, c.position, c.type, c.baseFrom, c.baseTo), c.cycle] #store the change details, keyed on the position in the chrom
    changes_fh.close()
    logging.info('Change Positions added to good aligns %s', str(cnt))
    
    return(chgPositions)


def load_read_qual_at_chg_pos(args, chgPositions, usableReads):
    qualScores ={}
    #create a dictionary by UID:seqUID:position to identify where we need to aggregate quality scores
    cnt = 0
    for UID in usableReads:
        UIDreads = usableReads[UID][1] #dictionary is in the tuple's second position
        for read in UIDreads:
            seqUID, family_cnt = UIDreads[read]
            if seqUID in chgPositions:
                changes = chgPositions[seqUID][1] #dictionary is in the tuple's second position
                #now we know there is at least one change for this read
                if UID not in qualScores:
                    qualScores[UID] = {}
                if seqUID not in qualScores[UID]:
                    qualScores[UID][seqUID] = [int(family_cnt), {}]
                for position in changes:
                    cnt += 1
                    qualScores[UID][seqUID][1][position] = [99999,0,0,0] #start with default values
    logging.info('Change Positions on Usable UIDs %s', str(cnt))

    reads_fh = open(args.reads,'r')
    for line in reads_fh:
        r = utilities.ReadRecord(*line.strip().split('\t'))
        #is this line for a usable UID with changes
        if r.uid in qualScores:
            #compress the read to compare with key for dictionary
            compRead = zlib.compress (r.read_seq.encode('utf-8'),9)
            seqUID = usableReads[r.uid][1][compRead][0]#second position at the UID level, first position at the read level
            #is this line for a read with changes
            if seqUID in qualScores[r.uid]:
                for position in qualScores[r.uid][seqUID][1]:
                    #cycle is the exact position in the read; python strings start at 0 so subtract 1
                    pos = int(chgPositions[seqUID][1][position][1]) -1
                    read_pos_qual = int(ord(r.read_qual[pos]) - args.ascii_adj)
                    qualScores[r.uid][seqUID][1][position] = aggregate_one(qualScores[r.uid][seqUID][1][position], read_pos_qual, 1)
    reads_fh.close()
    logging.info('Read Qualities aggregated on Usable UIDs Change Positions')
        
    return(qualScores)


def aggregate_scores(qualScores, chgPositions, usableReads):
    supMutTab = {}
    #aggregate quality scores, align stats by unique change details within UID
    for uid in qualScores:
        if uid not in supMutTab:
            supMutTab[uid] = {}
        for seqUID in qualScores[uid]: #each unique read with changes
            a_cnts = chgPositions[seqUID][0]
            for position in qualScores[uid][seqUID][1]:
                change = chgPositions[seqUID][1][position][0] #this will be a tuple with the change details
                if change not in supMutTab[uid]:
                    #key is tuple of change details: fam_cnt, primer, align counts -initial min, max, sum (weighted for cnt), cnt AND quality scores
                    supMutTab[uid][change] = [qualScores[uid][seqUID][0], a_cnts[0], \
                        a_cnts[1], a_cnts[1], (a_cnts[1]*qualScores[uid][seqUID][0]), qualScores[uid][seqUID][0], \
                        a_cnts[2], a_cnts[2], (a_cnts[2]*qualScores[uid][seqUID][0]), qualScores[uid][seqUID][0], \
                        a_cnts[3], a_cnts[3], (a_cnts[3]*qualScores[uid][seqUID][0]), qualScores[uid][seqUID][0], \
                        a_cnts[4], a_cnts[4], (a_cnts[4]*qualScores[uid][seqUID][0]), qualScores[uid][seqUID][0], \
                        a_cnts[5], a_cnts[5], (a_cnts[5]*qualScores[uid][seqUID][0]), qualScores[uid][seqUID][0]] + qualScores[uid][seqUID][1][position]
                else:
                    #aggregate this change instance with all others of its type in the UID
                    #add to count for how many reads were in this well_family_read count
                    supMutTab[uid][change][0] += qualScores[uid][seqUID][0]
                    #indel_cnt min, max, sum, cnt: positions 2-5
                    supMutTab[uid][change][2:6] = aggregate_one(supMutTab[uid][change][2:6], a_cnts[1], qualScores[uid][seqUID][0])
                    #mismatch_cnt min, max, sum, cnt: positions 6-9
                    supMutTab[uid][change][6:10] = aggregate_one(supMutTab[uid][change][6:10], a_cnts[2], qualScores[uid][seqUID][0])
                    #ins_bases min, max, sum, cnt: positions 10-13
                    supMutTab[uid][change][10:14] = aggregate_one(supMutTab[uid][change][10:14], a_cnts[3], qualScores[uid][seqUID][0])
                    #del_bases min, max, sum, cnt: positions 14-17
                    supMutTab[uid][change][14:18] = aggregate_one(supMutTab[uid][change][14:18], a_cnts[4], qualScores[uid][seqUID][0])
                    #corr_mismatch_cnt min, max, sum, cnt: positions 18-21
                    supMutTab[uid][change][18:22] = aggregate_one(supMutTab[uid][change][18:22], a_cnts[5], qualScores[uid][seqUID][0])
                    #quality scores min, max, sum, cnt: positions 22-25
                    supMutTab[uid][change][22:26] = aggregate_lists(supMutTab[uid][change][22:26], qualScores[uid][seqUID][1][position])
    return(supMutTab)

def aggregate_one(scores, number, cnt):
    if number < scores[0]: #minimum
        scores[0] = number
    if number > scores[1]: #maximum
        scores[1] = number
    scores[2] += (number * cnt) #sum - weighted for number of records with this value
    scores[3] += cnt #count

    return (scores)


def aggregate_lists(existing, new):
    if new[0] < existing[0]: #minimum
        existing[0] = new[0]
    if new[1] > existing[1]: #maximum
        existing[1] = new[1]
    existing[2] += new[2] #sum
    existing[3] += new[3] #count

    return (existing)


def create_sup_mut_tab_file(args, barcode, supMutTab, usableReads):
    #assemble collected tabulation info into format for output file
    output_fh = open(args.output,'w')
    for uid in supMutTab:
        for change in supMutTab[uid]:
            supMutTab[uid][change][0] = str(supMutTab[uid][change][0])
            #indel_cnt min, max, sum, cnt: positions 2-5
            indel_cnts = finalize_scores(supMutTab[uid][change][2:6])
            #mismatch_cnt min, max, sum, cnt: positions 6-9
            mismatch_cnts = finalize_scores(supMutTab[uid][change][6:10])
            #ins_bases min, max, sum, cnt: positions 10-13
            ins_bases = finalize_scores(supMutTab[uid][change][10:14])
            #del_bases min, max, sum, cnt: positions 14-17
            del_bases = finalize_scores(supMutTab[uid][change][14:18])
            #corr_mismatch_cnt min, max, sum, cnt: positions 18-21
            corr_mismatch_cnts = finalize_scores(supMutTab[uid][change][18:22])
            #quality scores min, max, sum, cnt: positions 22-25
            quality_scores = finalize_scores(supMutTab[uid][change][22:26])

            output_fh.write('\t'.join(supMutTab[uid][change][1:2] + [barcode, uid] + list(change) + supMutTab[uid][change][0:1] + [str(usableReads[uid][0])] + \
                                          mismatch_cnts + corr_mismatch_cnts + indel_cnts + ins_bases + del_bases + quality_scores) + '\n')
    output_fh.close()
    
    
def finalize_scores(scores):
    #convert scores into strings for output, calculate average from sum and count
    str_scores = []
    str_scores.append(str(scores[0])) #minimum
    str_scores.append(str(int(scores[2]/scores[3]))) #average
    str_scores.append(str(scores[1])) #maximum

    return (str_scores)


def main():
     
    args = get_args()
    
    try:
        perform_sup_mut_tab(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()