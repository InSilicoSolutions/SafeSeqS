import argparse
import os
import sys
import traceback
import logging
from safeseqs import utilities

#import cProfile
#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Well Families file.', required=True)
    parser.add_argument('-o', '--output', help='UIDstats file.', required=True)
    parser.add_argument('-g', '--goodReads', help='File containing compressed read seqs from good reads. To be used by Optical Dups and UIDstats processes.', required=True)
    parser.add_argument('-f', '--cFc', help='File containing corrected family good read counts by UID. Counts have been corrected for any optical duplicates detected.', required=False)
    parser.add_argument('-mg', '--min_good_reads', help='The minimum good reads a UIDfamily must contain to be usable.', required=True)
    parser.add_argument('-mf', '--min_frac_good_reads', help='The minimum percent of reads that must be good reads for a UIDfamily to be usable.', required=True)
    parser.add_argument('-un', '--mark_UIDs_with_Ns_UnUsable', help='flag indicating whether uids with Ns are usable.', required=True)

    args = parser.parse_args()
    return args    


def calculate_uid_stats(args):
    
    logfile = os.path.join(utilities.get_log_dir(args.output),"UIDStats"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('UIDSTATS PROCESSING STARTED for %s', args.input)

    args.min_good_reads = int(args.min_good_reads)
    args.min_frac_good_reads = float(args.min_frac_good_reads)/100

    try :
        filename = os.path.split(args.input)[1]
        barcode = filename.split('.')[0]

        #first create a dictionary of UIDs from Well Family Reads File with each of its seqUIDs and read count.
        #this will organize and sort them from the Well Family Reads file
        UIDs = load_UIDs(args)
        cFc = load_corrected_family_counts(args)
        good_reads = utilities.load_good_reads(args.goodReads)
         
        #now, step through the UIDs dictionary. assemble and write it's UIDstats
        output_fh = open(args.output,'w')

        for uid in UIDs:
            family_cnt = 0
            family_good_cnt = 0
            family_diversity = 0
            primers = []
            amplicon = ''
            usable = 0
            
            for seqUID in UIDs[uid]:
                family_diversity += 1
                family_cnt = family_cnt + UIDs[uid][seqUID][0]
                
                if seqUID in good_reads:
                    family_good_cnt = family_good_cnt + UIDs[uid][seqUID][0]
    
                if UIDs[uid][seqUID][1] != "No Match" and UIDs[uid][seqUID][1] not in primers:
                    primers.append(UIDs[uid][seqUID][1])
                    
            #keep the original Family Read Count and Family Good Read Count before checking for corrected counts
            orig_good_cnt = family_good_cnt
            orig_cnt = family_cnt
            #if this UID has a corrected family count, use that number on its UIDstat record
            #also adjust the total reads count for the family by the number of optical duplicates found
            if uid in cFc:
                family_good_cnt = cFc[uid][0]
                family_cnt = family_cnt - cFc[uid][1]               
                
            amplicon_diversity = len(primers)
            if amplicon_diversity == 1:
                amplicon = primers[0]
                if uid.find("N") > -1 and args.mark_UIDs_with_Ns_UnUsable: #if the UID had an N AND the flag says skip Ns, don't use this one
                    usable = 0
                elif family_good_cnt >= args.min_good_reads and (family_good_cnt/family_cnt) > args.min_frac_good_reads:
                    usable = 1
            
            output_fh.write('\t'.join([barcode, uid, str(family_cnt), str(family_diversity), str(family_good_cnt), str(amplicon_diversity), amplicon, str(usable), str(orig_cnt), str(orig_good_cnt)]) + '\n')
     
        output_fh.close()

        logging.info('UIDSTATS PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)


def load_UIDs(args):
    UIDs = {}
    input_fh = open(args.input,'r')

    #read through well families, creating a dictionary of UIDs with their seqUIds(ids for unique read sequences), read counts, and primer
    for line in input_fh:
        #For each line, create a data record with the value of the line (tab separated)
        w = utilities.WellFamilyRecord(*line.strip().split('\t'))
        if w.uid not in UIDs:
            #if this is the first we have seen of this UID, create a dictionary entry to hold all of its seqUIds
            UIDs[w.uid] = {}

        #store the read count and primer for the read sequence
        if w.seqUID not in UIDs[w.uid]:
            UIDs[w.uid][w.seqUID] = (int(w.read_cnt), w.primer)

    input_fh.close()
    logging.info('Loaded %s UIDs', str(len(UIDs)))

    return(UIDs)
        

def load_corrected_family_counts(args):
    cFc = {} 
    if args.cFc is not None:
        cFc_fh = open(args.cFc,'r')
        for line in cFc_fh:
            l = line.strip().split('\t')
            cFc[l[1]] = (int(l[2]), int(l[3]))

        cFc_fh.close()
    logging.info('Loaded %s Corrected Family Counts', str(len(cFc)))

    return(cFc)
        

def main():
     
    args = get_args()
    
    try:
        calculate_uid_stats(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()