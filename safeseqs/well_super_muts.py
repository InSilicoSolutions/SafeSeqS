import argparse
import os
import sys
import traceback
import logging
from safeseqs import utilities

#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-u','--uidStats', help='The UID stats file name', required=True)
    parser.add_argument('-s', '--supMutTabs', help='The Super Mutant Tabulation file name', required=True)
    parser.add_argument('-wa', '--wellAmpTabs', help='The Well Amplicon Tabulation file name', required=True)
    parser.add_argument('-b', '--barcodeMap', help='The barcode map file name', required=True)
    parser.add_argument('-o', '--output', help='The Well SuperMutant Tabulation file name', required=True)
    parser.add_argument('-sh', '--sm_homogeneity', help='SuperMutant Percent Homogeneity parameter.', required=True)
    parser.add_argument('-ir', '--indel_rate', help='Default Background Indel Rate (Percent - Used for all Reports)', required=True)
    parser.add_argument('-sr', '--sbs_rate', help='Default Background SBS Rate (Percent - Used for all Reports)', required=True)

    args = parser.parse_args()
    return args    


def perform_well_sm_tabs(args):
    ampTabs = {}
    wellSupMutTabs = {}

    args.sm_homogeneity = float(args.sm_homogeneity)/100
      
    logfile = os.path.join(utilities.get_log_dir(args.output),"WellSuperMut"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('WELL SUPER MUTANT TAB PROCESSING STARTED for %s', args.uidStats)

    try :
        #create a collection of primers with usable UIDs; write wellAmpTab file
        ampTabs, barcode = load_primers(args)

        #create a collection of changes that have the proper homogeneity
        wellSupMutTabs = load_chgs_by_primer(args, ampTabs)
        
        #write the dictionary to output file
        create_well_supMut_tab_file(args, barcode, wellSupMutTabs, ampTabs)

        logging.info('WELL SUPER MUTANT TAB PROCESSING COMPLETED for %s', args.uidStats)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)

def load_primers(args):
    ampTabs = {}
    barcode = ''
    #load a dictionary of primers found in usable UIDs from the UIDstats table
    uidStats_fh = open(args.uidStats,'r')
    for line in uidStats_fh:
        u = utilities.UidStatsRecord(*line.strip().split('\t'))
        barcode = u.barcode
        if u.usable == '1': # if the UID is usable
            # Keep a count of usable UIDs with this primer and a sum of their Family Good read Counts
            if u.primer not in ampTabs:
                ampTabs[u.primer] = [1, int(u.family_good_cnt)]
            else:
                ampTabs[u.primer][0] += 1
                ampTabs[u.primer][1] += int(u.family_good_cnt)

    uidStats_fh.close()        
    logging.info('Total primers on usable UIDs %s', str(len(ampTabs)))
    
    wellAmpTabs_fh = open(args.wellAmpTabs,'w')
    for primer in ampTabs:
        wellAmpTabs_fh.write('\t'.join([barcode, primer, str(ampTabs[primer][0]), str(ampTabs[primer][1])]) + '\n')
    wellAmpTabs_fh.close()

    return(ampTabs,barcode)


def load_chgs_by_primer(args, ampTabs):
    wellSupMutTabs = {}
    #aggregate SuperMutant Tabulation scores by primer
    supMutTabs_fh = open(args.supMutTabs,'r')
    for line in supMutTabs_fh:
        t = utilities.SuperMutantTabsRecord(*line.strip().split('\t'))
        if float(t.mut_cnt)/float(t.fam_good_reads) <= args.sm_homogeneity:
            continue
            #skip record
        
        #all superMutant Tabulations are for usable UIDs
        #key is tuple of change details: result is a list of stats
        key = (t.primer, t.chrom, t.position, t.mutType, t.baseFrom, t.baseTo)
        if key not in wellSupMutTabs:
            wellSupMutTabs[key] = [1] + [int(t.mut_cnt), int(t.mut_cnt), int(t.mut_cnt), 1] + \
                                [int(t.fam_good_reads), int(t.fam_good_reads), int(t.fam_good_reads), 1] + \
                                [int(t.minMis), float(t.avgMis), int(t.maxMis), 1] + [int(t.minCM), float(t.avgCM), int(t.maxCM), 1] + \
                                [int(t.minInd), float(t.avgInd), int(t.maxInd), 1] + [int(t.minInsB), float(t.avgInsB), int(t.maxInsB), 1] + \
                                [int(t.minDelB), float(t.avgDelB), int(t.maxDelB), 1] + [int(t.minQual), float(t.avgQual), int(t.maxQual), 1]
        else:
            #aggregate this change instance with all others of its type for this primer
            #increment the UID count (UIDs with this change)
            wellSupMutTabs[key][0] += 1
            #aggregate the mutant counts
            wellSupMutTabs[key][1:5] = aggregate_one(wellSupMutTabs[key][1:5], t.mut_cnt)
            #aggregate the family_good_read counts
            wellSupMutTabs[key][5:9] = aggregate_one(wellSupMutTabs[key][5:9], t.fam_good_reads)
            #mismatch_cnt min, avg, max, cnt 
            wellSupMutTabs[key][9:13] = aggregate_lists(wellSupMutTabs[key][9:13], t[10:13])
            #corr_mismatch_cnt min, avg, max, cnt
            wellSupMutTabs[key][13:17] = aggregate_lists(wellSupMutTabs[key][13:17], t[13:16])
            #indel_cnt min, avg, max, cnt
            wellSupMutTabs[key][17:21] = aggregate_lists(wellSupMutTabs[key][17:21], t[16:19])
            #ins_bases min, avg, max, cnt
            wellSupMutTabs[key][21:25] = aggregate_lists(wellSupMutTabs[key][21:25], t[19:22])
            #del_bases min, avg, max, cnt
            wellSupMutTabs[key][25:29] = aggregate_lists(wellSupMutTabs[key][25:29], t[22:25])
            #quality scores min, avg, max, cnt
            wellSupMutTabs[key][29:33] = aggregate_lists(wellSupMutTabs[key][29:33], t[25:28])

    supMutTabs_fh.close()        
    logging.info('Total changes on primers on usable UIDs %s', str(len(wellSupMutTabs)))

    return(wellSupMutTabs)


def aggregate_one(scores, number):
    #convert new number to an integer
    number = int(number)
    #scores is a list of Min, Sum, Max, Cnt for a set
    if number < scores[0]: #new minimum
        scores[0] = number
    if number > scores[2]: #new maximum
        scores[2] = number
    scores[1] += number #sum
    scores[3] += 1 #count

    return (scores)


def aggregate_lists(existing, new):
    #existing is a list of Min, Avg, Max, Cnt for a set 
    #new is a tuple of Min, Avg, Max from a new set to be aggregated, these are strings
    #first cast it to a list and convert to numbers
    new = list(new)
    new[0] = int(new[0])
    new[1] = float(new[1])
    new[2] = int(new[2])
    
    if new[0] < existing[0]: #minimum
        existing[0] = new[0]
    if new[2] > existing[2]: #maximum
        existing[2] = new[2]
    existing[1] += new[1] #sum of the averages
    existing[3] += 1 #count

    return (existing)


def create_well_supMut_tab_file(args, barcode, wellSupMutTabs, ampTabs):
    #assemble collected tabulation info into format for output file
    output_fh = open(args.output,'w')
    #If there are potential super mutants to write, load the reference data from files
    if len(wellSupMutTabs) > 0:
        barcodeNum, mapBarcode, template, purpose, GEs = utilities.get_barcode_details(args.barcodeMap, barcode)
        
    for change in wellSupMutTabs:
        change_details = [change[0], barcode, change[1], change[2], change[3], change[4], change[5]]
        #mut_cnts min, max, sum, cnt
        mut_cnts = finalize_scores(wellSupMutTabs[change][1:5])
        #family_good_read_cnts min, max, sum, cnt
        fam_good_cnts = finalize_scores(wellSupMutTabs[change][5:9])
        #mismatch_cnt min, max, sum, cnt
        mismatch_cnts = finalize_scores(wellSupMutTabs[change][9:13])
        #corr_mismatch_cnt min, max, sum, cnt
        corr_mismatch_cnts = finalize_scores(wellSupMutTabs[change][13:17])
        #indel_cnt min, max, sum, cnt
        indel_cnts = finalize_scores(wellSupMutTabs[change][17:21])
        #ins_bases min, max, sum, cnt
        ins_bases = finalize_scores(wellSupMutTabs[change][21:25])
        #del_bases min, max, sum, cnt
        del_bases = finalize_scores(wellSupMutTabs[change][25:29])
        #quality scores min, max, sum, cnt
        quality_scores = finalize_scores(wellSupMutTabs[change][29:33])

        #use primer from the change details to access the total uid count and sum of good reads from amplicon tabulations
        total_uids, sum_good_reads = ampTabs[change[0]]
        
        if sum_good_reads != 0:
            percent_mut_by_reads = str(wellSupMutTabs[change][2]/sum_good_reads * 100)
        else:
            percent_mut_by_reads = '0.0'

        if total_uids != 0:
            percent_mut_by_UIDs = str(wellSupMutTabs[change][0]/total_uids * 100)
        else:
            percent_mut_by_UIDs = '0.0'
        
        #use mut_type from the change details to access the tital uid cnt and sum of good reads from amplicon tabulations
        if change[3] == "SBS":
            background_rate = str(args.sbs_rate)
        else: # Indel
            background_rate = str(args.indel_rate)


        output_fh.write('\t'.join([barcodeNum, mapBarcode, template, purpose, GEs] + change_details  + \
                                  [str(wellSupMutTabs[change][0]), str(wellSupMutTabs[change][2]), str(wellSupMutTabs[change][6])] +\
                                  mut_cnts + fam_good_cnts + mismatch_cnts + corr_mismatch_cnts + \
                                  indel_cnts + ins_bases + del_bases + quality_scores +\
                                  [str(sum_good_reads), str(total_uids), percent_mut_by_reads, percent_mut_by_UIDs, background_rate]) + '\n')
    output_fh.close()
    
    
def finalize_scores(scores):
    #scores is a list of Min, sumOfAvgs, Max, Cnt
    #convert scores into strings for output, calculate average from sum and count
    str_scores = []
    str_scores.append(str(scores[0])) #minimum
    str_scores.append(str(float(scores[1])/int(scores[3]) )) #average
    str_scores.append(str(scores[2])) #maximum

    return (str_scores)


def main():
     
    args = get_args()
    
    try:
        perform_well_sm_tabs(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()