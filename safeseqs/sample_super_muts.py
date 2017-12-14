import argparse
import os
import sys
import traceback
import logging
from safeseqs import utilities

#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--barcodeMap', help='The barcode map file name', required=True)
    parser.add_argument('-t', '--template', help='The template (sample) to pull from barcode map file', required=True)
    parser.add_argument('-o', '--output', help='The Sample SuperMutant Tabulation file name', required=True)

    args = parser.parse_args()
    return args    


def perform_sample_sm_tabs(args):
    sample_details = ''
    barcodes = []
    ampTabs = {}
    sampleSupMutTabs = {}
      
    logfile = os.path.join(utilities.get_log_dir(args.output),"SampleSuperMut"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('SAMPLE SUPER MUTANT TAB PROCESSING STARTED for %s', args.template)

    try :
        #gather the list of barcodes in this sample
        barcodeMap_fh = open(args.barcodeMap,'r')
        for line in barcodeMap_fh:
            b = utilities.BarcodeMapRecord(*line.strip().split('\t'))
            if b.template == args.template:
                barcodes.append(b.barcode)
                sample_details = (args.template, b.purpose, b.gEsWellOrTotalULUsed, b.ampMatchName)
        barcodeMap_fh.close()

        #aggregate by primer all the WellAmpTabs for the barcodes in the sample
        ampTabs = load_primers(args, barcodes)

        #aggregate the mutants across the barcodes in the sample
        sampleSupMutTabs = load_chgs_by_sample(args, barcodes, ampTabs)
        
        #write the dictionary to output file
        create_sample_supMut_tab_file(args, sample_details, sampleSupMutTabs, ampTabs, len(barcodes))

        logging.info('SAMPLE SUPER MUTANT TAB PROCESSING COMPLETED for %s', args.template)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)

def load_primers(args, barcodes):
    ampTabs = {}
    
    directory = os.path.dirname(args.output)
    
    #for all barcodes in this sample, aggregate the the WellAmpTab file info on primer
    for barcode in barcodes:
        wat_filename = os.path.join(directory, os.pardir, "mutantTabs", barcode + ".wat")
          
        wat_fh = open(wat_filename,'r')
        for line in wat_fh:
            w = utilities.WellAmpTabRecord(*line.strip().split('\t'))
            if w.primer not in ampTabs: # initialize with 1 BC cnt; min, avg, max cnt and 1 cnt for Total UIDs and sumGoodReads
                ampTabs[w.primer] = [1, int(w.total_UIDs), int(w.total_UIDs),int(w.total_UIDs), 1, int(w.sumGoodReads), int(w.sumGoodReads), int(w.sumGoodReads), 1]
            else:  # if the primer has been encountered, aggregate totals
                ampTabs[w.primer][0] += 1 # add one to barcodes with primer count
                #Total_UIDs min, avg, max, cnt 
                ampTabs[w.primer][1:5] = aggregate_one(ampTabs[w.primer][1:5], w.total_UIDs)
                #sumGoodReadss min, avg, max, cnt 
                ampTabs[w.primer][5:9] = aggregate_one(ampTabs[w.primer][5:9], w.sumGoodReads)
        wat_fh.close() 
 
    return(ampTabs)


def load_chgs_by_sample(args, barcodes, ampTabs):
    sampleSupMutTabs = {}
    #aggregate WellSuperMutant Tabulation scores by primer across the barcodes in the sample
    directory = os.path.dirname(args.output)
    
    for barcode in barcodes:
        wsmt_filename = os.path.join(directory, os.pardir, "mutantTabs", barcode + ".wsmt")
          
        wsmt_fh = open(wsmt_filename,'r')
        for line in wsmt_fh:
            t = line.strip().split('\t')
            #key is tuple of change details: result is a list of stats
            key = (t[5], t[7], t[8], t[9], t[10], t[11])
            if key not in sampleSupMutTabs:
                sampleSupMutTabs[key] = [1, int(t[12]), float(t[12]), int(t[12]), 1, int(t[13]), int(t[14])] + \
                                        [int(t[15]), float(t[16]), int(t[17]), 1] + \
                                        [int(t[18]), float(t[19]), int(t[20]), 1] + \
                                        [int(t[21]), float(t[22]), int(t[23]), 1] + \
                                        [int(t[24]), float(t[25]), int(t[26]), 1] + \
                                        [int(t[27]), float(t[28]), int(t[29]), 1] + \
                                        [int(t[30]), float(t[31]), int(t[32]), 1] + \
                                        [int(t[33]), float(t[34]), int(t[35]), 1] + \
                                        [int(t[36]), float(t[37]), int(t[38]), 1] + \
                                        [float(t[41]), float(t[41]), float(t[41]), 1,] + \
                                        [float(t[42]), float(t[42]), float(t[42]), 1,] + \
                                        [int(t[39]), float(t[39]), int(t[39]), 1,] + \
                                        [int(t[40]), float(t[40]), int(t[40]), 1,] + \
                                        [t[43]]
            else:
                #aggregate this change instance with all others of its type for this primer
                sampleSupMutTabs[key][0] += 1       #increment the Mutant Well count (wells with this change)

                #aggregate the Distinct UID counts
                sampleSupMutTabs[key][1:5] = aggregate_one(sampleSupMutTabs[key][1:5], t[12])
                #aggregate the mutant read counts
                sampleSupMutTabs[key][5] += int(t[13])      #aggregate the Sum Mutant Read Counts
                sampleSupMutTabs[key][7:11] = aggregate_lists(sampleSupMutTabs[key][7:11], t[15:18])
                                                                
                #aggregate the family_good_read counts
                sampleSupMutTabs[key][6] += int(t[14])      #aggregate the Sum Family Good Reads Counts
                sampleSupMutTabs[key][11:15] = aggregate_lists(sampleSupMutTabs[key][11:15], t[18:21])
                
                #mismatch counts min, avg, max, cnt 
                sampleSupMutTabs[key][15:19] = aggregate_lists(sampleSupMutTabs[key][15:19], t[21:24])
                #corrected mismatch counts min, avg, max, cnt
                sampleSupMutTabs[key][19:23] = aggregate_lists(sampleSupMutTabs[key][19:23], t[24:27])
                #indel counts min, avg, max, cnt
                sampleSupMutTabs[key][23:27] = aggregate_lists(sampleSupMutTabs[key][23:27], t[27:30])
                #inserted bases min, avg, max, cnt
                sampleSupMutTabs[key][27:31] = aggregate_lists(sampleSupMutTabs[key][27:31], t[30:33])
                #deleted bases min, avg, max, cnt
                sampleSupMutTabs[key][31:35] = aggregate_lists(sampleSupMutTabs[key][31:35], t[33:36])
                #quality scores min, avg, max, cnt
                sampleSupMutTabs[key][35:39] = aggregate_lists(sampleSupMutTabs[key][35:39], t[36:39])

                #aggregate the % Mutant By All Reads
                sampleSupMutTabs[key][39:43] = aggregate_float(sampleSupMutTabs[key][39:43], t[41])
                #aggregate the % Mutant By UID Families
                sampleSupMutTabs[key][43:47] = aggregate_float(sampleSupMutTabs[key][43:47], t[42])
                #aggregate the Sum Good Reads
                sampleSupMutTabs[key][47:51] = aggregate_one(sampleSupMutTabs[key][47:51], t[39])
                #aggregate the Total Uids
                sampleSupMutTabs[key][51:55] = aggregate_one(sampleSupMutTabs[key][51:55], t[40])
    
        wsmt_fh.close()

    return(sampleSupMutTabs)


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


def aggregate_float(scores, number):
    #convert new number to an integer
    number = float(number)
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


def create_sample_supMut_tab_file(args, sample_details, sampleSupMutTabs, ampTabs, wells_per_sample):
    #assemble collected tabulation info into format for output file
    output_fh = open(args.output,'w')
    #write a header line
    output_fh.write('\t'.join(['Template', 'Purpose', 'GEsWellOrTotalULUsed', 'WellAmpMatchName', 'PrimerAmpMatchName', 'SampleWells', \
                               'SampleWellsPositiveForAmplimer', 'SumTotalUIDsForAmplimer', 'SumGoodReadsForAmplimer', 'MinTotalUIDsForAmplimerPerWell', \
                               'AverageTotalUIDsForAmplimerPerWell', 'MaxTotalUIDsForAmplimerPerWell', 'MinSumGoodReadsForAmplimerPerWell', 'AverageSumGoodReadsForAmplimerPerWell', \
                               'MaxSumGoodReadsForAmplimerPerWell', 'Chrom', 'Position', 'MutType', 'BaseFrom', 'BaseTo', 'MutantWellCount', \
                               'MutantDistinctUidCount', 'SumMutReadCounts', 'SumFamilyGoodReadCountFromMutantFamilies', 'MinMutantDistinctUidCountPerWell', \
                               'AverageMutantDistinctUidCountPerWell', 'MaxMutantDistinctUidCountPerWell', 'MinMutCount', 'AverageMutCount', 'MaxMutCount', \
                               'MinGoodReadCountInMutantWells', 'AverageFamilyGoodReadCountInMutantWells', 'MaxFamilyGoodReadCountInMutantWells', \
                               'MinMismatchCountInMutantReads', 'AverageMismatchCountInMutantReads', 'MaxMismatchCountInMutantReads', \
                               'MinCorrectedMismatchCountInMutantReads', 'AverageCorrectedMismatchCountInMutantReads', 'MaxCorrectedMismatchCountInMutantReads', \
                               'MinIndelCountInMutantReads', 'AverageIndelCountInMutantReadsInMutantReads', 'MaxIndelCountInMutantReadsInMutantReads', \
                               'MinInsertedBasesInMutantReads', 'AverageInsertedBasesInMutantReads', 'MaxInsertedBasesInMutantReads', \
                               'MinDeletedBasesInMutantReads', 'AverageDeletedBasesInMutantReads', 'MaxDeletedBasesInMutantReads', \
                               'MinQualityScoreForMutation', 'AverageQualityScoreForMutation', 'MaxQualityScoreForMutation', \
                               'Min%MutantByAllReadsByWellInMutantWells', 'Average%MutantByAllReadsByWellInMutantWells', 'Max%MutantByAllReadsByWellInMutantWells', \
                               'Min%MutantByUidFAmiliesByWellInMutantWells', 'Average%MutantByUidFAmiliesByWellInMutantWells', 'Max%MutantByUidFAmiliesByWellInMutantWells', \
                               'MinSumGoodReadsByWellInMutantWells', 'AverageSumGoodReadsByWellInMutantWells', 'MaxSumGoodReadsByWellInMutantWells', \
                               'MinTotalUIDsByWellInMutantWells', 'AverageTotalUIDsByWellInMutantWells', 'MaxTotalUIDsByWellInMutantWells', \
                               'SumGoodReadsFromMutantWells', 'TotalUIDsFromMutantWells', '%WellsWithMutation', '%MutantByAllReads', '%MutantByUidFAmilies', '%BackgroundRate']) + '\n')

    sample_details = list(sample_details)
        
    for change in sampleSupMutTabs:
        change_details = list(change)
        #distinct UID cnts min, max, sum, cnt
        distinct_UID_cnts = finalize_scores(sampleSupMutTabs[change][1:5])
        #mutant read cnts min, max, sum, cnt
        mut_cnts = finalize_scores(sampleSupMutTabs[change][7:11])
        #family_good_read_cnts min, max, sum, cnt
        fam_good_cnts = finalize_scores(sampleSupMutTabs[change][11:15])
        #mismatch_cnt min, max, sum, cnt
        mismatch_cnts = finalize_scores(sampleSupMutTabs[change][15:19])
        #corr_mismatch_cnt min, max, sum, cnt
        corr_mismatch_cnts = finalize_scores(sampleSupMutTabs[change][19:23])
        #indel_cnt min, max, sum, cnt
        indel_cnts = finalize_scores(sampleSupMutTabs[change][23:27])
        #ins_bases min, max, sum, cnt
        ins_bases = finalize_scores(sampleSupMutTabs[change][27:31])
        #del_bases min, max, sum, cnt
        del_bases = finalize_scores(sampleSupMutTabs[change][31:35])
        #quality scores min, max, sum, cnt
        quality_scores = finalize_scores(sampleSupMutTabs[change][35:39])
        #percent mutant by all reads min, max, sum, cnt
        mut_by_reads = finalize_scores(sampleSupMutTabs[change][39:43])
        #percent mutant by UID families min, max, sum, cnt
        mut_by_UIDs = finalize_scores(sampleSupMutTabs[change][43:47])
        #sum good reads by well min, max, sum, cnt
        good_reads_by_well = finalize_scores(sampleSupMutTabs[change][47:51])
        #total UIds by well min, max, sum, cnt
        UIDs_by_well = finalize_scores(sampleSupMutTabs[change][51:55])

        #use primer from the change details to access the total uid count and sum of good reads from amplicon tabulations
        primer = change[0]
        #sample_UID_cnts min, max, sum, cnt
        sample_UID_cnts = finalize_scores(ampTabs[primer][1:5])
        #sample Sum Good Reads cnts min, max, sum, cnt
        sample_SGR_cnts = finalize_scores(ampTabs[primer][5:9])
        sample_sums = [str(ampTabs[primer][0]), str(ampTabs[primer][2]), str(ampTabs[primer][6])]
        
        percent_wells_with_mut = 100*(sampleSupMutTabs[change][0]/wells_per_sample)
        percent_mut_by_reads = 100*(sampleSupMutTabs[change][5]/ampTabs[primer][6])
        percent_mut_by_UID_fams = 100*(sampleSupMutTabs[change][2]/ampTabs[primer][2])
        final_sums = [str(sampleSupMutTabs[change][48]), str(sampleSupMutTabs[change][52]),  \
                      str(percent_wells_with_mut), str(percent_mut_by_reads), str(percent_mut_by_UID_fams), \
                      str(sampleSupMutTabs[change][55])]

        output_fh.write('\t'.join(sample_details + [primer, str(wells_per_sample)] + sample_sums + \
                                  sample_UID_cnts + sample_SGR_cnts + change_details[1:6]  + \
                                  [str(sampleSupMutTabs[change][0]), str(sampleSupMutTabs[change][2]), str(sampleSupMutTabs[change][5]), str(sampleSupMutTabs[change][6])] + \
                                  distinct_UID_cnts + mut_cnts + fam_good_cnts + mismatch_cnts + corr_mismatch_cnts + indel_cnts + \
                                  ins_bases + del_bases + quality_scores + mut_by_reads + mut_by_UIDs + good_reads_by_well + UIDs_by_well +\
                                  final_sums) + '\n')
    output_fh.close()
    
    
def finalize_scores(scores):
    #scores is a list of Min, sum, Max, Cnt
    #convert scores into strings for output, calculate average from sum and count
    str_scores = []
    str_scores.append(str(scores[0])) #minimum
    str_scores.append(str(float(scores[1])/int(scores[3]) )) #average
    str_scores.append(str(scores[2])) #maximum

    return (str_scores)


def main():
     
    args = get_args()
    
    try:
        perform_sample_sm_tabs(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()