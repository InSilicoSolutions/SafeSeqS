import argparse
import os
import sys
import traceback
import logging
import binascii
import math
from scipy.cluster.hierarchy import fclusterdata
from safeseqs import utilities


#Read the command line arguments and return them in args
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Reads file.', required=True)
    parser.add_argument('-o', '--output', help='Corrected Family Counts File.', required=True)
    parser.add_argument('-g', '--goodreads', help='File containing compressed read seqs from bad reads. To be used by Optical Dups and UIDstats processes.', required=True)
    parser.add_argument('-d', '--distance', help='Minimum Distance between records that are NOt optical duplicates.', required=True)

    args = parser.parse_args()
    return args    


def load_good_readseqs(filename):
    good_reads = {}
    gr_fh = open(filename,'r')
    #loop through the list of good reads for the barcode and store them in a list
    for line in gr_fh:
        l = line.strip().split('\t')
        good_reads[l[1]] = 1
    gr_fh.close()    
    logging.info('Loaded %s good reads', str(len(good_reads)))


    return(good_reads)


def load_UID_coordinates(args, good_reads):
    UID_coords = {}
    gr_cnt = 0
    br_cnt = 0

    input_fh = open(args.input,'r')

    #read through raw read data, creating a dictionary of UIDs with their read header information
    for line in input_fh:
        r = utilities.ReadRecord(*line.strip().split('\t'))
        if r.read_seq in good_reads:    
            uid = r.uid
            readid = r.read_hdr.split(' ')[0].split(':')
            instrument = readid[0]
            instrument_int = int(binascii.hexlify(instrument.encode('utf-8')), 16) # convert instrument string to int
            run = int(readid[1])
            flowcell = readid[2]
            flowcell_int = int(binascii.hexlify(flowcell.encode('utf-8')), 16) # convert flowcell string to int
            lane = int(readid[3])
            tile = int(readid[4])
            x = int(readid[5])
            y = int(readid[6])
            cluster = (instrument_int, run, flowcell_int, lane, tile, x, y) # (Instrument, Run, Flowcell, Lane, Tile, X, Y)
            try:
                UID_coords[uid].append(cluster)
            except:
                UID_coords[uid] = [cluster]

            gr_cnt += 1
        else:
            br_cnt += 1

    input_fh.close()
    logging.info('Number of reads processed. good: %s vs. bad: %s' %(str(gr_cnt), str(br_cnt)))
    
    return(UID_coords) 


def distance(p0, p1):
    if p0[0] == p1[0] and p0[1] == p1[1] and p0[2] == p1[2] and p0[3] == p1[3] and p0[4] == p1[4]: # same instrument, run, flowcell, lane, tile
        return math.sqrt((p0[5] - p1[5])**2 + (p0[6] - p1[6])**2)
    else:
        return 100000
    
        
def cFC(coords, dist_parm):
    return max(fclusterdata(coords, dist_parm, criterion = 'distance', metric = distance))
    return 


def find_optical_duplicates(args):
    
    logfile = os.path.join(utilities.get_log_dir(args.output),"OptDups"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('OPTICAL DUPLICATES PROCESSING STARTED for %s', args.input)

    try :
        args.distance = int(args.distance)
        
        filename = os.path.split(args.input)[1]
        barcode = filename.split('.')[0]

        good_reads = load_good_readseqs(args.goodreads)
        UID_coords = load_UID_coordinates(args, good_reads)

        #now, step through the UIDs dictionary. calculate and write it's corrected family count
        output_fh = open(args.output,'w')

        for uid in UID_coords:

            family_good_cnt = len(UID_coords[uid])
            if family_good_cnt > 1:
                corrected_cnt = cFC(UID_coords[uid], args.distance)
            else:
                corrected_cnt = family_good_cnt 
            
            #only record the uids with optical duplicates
            if corrected_cnt != family_good_cnt:
                output_fh.write('\t'.join([barcode, uid, str(corrected_cnt), str(family_good_cnt - corrected_cnt)]) + '\n')
     
        output_fh.close()

        logging.info('OPTICAL DUPLICATES PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)


def main():
     
    args = get_args()
    
    try:
        find_optical_duplicates(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()