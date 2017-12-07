import argparse
import os
import sys
import traceback
import logging
import zlib
from safeseqs import utilities

# unique     This program takes file of fastq reads that have been split by barcode and  
#            compresses it further into a list of the unique read sequences. The output 
#            is a tab delimited file with one line per unique read sequence with a counter
#            for how many records were compressed.
#
#            It performs the Read1 and Read2 searches. It only looks for and reports on 
#            perfect matches. Read1 must be in the first position of the read sequence. If 
#            Read1 is found, it searches for Read2. 
#
#            It also counts the Well Families for future reporting use later.
#

 
#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='File with reads', required=True)
    parser.add_argument('-o', '--output', help='File containing unique read sequences for the barcode.', required=True)
    parser.add_argument('-f', '--family', help='The well family file name', required=True)
    parser.add_argument('-p', '--primerset', help='The primer set file name.', required=True)

    args = parser.parse_args()
    return args    


def perform_unique(args):
    
    logfile = os.path.join(utilities.get_log_dir(args.output),"Unique"+ str(os.getpid()) + ".log")
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s', 
        datefmt='%m/%d/%y %I:%M:%S %p',
        filename=logfile,
        filemode = 'w',
        level=logging.DEBUG)
        
    logging.info('UNIQUE PROCESSING STARTED for %s', args.input)

    try :
        filename = os.path.split(args.input)[1]
        barcode = filename.split('.')[0]
        
        #load primers from file into dictionary
        primer_dict = utilities.load_primers(args.primerset)
        
        unique_dict = load_reads(args)
        
        output_fh = open(args.output,'w')
        family_fh = open(args.family,'w')

        unique_id = 0 # begin a unique counter to use as identifier for the unique read
        for key in unique_dict:
            unique_id += 1
            read = zlib.decompress(key).decode("utf-8") 
            #look for primer match for this unique read sequence (well family will hold primer for UIDstats)
            # initialize file output for No Match scenario
            primer = 'No Match'
            r1_match_type = 'No Match'
            r2_match_type = 'No Match'
            r1_primer_pos = 0
            r2_primer_pos = 0

            p, r1_match_type, r1_primer_pos = find_read1(read, primer_dict)
            if r1_match_type == 'Perfect Match':
                primer = p.ampMatchName
                #if there is a perfect match for the read1, look for the read2 AFTER read1
                r2_primer_pos = read[len(p.read1):len(read)].find(p.read2)
                if r2_primer_pos == -1: # read2 was not found, reset position
                    r2_primer_pos = 0
                else: #read2 was found
                    r2_match_type = 'Perfect Match'
                    r2_primer_pos = r2_primer_pos + len(p.read1) + 1 #add length of Read1, increment python string counter by one to reflect actual position of read2

            #loop through the UIDs for this read sequence to get well families, and sum total read count            
            read_count = 0
            for uid in unique_dict[key]:
                family_fh.write('\t'.join([str(unique_id), read, barcode, uid, str(unique_dict[key][uid]), primer]) + '\n') 
                #count for this read will be sum of all family counts within it 
                read_count = read_count + unique_dict[key][uid]
                
            output_fh.write('\t'.join([str(unique_id), read, str(read_count), primer, r1_match_type, str(r1_primer_pos), r2_match_type , str(r2_primer_pos)]) + '\n')
 
        output_fh.close()
        
        logging.info('UNIQUE PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)

        
def load_reads(args):
    #load reads for this barcode into a dictionary to group them by read sequence and UID
    unique_dict = {}
    input_fh = open(args.input,'r')
    counter = 0
    for line in input_fh:
        counter += 1
#       ReadRecord = namedtuple('ReadRecord', ['read_hdr', 'read_seq', 'read_qual', 'barcode', 'bc_quality', 'uid', 'uid_qual'])
        r = utilities.ReadRecord(*line.strip().split('\t'))
        #compress the read to use it as key for dictionary
        key = zlib.compress (r.read_seq.encode('utf-8'),9) 
        if key not in unique_dict:
            #if this is the first we have seen of this read, create a dictionary entry to hold all of its families
            unique_dict[key] = {}

        #set initial value to 1, then increment it on future instances of the UID within the same read sequence
        if r.uid not in unique_dict[key]:
            unique_dict[key][r.uid] = 1
        else:
            unique_dict[key][r.uid] += 1
            
        if counter % 1000000 == 0:
            logging.debug('processed %s', counter)
            
    input_fh.close()
    logging.debug('total reads processed: %s', str(counter))
    
    return(unique_dict)


def find_read1(read, primer_dict):
    primer_match = utilities.PrimerRecord
    r1_match_type = 'No Match'
    r1_primer_pos = 0

    for primer in primer_dict:
        p = primer_dict[primer]
        #look for perfect match
        if read.startswith(p.read1):
            primer_match = p
            r1_match_type = 'Perfect Match'
            r1_primer_pos = 1
            #stop at first primer match
            break
            
    return primer_match, r1_match_type, r1_primer_pos
    
    
def main():
     
    args = get_args()
    
    try:
        perform_unique(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()