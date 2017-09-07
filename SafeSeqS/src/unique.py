import argparse
import os
import sys, traceback
import logging
import uuid

# unique          - This program takes file of fastq reads split by barcode and 
#                   counts the occurances of unique read sequences. The output is 
#                   a tab delimited file with one line per unique read sequence.
#
        
 
#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='File with reads', required=True)
    parser.add_argument('-o', '--output', help='The output file name', required=True)
    parser.add_argument('-f', '--family', help='The well family file name', required=True)

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


def perform_unique(args):
    
    logfile = os.path.join(get_log_dir(args.output),"Unique"+ str(os.getpid()) + ".log")
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
        input_fh = open(args.input,'r')
        output_fh = open(args.output,'w')
        family_fh = open(args.family,'w')
        unique_dict = {}
    
        for line in input_fh:
            split_line = line.split('\t')
            if split_line[1] not in unique_dict:
                unique_dict[split_line[1]] = {}
                
            if split_line[5] not in unique_dict[split_line[1]]:
                unique_dict[split_line[1]][split_line[5]] = 1
            else:
                unique_dict[split_line[1]][split_line[5]] += 1
                
        for read in unique_dict:
            read_count = 0
            for uid in unique_dict[read]:
                family_fh.write('\t'.join([str(uuid.uuid4()), read, barcode, uid, str(unique_dict[read][uid])]) + '\n')  
                read_count = read_count + unique_dict[read][uid]
            output_fh.write('\t'.join([str(uuid.uuid4()), read, str(read_count), 'No Match', 'No Match', '0']) + '\n')

    

        input_fh.close()
        output_fh.close()
        logging.info('UNIQUE PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
        

def main():
     
    args = get_args()
    
    try:
        perform_unique(args)

    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)
    
       
if __name__ == "__main__": main()