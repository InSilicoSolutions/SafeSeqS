import argparse
import os
import sys, traceback
import logging
import zlib
import utilities
from utilities import PrimerRecord

#import cProfile

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


#create a global dictionary to hold the primer info from the primer file
primer_dict = {}
        
 
#Read the command line arguments and return them in args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='File with reads', required=True)
    parser.add_argument('-o', '--output', help='The output file name', required=True)
    parser.add_argument('-f', '--family', help='The well family file name', required=True)
    parser.add_argument('-p', '--primerset', help='The primer set file name.', required=True)

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
    global primer_dict
    
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

        #load primers from file into dictionary
        primer_dict = utilities.load_primers(args.primerset)
        
        unique_dict = {}
        counter = 0
        for line in input_fh:
            counter += 1
            split_line = line.split('\t')
            #compress the read to use it as key for dictionary
            key = zlib.compress (split_line[1].encode('utf-8'),9) 
            if key not in unique_dict:
                #if this is the first we have seen of this read, create a dictionary entry to hold all of its families
                unique_dict[key] = {}

            #set initial value to 1, then increment it on future instances of the UID within the same read sequence
            if split_line[5] not in unique_dict[key]:
                unique_dict[key][split_line[5]] = 1
            else:
                unique_dict[key][split_line[5]] += 1
                
            if counter % 100000 == 0:
                logging.debug('processed %s', counter)
                
        input_fh.close()
        logging.debug('total reads: %s', str(counter))

        unique_id = 0 # begin a unique counter to use as identifier for the unique read
        debug_family_cnt = 0
        for key in unique_dict:
            unique_id += 1
            read = zlib.decompress(key).decode("utf-8") 
            read_count = 0
            for uid in unique_dict[key]:
                debug_family_cnt += 1
                family_fh.write('\t'.join([str(unique_id), read, barcode, uid, str(unique_dict[key][uid])]) + '\n') 
                #count for this read will be sum of all family counts within it 
                read_count = read_count + unique_dict[key][uid]
                
            # initialize file output for No Match scenario
            primer = 'No Match'
            r1_match_type = 'No Match'
            r2_match_type = 'No Match'
            r1_primer_pos = 0
            r2_primer_pos = 0

            p, r1_match_type, r1_primer_pos = find_read1(read)
            if r1_match_type == 'Perfect Match':
                primer = p.ampMatchName
                #if there is a perfect match for the read1, look for the read2 AFTER read1
                r2_primer_pos = read[len(p.read1):len(read)].find(p.read2)
                if r2_primer_pos == -1: # read2 was not found, reset position
                    r2_primer_pos = 0
                else: #read2 was found
                    r2_match_type = 'Perfect Match'
                    r2_primer_pos = r2_primer_pos + len(p.read1) + 1 #add length of Read1, increment python string counter by one to reflect actual position of read2

            output_fh.write('\t'.join([str(unique_id), read, str(read_count), primer, r1_match_type, str(r1_primer_pos), r2_match_type , str(r2_primer_pos)]) + '\n')
 
        output_fh.close()
        
        logging.debug('unique reads: %s', str(unique_id))
        logging.debug('families: %s', str(debug_family_cnt))
        
        logging.info('UNIQUE PROCESSING COMPLETED for %s', args.input)
        return
    
    except Exception as err:
        logging.exception(err)
        traceback.print_exc()
        sys.exit(1)

        
#compress the read in 3 character chunks using the dictionary. don't read past end of read string. dictionary handles spaces at he end of read.
def compress_read(read):
    comp_dic = {'aaa': 'A',
                'aat': 'B',
                'aac': 'C',
                'aag': 'D',
                'ata': 'E',
                'att': 'F',
                'atc': 'G',
                'atg': 'H',
                'aca': 'I',
                'act': 'J',
                'acc': 'K',
                'acg': 'L',
                'aga': 'M',
                'agt': 'N',
                'agc': 'O',
                'agg': 'P',
                'taa': 'Q',
                'tat': 'R',
                'tac': 'S',
                'tag': 'T',
                'tta': 'U',
                'ttt': 'V',
                'ttc': 'X',
                'ttg': 'Y',
                'tca': 'Z',
                'tct': 'a',
                'tcc': 'b',
                'tcg': 'c',
                'tga': 'd',
                'tgt': 'e',
                'tgc': 'f',
                'tgg': 'g',
                'caa': 'h',
                'cat': 'i',
                'cac': 'j',
                'cag': 'k',
                'cta': 'l',
                'ctt': 'm',
                'ctc': 'n',
                'ctg': 'o',
                'cca': 'p',
                'cct': 'q',
                'ccc': 'r',
                'ccg': 's',
                'cga': 't',
                'cgt': 'u',
                'cgc': 'v',
                'cgg': 'w',
                'gaa': 'x',
                'gat': 'y',
                'gac': 'z',
                'gag': '1',
                'gta': '2',
                'gtt': '3',
                'gtc': '4',
                'gtg': '5',
                'gca': '6',
                'gct': '7',
                'gcc': '8',
                'gcg': '9',
                'gga': '!',
                'ggt': '~',
                'ggc': '#',
                'ggg': '$',
                'aa': '.',
                'at': '^',
                'ac': '&',
                'ag': '*',
                'ta': '(',
                'tt': ')',
                'tc': '_',
                'tg': '+',
                'ca': '|',
                'ct': '-',
                'cc': '=',
                'cg': ',',
                'ga': '>',
                'gt': '[',
                'gc': ']',
                'gg': '{',
                'a': '}',
                't': '?',
                'c': ':',
                'g': '<',
                }

    read = read.lower()
    packed_read = ""
    
    for i in range(0,len(read),3):
        if (i+3) < len(read):
            key = read[i:i+3]
        else:
            key = read[i:len(read)]
            
        packed_read = packed_read + comp_dic[key]

    return packed_read
    
        
#compress the read in 3 character chunks using the dictionary. don't read past end of read string. dictionary handles spaces at he end of read.
def decompress_read(string):
    comp_dic = {'A': 'aaa',
                'B': 'aat',
                'C': 'aac',
                'D': 'aag',
                'E': 'ata',
                'F': 'att',
                'G': 'atc',
                'H': 'atg',
                'I': 'aca',
                'J': 'act',
                'K': 'acc',
                'L': 'acg',
                'M': 'aga',
                'N': 'agt',
                'O': 'agc',
                'P': 'agg',
                'Q': 'taa',
                'R': 'tat',
                'S': 'tac',
                'T': 'tag',
                'U': 'tta',
                'V': 'ttt',
                'X': 'ttc',
                'Y': 'ttg',
                'Z': 'tca',
                'a': 'tct',
                'b': 'tcc',
                'c': 'tcg',
                'd': 'tga',
                'e': 'tgt',
                'f': 'tgc',
                'g': 'tgg',
                'h': 'caa',
                'i': 'cat',
                'j': 'cac',
                'k': 'cag',
                'l': 'cta',
                'm': 'ctt',
                'n': 'ctc',
                'o': 'ctg',
                'p': 'cca',
                'q': 'cct',
                'r': 'ccc',
                's': 'ccg',
                't': 'cga',
                'u': 'cgt',
                'v': 'cgc',
                'w': 'cgg',
                'x': 'gaa',
                'y': 'gat',
                'z': 'gac',
                '1': 'gag',
                '2': 'gta',
                '3': 'gtt',
                '4': 'gtc',
                '5': 'gtg',
                '6': 'gca',
                '7': 'gct',
                '8': 'gcc',
                '9': 'gcg',
                '!': 'gga',
                '~': 'ggt',
                '#': 'ggc',
                '$': 'ggg',
                '.': 'aa',
                '^': 'at',
                '&': 'ac',
                '*': 'ag',
                '(': 'ta',
                ')': 'tt',
                '_': 'tc',
                '+': 'tg',
                '|': 'ca',
                '-': 'ct',
                '=': 'cc',
                ',': 'cg',
                '>': 'ga',
                '[': 'gt',
                ']': 'gc',
                '{': 'gg',
                '}': 'a',
                '?': 't',
                ':': 'c',
                '<': 'g',
                }
    unpacked_read = ""
    
    for i in string:
        unpacked_read = unpacked_read + comp_dic[i]

    return unpacked_read.upper()
    
        
def find_read1(read):
    primer_match = PrimerRecord
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
#if __name__ == "__main__": cProfile.run('main()')