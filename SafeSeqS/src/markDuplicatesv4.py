import pyodbc
import math
import numpy as np
from scipy.cluster.hierarchy import fclusterdata
from tqdm import *
from joblib import Parallel, delayed
from collections import Counter
from itertools import izip
import time
import sys
import binascii

def distance(p0, p1):
	if p0[0] == p1[0] and p0[1] == p1[1] and p0[2] == p1[2] and p0[3] == p1[3] and p0[4] == p1[4]: # same instrument, run, flowcell, lane, tile
		return math.sqrt((p0[5] - p1[5])**2 + (p0[6] - p1[6])**2)
	else:
		return 5001
		
def cFC(coords):
	return max(fclusterdata(coords, 5000, criterion = 'distance', metric = distance))
	
def markDuplicates(coords):
	clusters = fclusterdata(coords, 5000, criterion = 'distance', metric = distance)
	counts = Counter(clusters)
	result = ['Unique' if counts[c] == 1 else 'Duplicate' for c in clusters]
	return result

db = sys.argv[1] # name of SQL database
cnxn = pyodbc.connect('DSN=ludseq8;uid=safe_seqs_processing;pwd=LudwigCluster1@@;database=' + db) # connect to DB
cnxn.setdecoding(pyodbc.SQL_CHAR, encoding='latin1', to=str)
cnxn.setencoding(str, encoding='latin1')
print('Connected to ' + db)

print('Fetching SM reads...')
start = time.time()
cursor = cnxn.cursor()
cursor.execute("SELECT DISTINCT A.IndexSequence, A.UID, A.FamilyGoodReadCount, B.*, C.* FROM [SuperMutantTabulations] A INNER JOIN [Reads] B ON A.[IndexSequence]=B.[IndexSequence] AND A.[UID]=B.[UID] INNER JOIN [SequenceAlignment] C ON B.[ReadSequence]=C.[ReadSequence] WHERE cast(A.MutCount as decimal)/cast(A.FamilyGoodReadCount as decimal) > 0.9 AND C.CorrectedMismatchCount <= 3 AND C.indelcount <= 1 AND charindex('Perfect', C.[Read1PrimerMatch]) > 0 AND charindex('Perfect', C.[Read2PrimerMatch]) > 0")
reads = cursor.fetchall()
end = time.time()
print('Fetched SM reads in ' + str(round(end - start)) + ' seconds.')

print('Organizing read information into UID families...')
# Read UID coordinates into memory for parallel processing--20 min
start = time.time()
UID_coords = {}
for read in tqdm(reads):
	wbc = read[0]
	uid = read[1]
	readid = read[3].split(' ')[0].split(':')
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
		UID_coords[uid + '-' + wbc].append(cluster)
	except:
		UID_coords[uid + '-' + wbc] = [cluster]

end = time.time()
print('Read organization completed in ' + str(round(end - start)) + ' seconds.')
				
print('Calculating corrected FC...')
start = time.time()
results = Parallel(n_jobs=-1, verbose=5)(delayed(cFC)(i) for i in UID_coords.values())
end = time.time()
print Counter(results)[1]/float(len(results))
print('Corrected FC completed in ' + str(round(end - start)) + ' seconds.')

print('Insert into new corrected FC table...')
start = time.time()
cursor.execute("CREATE TABLE CorrectedFCv2 (IndexSequence varchar(250), UID char(14), cFC int)")
cnxn.commit()

for i,j in izip(UID_coords.keys(), results):
	IDcFC = i.split('-')
	UIDcFC = IDcFC[0]
	WBCcFC = IDcFC[1]
	cursor.execute("INSERT INTO CorrectedFCv2 values (?, ?, ?)", WBCcFC, UIDcFC, int(j))

cnxn.commit()
cnxn.close()
end = time.time()
print('Inserted corrected FC values in ' + str(round(end - start)) + ' seconds.')
print('Done!')
