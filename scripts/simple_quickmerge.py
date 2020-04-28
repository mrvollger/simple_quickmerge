#!/usr/bin/env python
import argparse
import os 
import sys
import pysam
import pandas as pd
pd.set_option("display.width", 120)
COMPLEMENT = str.maketrans("ACTG", "TGAC")

# global var for inputs
args=None 
small = ["QUERY", "REF_START", "REF_END", "REF_LEN", "Q_START", "Q_END", "Q_LEN", "INNIE", "ORIENTATION", "OVERHANG"]


def intersect(a1, a2, b1, b2):
	inter = (a1 <= b2) and (a2 >= b1)
	if(inter):
		return(a2-b1)
	return(inter)

def check_all_intersections(group):
	for i in range( group.shape[0] ):
		for j in range(i+1, group.shape[0] ):
			inter = intersect(
						group.REF_START.iloc[i], group.REF_END.iloc[i], 
						group.REF_START.iloc[j], group.REF_END.iloc[j] )
			if(inter):
				return(True)
	return(False)

def get_frac_coverage(group, r):
	#length = r.get_reference_length(group["REF"].iloc[0])
	length = group.REF_END.max() - group.REF_START.min()
	covered = sum(group.REF_END - group.REF_START  )

	#for query, qgroup in group.groupby("QUERY"):
	#	qlength = sum(qgroup.Q_END - qgroup.QSTART)
	#	qcoveraed = qlength / qgroup.
	#	sys.stderr.write(query, qcovered)

	return(covered/length)


def get_query(row, q, RC):
	seq = q.fetch(row.QUERY)
	START = row.Q_START
	END = row.Q_END
	
	if(RC):
		seq = seq.translate(COMPLEMENT)[::-1]
		START = row.Q_LEN - START
		END = row.Q_LEN - END
			
	if(not row.INNIE and row.ORIENTATION == "R"):
		END = row.Q_LEN - 1
	if(not row.INNIE and row.ORIENTATION == "L"):
		START = 0
	
	return(seq[START:END])

def merge_pair(row, seq, q):
	RC = False
	if(row.Q_START > row.Q_END): RC=True

	if(row.INNIE): # repalce innner contigs with query sequence 
		seq = seq[0:row.REF_START] + get_query(row, q, RC) + seq[row.REF_END:]
	elif(row.ORIENTATION == "R"):
		seq = seq[0:row.REF_START] + get_query(row, q, RC)
	elif(row.ORIENTATION == "L"):
		seq = get_query(row, q, RC) + seq[row.REF_END:]
	else:
		raise("Should be either innie L or R")

	return(seq)

def make_merge(group, r, q):
	# limit REF seq to be at most the edges of the QUERY sequence
	group = group.copy()
	MIN = group.REF_START.min()
	MAX = group.REF_END.max()
	group.REF_START = group.REF_START - MIN
	group.REF_END = group.REF_END - MIN

	seq = r.fetch(group["REF"].iloc[0])[MIN:MAX+1]
	for idx, row in group.iterrows():
		seq = merge_pair(row, seq, q)
	return(seq)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("tbl", help="positional input")
	parser.add_argument("reference", help="positional input")
	parser.add_argument("query", help="positional input")
	parser.add_argument("-o", help="output fasta file", required=True)
	parser.add_argument("--ml", help="minimum overlap length (TODO)", type=int, default=5000)
	parser.add_argument("--mincoverage", help="numeric option", type=float, default=0.98)
	parser.add_argument("-l", help="minimum length for a reference anchor contig", type=int, default=10*10**6)
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	
	#REF QUERY REF_START REF_END Q_START Q_END ORIENTATION INNIE(1/0) OVERLAP_LEN OVERLAP_PROP NO_OVERLAP_AT_ENDS OVERHANG
	merge = pd.read_csv(args.tbl, sep=r"\s+")
	merge.sort_values(by=['REF', 'REF_END'], ascending=[True, False], inplace=True)
	merge.rename(columns={"INNIE(1/0)":"INNIE"}, inplace=True)
	# make zero based coords
	for col in ["REF_START", "REF_END", "Q_START", "Q_END"]:
		merge[col] = merge[col] - 1

	# read in fasta files
	r = pysam.FastaFile(args.reference)
	q = pysam.FastaFile(args.query)
	
	# add contig lengths 
	merge["REF_LEN"] = merge.REF.map(r.get_reference_length)
	merge["Q_LEN"] = merge.QUERY.map(q.get_reference_length)
	
	# remove references that are internal to the query OVERHANG==-1
	merge = merge.loc[ merge.OVERHANG != -1 ]

	out = open(args.o, "w+")
	USED = set() # these contigs from the query have now been merged
	counter = 1
	for ref, group in merge.groupby("REF"):
		#print(group[small])
		NINNIE = sum(group.INNIE)
		N = group.shape[0]
		cov = get_frac_coverage(group, r) 
		length = r.get_reference_length(ref)
		
		if( length < args.l):
			if(args.d): sys.stderr.write(f"Skipped {ref} because {length} is less than {args.l} the min length.\n")
		elif( N-NINNIE > 2  or  check_all_intersections(group) ):
			if(args.d): sys.stderr.write(f"Skipped {ref} because of overlapping query contigs.\n")
		elif( N < 2):
			if(args.d): sys.stderr.write(f"Skipped {ref} because it is not merging two query contigs.\n")
		elif( cov < args.mincoverage ):
			if(args.d): sys.stderr.write(f"Skipped {ref} because {cov:.2f} is less than {args.mincoverage} coverage.\n")
		else: # these are the valid overlaps
			USED = USED.union(group.QUERY)	
			seqns = ",".join(group.QUERY.iloc[::-1])
			sys.stderr.write("\rMaking merge of: " + seqns)

			name = ">merged{:08}\tmerged:".format(counter)  + seqns  + "\n"
			seq = make_merge(group, r, q)
			out.write(name + seq + "\n")
			#if(counter == 5): break
			counter += 1

	sys.stderr.write("\n")

	for _, name in sorted(zip(q.lengths, q.references), reverse=True):
		if(name not in USED):
			out.write( ">{}\n{}\n".format(name, q.fetch(name)) )

	out.close()
