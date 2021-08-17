from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import subprocess
import linecache
import argparse
import sys
import os

def hmmscan_type(hmmscan_infile, path_hmm, profile_name):
	"""
	Runs HMMER hmmscan.
	hmmscan searches each sequence against the A and B merged profile.
	"""

	hmmscan_out_base = os.path.splitext(hmmscan_infile)[0] # FASTA file containing trimmed variants of the input sequence(s) w/o extension.

	# Use N-terminal Domain HMM Profile to find and type N-domain sequence
	Type_file = hmmscan_out_base + "_Type.out" # Name for hmmscan table output file
	subprocess.call(['hmmscan --noali -o '+hmmscan_out_base+'_Type_hmmscan.txt --domE 0.00001 --domtblout '+Type_file+' '+path_hmm+profile_name+'.hmm '+hmmscan_infile+''], shell=True)

	hmmscan_out_files = []

#	for file in os.listdir(os.getcwd()):
#		if file.endswith("_Type.out"):
#			hmmscan_out_files.append(file)
#	if len(hmmscan_out_files) != 1:
#		sys.exit("HMMER hmmscan output files were not created or more output files were found than created.")

def Nterm_trim_hmm(hmmscan_out, original_seq_file):
	"""
	parses HMMscan output, creates dictionaries of most probable N-term sequence.
	Reports all VSG sequences for which N-termini could not be found.
	"""
	#hmmscan.out files report each hit in ranked order, with best scoring hit first
	#we noticed that hmmscan very often reported two high quality domain hits when set to an e-value threshold of 1e-5
	#these hits often overlapped. therefore, we define the N-domain boundary by the rightmost coordinate
	type_dict = {}
	with open(hmmscan_out, 'r') as f:
	    for line in f.readlines():
	        filler = "#"
	        if filler not in line:
	            rows = line.split()
	            VSG = rows[3]
	            length = int(rows[5])
	            Profile = str(rows[0])
	            domscore = float(rows[13])
	            coord = int(rows[20])
	            if VSG in type_dict:
	                if coord > type_dict[VSG][1]:
	                    type_dict[VSG] = [length, coord, Profile]
	            else:
	                type_dict[VSG] = [length, coord, Profile]

	seqs = {}
	with open(original_seq_file, "r") as s: #This is the original fasta that you also plugged into hmmscan.
		for record in SeqIO.parse(s, "fasta"):
			ID = str(record.id)
			seq = str(record.seq)
			seqs[ID] = seq

	trim_dict = {}
	fail_dict = {}
	for ID, seq in seqs.items():
		if ID in type_dict:
			trim = seq[:type_dict[ID][1]]
			trim_dict[ID] = trim
		else:
			fail_dict[ID] = seq
	return trim_dict, type_dict, fail_dict

def hmmscan_subtype(hmmscan_infile, path_hmm):
	"""
	Runs HMMER hmmscan.
	hmmscan searches each sequence variant against the VSG N-terminal Type A(1-3) and B(1-2) HMM profile databases.
	"""

	hmmscan_out_base = os.path.splitext(hmmscan_infile)[0] # FASTA file containing the identified N-terminal domain(s).

	# TypeA1-3 and TypeB1-2 N-terminal Domain HMM Profiles
	ABsubtype_file = hmmscan_out_base + "_TypeABsubtype.out" # Name for hmmscan table output file (TypeA1-3/TypeB1-2).
	subprocess.call(['hmmscan --noali -o '+hmmscan_out_base+'_TypeABsubtype_hmmscan.txt --tblout '+ABsubtype_file+' '+path_hmm+'VSG-N-TypeSubtype.hmm '+hmmscan_infile+''], shell=True)

	hmmscan_out_files = []

	for file in os.listdir(os.getcwd()):
		if file.endswith("_TypeABsubtype.out"):
			hmmscan_out_files.append(file)
	if len(hmmscan_out_files) != 1:
		sys.exit("HMMER hmmscan output file was not created or more output files were found than created.")

def nterm_subtype_ID(hmmscan_out):
	"""
	Identifies the most probable N-terminal domain subtype of the N-terminal TypeA and/or TypeB domain sequence(s).
	Returns a dictionary containing the N-terminal domain sequence ID(s) and the corresponding type-subtype.
	"""
	subtype_dict = {}
	with open(hmmscan_out, 'r') as f:
	    for line in f.readlines():
	        filler = "#"
	        if filler not in line:
	            rows = line.split()
	            VSG = rows[2]
	            subtype = str(rows[0])[-2:]
	            domscore = rows[7]
	            if not VSG in subtype_dict: #take only the first hmmscan hit in list, this is highest quality
	                subtype_dict[VSG] = subtype
	return subtype_dict


# Arguments entered on command line.
parser = argparse.ArgumentParser(description="Identify the N-terminal domain type and subtype of one or more VSG sequences.")
parser.add_argument('file', help="FASTA file containing VSG protein sequences.", action="store")
parser.add_argument('path', help="Path to directory containing VSG N-terminal TypeA/TypeB/TypeSubtype HMM profiles.", action="store")

args = parser.parse_args()
infile = args.file
path_hmm = args.path


# Confirm that HMMER is installed in PATH.
try:
	subprocess.check_output(["hmmscan","-h"])
except subprocess.CalledProcessError:
	sys.exit("HMMER not installed in PATH.")


infile_base = os.path.splitext(infile)[0] # Input filename as string w/o extension.

#sys.stdout = open(infile_base + "_Nterm_stdout.txt", 'w')  # write all unix output to txt file

if not os.path.exists(infile):
	sys.exit("%r does not exist." % infile)
else:
	print "%r found.\n" % infile
	if not os.path.exists(path_hmm):
		sys.exit("Path to HMM profiles does not exist.")
	else:
		if not path_hmm.endswith("/"):
			path_hmm += "/"

# halt pipeline if output already exists
nterms_final_file = infile_base + "_Nterms.fa"
if os.path.exists(nterms_final_file):
	sys.exit("%r already exists." % nterms_final_file)

print "\nProceeding to HMMER hmmscan analysis that will determine the type of the N-terminal domain(s)..."
hmmscan_type(infile, path_hmm, str('VSG-N-mergeAB'))
hmmscan_first_pass = os.path.splitext(infile)[0] + "_Type.out"
first_pass_seqs, first_pass_types, first_pass_fails = Nterm_trim_hmm(hmmscan_first_pass, infile)

print "\nFailed to determine N-terminal domain sequences in input file:"
print bool(first_pass_fails)
if bool(first_pass_fails) == True:
	print "\nProceeding to HMMER recovery hmmscan analysis on sequence(s) that failed first pass..."
	fail_file = os.path.splitext(infile)[0] + "_first_pass_fails.fa"
	with open(fail_file, "w") as fail:
		for ID, seq in first_pass_fails.items():
			fail.write(">"+ID+'\n'+str(seq)+'\n')
	hmmscan_type(fail_file, path_hmm, str('VSG-N-DomainBoundary'))
	hmmscan_second_pass = os.path.splitext(fail_file)[0] + "_Type.out"
	second_pass_seqs, second_pass_types, second_pass_fails = Nterm_trim_hmm(hmmscan_second_pass, fail_file)
	ultimate_fail = os.path.splitext(infile)[0] + "_undefinedNdomain.fa"
	with open(ultimate_fail, "w") as final_fail:
		for ID, seq in second_pass_fails.items():
			final_fail.write(">"+ID+'\n'+str(seq)+'\n')
	All_N_seqs = first_pass_seqs.copy()
	All_N_seqs.update(second_pass_seqs)
	All_N_types = first_pass_types.copy()
	All_N_types.update(second_pass_types)
	print "%r contains sequence(s) with undefined N-terminal domain(s).\n" % (infile_base+"_undefinedNdomain.fa")

else:
	All_N_seqs = first_pass_seqs
	All_N_types = first_pass_types

with open(nterms_final_file, "w") as final:
	for ID, seq in All_N_seqs.items():
		final.write(">"+ID+'\n'+str(seq)+'\n')

print "N-terminal domains identified.\n"

print "Proceeding to second HMMER hmmscan analysis that will determine the subtype of the identified N-terminal domain(s)..."
hmmscan_subtype(nterms_final_file, path_hmm)
print "HMMER hmmscan analysis complete.\n"

all_original_seqIDs = []

for sequence in SeqIO.parse(infile, "fasta"):
	all_original_seqIDs.append(sequence.id)

for file in os.listdir(os.getcwd()):
	if file.endswith("TypeABsubtype.out"):
		empty = "#" # If there are no hmmscan hits, then the fourth line in hmmscan .out file contains "#".
		fourth_line = str(linecache.getline(file, 4))
		if not empty in fourth_line:
			typed_nterms = nterm_subtype_ID(file) # Dictionary - {N-terminal domain sequence ID: type-subtype}.
			print "Found most probable type and subtype of the identified N-terminal domain(s) in %r.\n" % file

print "Preparing summary report..."

with open(infile_base+"_NtermSummary.csv", "w") as summary:
	summary.write('original_seqID,original_seq_length,nterm_seq_length,HMM_profile,nterm_typesubtype\n')
	for key in All_N_types:
		for ID in typed_nterms:
			if key == ID:
				summary.write(ID+','+str(All_N_types[ID][0])+','+str(All_N_types[ID][1])+','+str(All_N_types[ID][2])+','+typed_nterms[ID]+'\n')
	no_sub = {} # find all untyped seq IDs
	for diff in set(All_N_types) - set(typed_nterms): #the set difference, iDs present in All_N_types but not the subtype dictionary
		no_sub[diff] = All_N_types[diff]
		summary.write(diff+','+str(no_sub[diff][0])+','+str(no_sub[diff][1])+','+str(no_sub[diff][2])+',NA\n')
	if os.path.exists(os.path.splitext(infile)[0] + "_undefinedNdomain.fa"):
		with open(os.path.splitext(infile)[0] + "_undefinedNdomain.fa", "r") as undefined:
			for record in SeqIO.parse(undefined, "fasta"):
				ID = str(record.id)
				length = str(len(record.seq))
				summary.write(ID+','+length+',NA,None,NA\n')
	else:
		summary.close()

# this was never a me problem. the missing VSG are not making it into the nterm_subtype_ID dictionary but they ARE IN THE HMMSCAN OUTPUT!!!!
# I am angry

print "%r contains information about the identified N-terminal domain sequence(s).\n" % (infile_base+"_NtermSummary.csv")

#clean up output files
#Move all HMMER hmmscan input and output files to a subdirectory.
subprocess.call(['mkdir -p '+infile_base+"_HMMSCAN"], shell=True)
subprocess.call(['mv *.out *hmmscan.txt '+infile_base+"_HMMSCAN"], shell=True)
#move fasta raw sequence files to subdirectory
subprocess.call(['mkdir -p '+infile_base+"_Sequences"], shell=True)
subprocess.call(['mv *.fa '+infile_base+"_Sequences"], shell=True)

print "Finished with %r." % infile
