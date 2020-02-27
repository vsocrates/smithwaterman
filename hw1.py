#!/usr/bin/python
__author__ = "Vimig Socrates"
__email__ = "vimig.socrates@yale.edu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage:      python hw1.py -i <input file> -s <score file>
### Example:    python hw1.py -i input.txt -s blosum62.txt
### Note:       Smith-Waterman Algorithm

### Scripting must be done from scratch, without the use of any pre-existing packages.
### Python standard library (I/O) and numpy are allowed.

import argparse
import numpy as np
import pandas as pd
import sys 
import StringIO

### This is one way to read in arguments in Python.
### We need to read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

SUBSTITUTION_MAT_NUM_COLS = 24
SUBSTITUTION_MAT_NAMES = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y", "Z"] 

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
	### Print input and score file names. You can comment these out.
	print ("input file : %s" % inputFile)
	print ("score file : %s" % scoreFile)
	print ("open gap penalty : %s" % openGap)
	print ("extension gap penalty : %s" % extGap)

	with open(inputFile) as input_file:
		lines = input_file.read().splitlines()
		seq1 = lines[0]
		seq2 = lines[1]

	# Sequence 1 is the cols
	# Sequence 2 is the rows
	# print("seq2" , seq2)
	# print("seq1", seq1)

	sub_mat = np.loadtxt(scoreFile, skiprows=1, usecols=tuple(range(1, SUBSTITUTION_MAT_NUM_COLS)))
	sub_mat = pd.DataFrame(sub_mat, columns=SUBSTITUTION_MAT_NAMES, index=SUBSTITUTION_MAT_NAMES)
	
	scoring_mat = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=np.int32)
	traceback_mat = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=np.int32)

	LEFT, DIAGONAL, UP = range(3)

	def calculateGapPenalty(opening_penalty, extension_penalty, length):
		return opening_penalty + (extension_penalty * (length - 1))

	for i in range(1, len(seq2) + 1):

		for j in range(1, len(seq1) + 1):

			# diagonal score
			diag_score = scoring_mat[i - 1, j -1] + sub_mat.loc[seq2[i - 1], seq1[j - 1]]
			diag_score = (diag_score, DIAGONAL)

			# vertical score
			vert_score_with_new_gaps = scoring_mat[:i, j] + [calculateGapPenalty(openGap, extGap, length) for length in range(i, 0, -1)]
			if vert_score_with_new_gaps.size != 0:
				max_vert_score = max(scoring_mat[:i, j] + [calculateGapPenalty(openGap, extGap, length) for length in range(i, 0, -1)])
				max_vert_score = (max_vert_score, UP)
			else:
				max_vert_score = (0, UP)
			
			# horizontal score
			horz_score_with_new_gaps = scoring_mat[i, :j] + [calculateGapPenalty(openGap, extGap, length) for length in range(j, 0, -1)]
			if horz_score_with_new_gaps.size != 0:
				max_horz_score = max(scoring_mat[i, :j] + [calculateGapPenalty(openGap, extGap, length) for length in range(j, 0, -1)])
				max_horz_score = (max_horz_score, LEFT)
			else:
				max_horz_score = (0, LEFT)

			# find overall max
			# print("[diag_score, max_vert_score, max_horz_score, 0]", (seq2[i - 1], seq1[j - 1], [diag_score, max_vert_score, max_horz_score, 0]))
			scoring_mat[i, j], traceback_mat[i, j] = max([max_horz_score, diag_score, max_vert_score, (0, 0)])


	with open("smaple-traceback-example.txt", "w+") as testout:
		np.savetxt(testout, traceback_mat, fmt="%d", delimiter=" ")

	# backtrack time
	alignment_score = np.amax(scoring_mat)
	xTrace, yTrace = np.unravel_index(np.argmax(scoring_mat, axis=None), scoring_mat.shape)
	# print("maxs" , (alignment_score, xTrace, yTrace))
	# print("scoring_mat", scoring_mat[xTrace, yTrace])

	initXTrace = xTrace
	initYTrace = yTrace

	seq2_align = ""
	seq1_align = ""
	seq_identity = ""

	if seq1[yTrace - 1] == seq2[xTrace - 1]:
		seq2_align = seq2[xTrace - 1] + seq2_align
		seq1_align = seq1[yTrace - 1] + seq1_align
		seq_identity = "|" + seq_identity
		yTrace -= 1
		xTrace -= 1

	else:
		seq2_align = seq2[xTrace - 1] + seq2_align
		seq1_align = seq1[yTrace - 1] + seq1_align
		seq_identity = " " + seq_identity
		yTrace -= 1
		xTrace -= 1

	lastAlignIdx = 0

	while scoring_mat[xTrace, yTrace] != 0 and (xTrace, yTrace) != (0,0):

		trace_max_idx = traceback_mat[xTrace, yTrace]

		lastAlignIdx = trace_max_idx
		# go up: gap in seq1
		if trace_max_idx == UP:
			seq2_align = seq2[xTrace - 1] + seq2_align
			seq1_align = "-" + seq1_align
			seq_identity = " " + seq_identity
			xTrace -= 1

		# go diagonal
		elif trace_max_idx == DIAGONAL:
			if seq1[yTrace - 1] == seq2[xTrace - 1]:
				seq2_align = seq2[xTrace - 1] + seq2_align
				seq1_align = seq1[yTrace - 1] + seq1_align
				seq_identity = "|" + seq_identity
				yTrace -= 1
				xTrace -= 1
			else:
				seq2_align = seq2[xTrace - 1] + seq2_align
				seq1_align = seq1[yTrace - 1] + seq1_align
				seq_identity = " " + seq_identity
				yTrace -= 1
				xTrace -= 1

		# go left: gap in seq2
		elif trace_max_idx == LEFT:
			seq2_align = "-" + seq2_align
			seq1_align = seq1[yTrace - 1] + seq1_align 
			seq_identity = " " + seq_identity
			yTrace -= 1

	# 	print(trace_max_idx)
	# 	print("seq1", seq1_align)
	# 	print("mid1", seq_identity)
	# 	print("seq2", seq2_align)
	# 	print("(xTrace, yTrace)", (xTrace, yTrace))
	# 	print("(scoring_mat[xTrace, yTrace])", scoring_mat[xTrace, yTrace])

	# 	print("\n\n")

	# print("final xtrace", xTrace)
	# print("final ytrace", yTrace)
	# print("init xtrace", initXTrace)
	# print("init ytrace", initYTrace)
	# print("trace_max_idx", trace_max_idx)
	# figure out which seq is longer and add in the parantheses
	seq2_align = "(" + seq2_align + ")"
	seq1_align = "(" + seq1_align + ")"

	if trace_max_idx == 1:
		seq1_align = seq1[:yTrace] + seq1_align + seq1[initYTrace:]
		seq2_align = " " * (len(seq1[:yTrace])) + seq2[:xTrace] + seq2_align + seq2[initXTrace:]
		seq_identity = " " * (max([xTrace, yTrace]) + 1) + seq_identity + " " * max([len(seq1) - initYTrace, len(seq2) - initXTrace])
		# print("max([len(seq1) - initYTrace, len(seq2) - initXTrace]", max([len(seq1) - initYTrace, len(seq2) - initXTrace]))
	elif trace_max_idx == 2:
		seq1_align = seq1[:yTrace] + seq1_align + seq1[initYTrace - 1:]
		seq2_align = seq2[:xTrace] + seq2_align + seq2[initXTrace:]
		seq_identity = " " * max([xTrace, yTrace]) + seq_identity + " " * max([len(seq1) - initYTrace, len(seq2) - initXTrace])
	else:
		seq1_align = seq1[:yTrace] + seq1_align + seq1[initYTrace:]
		seq2_align = seq2[:xTrace] + seq2_align + seq2[initXTrace:]

		seq_identity = " " * max([xTrace, yTrace]) + seq_identity + " " * max([len(seq1) - initYTrace, len(seq2) - initXTrace])

	n = 100
	seq1_align = [seq1_align[i:i+n] for i in range(0, len(seq1_align), n)]
	seq2_align = [seq2_align[i:i+n] for i in range(0, len(seq2_align), n)]
	seq_identity = [seq_identity[i:i+n] for i in range(0, len(seq_identity), n)]

	# Write out the sequences, score matrix, sequences alignment in human-readable form, and score
	# OUTPUT_FILE = "test_output_mat.txt"
	OUTPUT_FILE = "output.txt"
	with open(OUTPUT_FILE, "w+") as out:

		# write out sequences
		out.write('----------------------------------------------------\n')
		out.write('|	Sequences                                       |\n')
		out.write('----------------------------------------------------\n\n')

		out.write("sequence1\t" + seq1 + "\n")

		out.write("sequence2\t" + seq2 + "\n\n")

		# write out score matrix
		out.write('----------------------------------------------------\n')
		out.write('|	Score Matrix                                    |\n')
		out.write('----------------------------------------------------\n\n')

		# print out the score matrix to a String IO object so we can modify it
		sio = StringIO.StringIO()
		np.savetxt(sio, scoring_mat, fmt="%d", delimiter="\t")
		mystr = sio.getvalue()
		mat_list = mystr.splitlines()

		# write seq1 as headers
		out.write("\t\t" + "\t".join(seq1) + "\t\n")

		# write the 0 row
		out.write("\t" + mat_list[0] + "\t\n")

		# write the rest of the rows with seq2 added as row names
		for idx, line in enumerate(mat_list[1:]):
			out.write(seq2[idx] + "\t" + line + "\t\n")

		out.write("\n")

		# write out the alignment and score
		out.write('----------------------------------------------------\n')
		out.write('|	Best Local Alignment                            |\n')
		out.write('----------------------------------------------------\n\n')

		out.write("Alignment Score:\t" + str(alignment_score) + "\n")

		for seq1, seq_iden, seq2 in zip(seq1_align, seq_identity, seq2_align):
			out.write(seq1 + "\n")
			out.write(seq_iden + "\n")
			out.write(seq2 + "\n")

		out.write("\n")


### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, int(args.opengap), int(args.extgap))

