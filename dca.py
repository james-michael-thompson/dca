#!/usr/bin/python

import os
import re
import sys
import pickle

# free parameters
filename = 'data/T0594.fasta.chkali-5.sto'
alphabet = '-ARNDCEQGHILKMFPSTWYV'
id_cutoff = 0.8
pseudo_lambda = 10

# assertions on inputs
assert( id_cutoff > 0 and id_cutoff <= 0.8 )
assert( os.path.exists(filename) )

def is_aa_sequence(str):
	aa_regex = '^([%s]+)$' % alphabet
	m = re.match(aa_regex, str)
	return m != None

def read_MSA(filename):
	# create list of strings representing the aligned sequences
	MSA = []
	for line in open(filename,'r'):
		if line[0] == '#': continue
		line   = line.rstrip()
		tokens = line.split()
		if len(tokens) == 2 and is_aa_sequence(tokens[1]): MSA.append(tokens[1])
	return MSA

def pairwise_seq_id(MSA,idx1,idx2):
	assert( len(MSA[idx1]) == len(MSA[idx1]) )
	ident = 0
	for idx in xrange(0,len(MSA[idx1])):
		if MSA[idx1][idx] == MSA[idx2][idx]:
			ident += 1
	return float(ident) / len(MSA[idx1])

def calc_seq_weights(MSA,seq_idx):
	N = 0
	for idx in xrange(0,len(MSA)):
		idx1 = min(idx,seq_idx)
		idx2 = max(idx,seq_idx)
		if pairwise_seq_id(MSA,idx1,idx2) >= id_cutoff: N += 1
	return 1/float(N)

def calc_seq_weights(MSA):
	nrow = len(MSA)
	M = [ 0 for i in xrange(nrow) ]
	for i in xrange(nrow):
		for j in xrange(i,nrow): # intentionally calculate numbers for diagonal
			if i == j:
				M[i] += 1
			else:
				N_ident = pairwise_seq_id(MSA,i,j)
				if N_ident > id_cutoff:
					M[i] += 1
					M[j] += 1
	return [ (1.0/M[i]) for i in xrange(nrow) ]

def norm_freq(sym,pos,MSA,weights,lam,alphabet_size):
	Meff  = sum(weights)
	total = lam/alphabet_size
	for seq_idx in xrange(len(MSA)):
		if MSA[seq_idx][pos] == sym:
			total += weights[seq_idx]
	# from the text this should be:
	# 1/(lam+Meff) * total
	# but that never gives probabilities that sum to 1
	return 1/Meff * total

def norm_pair_freq(sym1,pos1,sym2,pos2,MSA,weights,lam,alphabet_size):
	Meff  = sum(weights)
	total = lam/(alphabet_size**2)
	for seq_idx1 in xrange(len(MSA)):
		for seq_idx2 in xrange(len(MSA)):
			if MSA[seq_idx1][pos2] == sym1 and MSA[seq_idx2][pos2] == sym2:
				total += weights[seq_idx]
	# from the text this should be:
	# 1/(lam+Meff) * total
	# but that never gives probabilities that sum to 1
	return 1/Meff * total

# read in the MSA (matrix A from equation 1 Morcos et al, 2011)
MSA = read_MSA(filename)

# calculate sequence weights for each sequence in the MSA (1/m_a) from in equation 3
weights = []
weights_fn = 'seq.weights'
if os.path.exists(weights_fn):
	weights = pickle.load(open(weights_fn,'r'))
else:
	weights = calc_seq_weights(MSA)
	pickle.dump(weights,open(weights_fn,'w'))

# parameters for the model (two-body = e, one-body = h)
e = [ [ 0 for i in xrange(len(MSA)) ] for j in xrange(len(MSA)) ]
h = [ 0 for i in xrange(len(MSA)) ]

