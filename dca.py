#!/usr/bin/python

import os
import re
import sys
import math
import pickle

# free parameters
filename = 'data/T0594.fasta.chkali-5.sto'
alphabet = '-ARNDCEQGHILKMFPSTWYV'
id_cutoff = 0.8
# lambda value for pseudo-counting, optimally set to Meff in paper
pseudo_lambda = 443

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

nf_cache = {}
def norm_freq(sym,pos,MSA,weights,lam,alphabet_size):
	if not sym in nf_cache:
		nf_cache[sym] = {}
	if not pos in nf_cache[sym]:
		Meff  = sum(weights)
		total = lam/alphabet_size

		for seq_idx in xrange(len(MSA)):
			if MSA[seq_idx][pos] == sym:
				total += weights[seq_idx]
		nf_cache[sym][pos] = 1/(lam+Meff) * total

	return nf_cache[sym][pos]

def norm_pair_freq(sym1,pos1,sym2,pos2,MSA,weights,lam,alphabet_size):
	Meff  = sum(weights)
	total = lam/(alphabet_size**2)
	for seq_idx in xrange(len(MSA)):
		if MSA[seq_idx][pos2] == sym1 and MSA[seq_idx][pos2] == sym2:
			total += 2*weights[seq_idx]
	return 1/(lam+Meff) * total

def calc_mi(pos1,pos2,MSA,weights,lam):
	mi = 0
	for sym1 in alphabet:
		for sym2 in alphabet:
			f1 = norm_freq(sym1,pos1,MSA,weights,lam,len(alphabet))
			f2 = norm_freq(sym2,pos2,MSA,weights,lam,len(alphabet))
			f12_obs = norm_pair_freq(sym1,pos1,sym2,pos2,MSA,weights,lam,len(alphabet))
				#print '(%s,%d,%s,%d,%f,%f,%f)' % (sym1,pos1,sym2,pos2,f1,f2,f12_obs)
			mi += f12_obs * math.log(f12_obs/(f1*f2))
	return mi

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

tot = 0
for sym in alphabet:
	f = norm_freq(sym,1,MSA,weights,pseudo_lambda,len(alphabet))
	tot += f
	print 'f(%s) = %f' % (sym,f)
print 'tot = %f' % tot

mi_values = [ [0 for resi in xrange(len(MSA)) ] for resj in xrange(len(MSA)) ]
n_cols = len(MSA[0])
#n_cols = 10
for resi in xrange(n_cols):
	for resj in xrange(resi+1,n_cols):
		mi_values[resi][resj] = calc_mi(resi,resj,MSA,weights,pseudo_lambda)
		print '%d %d %f' % (resi,resj,mi_values[resi][resj])

# calculate the invertible matrix C
# C[i][j][A][B] = norm_pair_freq(A,i,B,j,...) - (norm_freq(A,i) * norm_freq(B,j))
# it's important to note that C is a (q-1)*L x (q-1)*L matrix, so the dual
# indices (A,i) and (B,j) must be combined into a single index.
# eg (for DNA), the matrix C could look like this:
