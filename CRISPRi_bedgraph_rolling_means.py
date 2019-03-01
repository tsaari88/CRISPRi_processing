#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy import interpolate


def nearest_nice(x, nice):
	'''
	Round to the closest specified number
	'''
	quotient = x//nice
	remainder = x%nice
	if np.round(remainder/nice) > 0:
		quotient = quotient + 1
	return(nice*quotient)
 

def rolling_means_bedgraph(df, nguides, bgint):
	'''
	Calculate rolling means for N-guide windows,
	then use scipy to interpolate from these onto evenly-spaced X's,
	and generate a bedgraph interval (without chr info)
	'''
	df.sort_values(by='position', inplace=True)
	#Position median allows for consecutive evaluated window positions to roll as intended
	nguide_pos = df['position'].rolling(nguides).median().dropna()
	nguide_score = df_rawdata['score'].rolling(nguides).mean().dropna()
	#
	fxn = interpolate.interp1d(nguide_pos, nguide_score, bounds_error=False)
	#
	minx = nearest_nice(df['position'].min(), bgint)
	maxx = nearest_nice(df['position'].max(), bgint)
	#
	bg_starts = pd.Series(range(minx, maxx, bgint), name='start')
	bg_ends = bg_starts + bgint
	bg_ends.rename('end', inplace=True)
	range_center = ( bg_starts + (bgint/2) )
	#
	interpol_y = pd.Series(fxn(range_center), name='score')
	return(pd.concat([bg_starts, bg_ends, interpol_y], axis=1))


if __name__ == '__main__':

	parser = argparse.ArgumentParser(prog='python CRISPRi_bedgraph_rolling_means.py', description = "Use rolling means to estimate regional signal, and generate a valid bedgraph from this information.")
	parser.add_argument('-b', '--bed_input', required=True, help="6-column bed file with genomic coordinates and OligoIDs from one tiling region. Do not combine tiling regions")
	parser.add_argument('-r', '--results_input', required=True, help="Results file containing OligoID and mean_log2FC columns")
	parser.add_argument('-o', '--output', required=True, help="Output bedgraph of smoothed signal")
	parser.add_argument('-g', '--n_guides', required=True, help="Number of guides per window")
	parser.add_argument('-i', '--interval_size', default=25, help="Interval size to use in output bedgraph")
	parser.add_argument('-s', '--subtract_background', default=False, action='store_true', help="Boolean flag - if set, subtract median background signal")

	args = parser.parse_args()

	#prepare bed_input
	bed_6col = pd.read_csv(args.bed_input, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand"])
	bed_chr = bed_6col['chr'].unique()
	if not len(bed_chr) == 1:
		raise ValueError("Bed file appears to have multiple chromosomes. This shouldn't be the case; it should contain a single tiling region")
	else:
		region_chr = "".join(bed_chr)
	#prepare results_input
	expt = pd.read_csv(args.results_input, sep="\t")
	if (not 'OligoID' in expt.columns.tolist()) or (not 'mean_log2FC' in expt.columns.tolist()):
		raise ValueError("Results file must contain columns labeled OligoID and mean_log2FC")

	#create df of position and score of each guide
	bed_6col['position'] = pd.to_numeric(( (bed_6col['start'] + bed_6col['end']) / 2 ).round(), downcast='integer')
	df_position = bed_6col[['name', 'position']]
	df_scores = expt[['OligoID', 'mean_log2FC']]
	df_rawdata = df_position.merge(df_scores, how='left', left_on='name', right_on='OligoID').dropna().drop(['name', 'OligoID'], 1)
	df_rawdata.rename({'mean_log2FC' : 'score'}, axis=1, inplace=True)

	#Calculate rolling means and generate bedgraph from these values
	df_bedgraph = rolling_means_bedgraph(df_rawdata, nguides=int(args.n_guides), bgint=int(args.interval_size))
	df_bedgraph['chr'] = region_chr
	df_bedgraph = df_bedgraph[['chr', 'start', 'end', 'score']]

	if args.subtract_background:
		df_bedgraph['score'] = df_bedgraph['score'] - df_bedgraph['score'].median()

	df_bedgraph.to_csv(args.output, sep='\t', header=False, index=False)

