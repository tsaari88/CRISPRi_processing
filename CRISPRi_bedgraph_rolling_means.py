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


def exclude_sparse_regions(df_bg, df_raw, bgint, introll):
	'''
	Takes a bedgraph df with evenly-spaced intervals, and the df with guides' positional values,
	counts # guides within a bedgraph interval, and rolls over several intervals to get rolling counts.
	Then large regions (spanning a number of intervals) without guides can be filtered out.
	'''
	#Bin and count
	bg_bins = df_bg['start']
	bg_bins = bg_bins.append(pd.Series([bg_bins.max() + bgint])).reset_index(drop=True)
	count = df_raw['position'].groupby(pd.cut(df_raw.position, bg_bins)).count()
	#Roll the counts to consider multiple intervals at a time
	rolling_count = count.rolling(introll, min_periods=1).sum()
	rolling_count.reset_index(inplace=True, drop=True)
	#Add rolling counts to bedgraph and exclude based on these
	df_bg['rolling_guidecount'] = rolling_count
	df_filtered = df_bg.loc[df_bg.rolling_guidecount > 0].copy()
	df_filtered.drop(columns=['rolling_guidecount'], inplace=True)
	return(df_filtered)


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
	df_rawdata['score'] *= -1 #Negate values so that essential areas are positive

	#Calculate rolling means and generate bedgraph from these values
	df_bedgraph = rolling_means_bedgraph(df_rawdata, nguides=int(args.n_guides), bgint=int(args.interval_size))
	df_bedgraph['chr'] = region_chr
	df_bedgraph = df_bedgraph[['chr', 'start', 'end', 'score']]

	#Exclude sparse regions (no guides) from bedgraph
	df_bedgraph = exclude_sparse_regions(df_bedgraph, df_rawdata, bgint=int(args.interval_size), introll=10)

	if args.subtract_background:
		df_bedgraph['score'] = df_bedgraph['score'] - df_bedgraph['score'].median()

	df_bedgraph.to_csv(args.output, sep='\t', header=False, index=False)

