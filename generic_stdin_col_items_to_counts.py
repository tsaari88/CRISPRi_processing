#!/usr/bin/env python

import sys
from collections import defaultdict

counts = defaultdict(int)

if len(sys.argv) < 2:
	sys.exit('Usage: %s INT.\n   Where INT is index of column number to count (0-based).' % sys.argv[0])

#Which column
indx = int(sys.argv[1])

#Gather counts
for line in sys.stdin:
	linedata = line.split("\t")
	item = linedata[indx]
	counts[item] += 1

#Print counts to stdout
for key, value in counts.items():
	count = str(value)
	print(key + '\t' + count)