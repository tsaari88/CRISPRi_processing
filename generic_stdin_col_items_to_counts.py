#!/usr/bin/env python

import sys
from collections import defaultdict

counts = defaultdict(int)

if len(sys.argv) < 2:
	sys.exit('Usage: %s INT.\n   Where INT is index of column number to count (0-based).' % sys.argv[0])

indx = int(sys.argv[1])

for line in sys.stdin:
	linedata = line.split("\t")
	
	item = linedata[indx]
	counts[item] += 1

#Convert to dict - avoid defaultdict's strange behavior when converting
#from int to string (prints to console during conversion, not sure why)
countsdict = dict(counts)

for key, value in countsdict.items():
	count = str(value)
	print(key + '\t' + count)