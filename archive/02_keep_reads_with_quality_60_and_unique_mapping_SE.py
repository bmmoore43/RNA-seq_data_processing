import sys,os
sam = sys.argv[1]
out = open(sam.split('.sam')[0] + '_quality_60_unique.sam','w')
inp = open(sam,'r')
inl = inp.readlines()
#read = []
for lines in inl:
	lines = lines.strip(' ')
	if lines.startswith('@'):
		out.write(lines)
		print(lines)
	else:
		if lines.split('\t')[4] >= '60':  ### this column shows the Mapping quality of HISAT2. 
			out.write(lines)
out.close()



