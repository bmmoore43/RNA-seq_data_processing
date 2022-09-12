import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is to wirte sh and sub files for HiSAT2 and HTseq')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-trimmed_list', help='tab-delimited file with list of trimmed fastq files to input; for PE R1 in first column, R2 in second column', required=True)
	req_group.add_argument('-genome_seq', help='genome sequence fasta file', required=True)
	req_group.add_argument('-gff', help='gff3 file', required=True)
	req_group.add_argument('-layout', help='PE: paired end, SE: single end', required=True)
	req_group.add_argument('-dir', help='the path to your staging directory on chtc', required=True)
	req_group.add_argument('-base', help='base name for genome index', required=True)	
	
	# optional
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-threads', help='how many threads to be used',default=4)
	inp_group.add_argument('-mem', help='mem asked for job running',default=60)
	inp_group.add_argument('-disk', help='disk space asked for job running',default=100)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	
	
	# write input list
	print("writing trimmed fastq list for HiSAT2 and HTSeq")
	trim_file = args.trimmed_list
	trim_file_open= open(trim_file,"r")
	trim_file_out= open("hisat2-htseq_filelist.csv","w")
	for line in trim_file_open:
		L=line.strip().split("\t")
		R1z=L[0]
		R2z=L[1]
		R1= R1z.replace('.gz','')
		R2= R2z.replace('.gz','')
		sam= R1.split("_R1")[0]+"_sam"
		sams= sam+"_sorted.sam"
		htseq= "HTSeqCount_"+sams+".out"
		trim_file_out.write("%s,%s,%s,%s,%s,%s,%s\n" %(R1z,R2z,R1,R2,sam,sams,htseq))
		
	trim_file_open.close()
	trim_file_out.close()
	
	print("write hisat2-htseq.sh file")
	hisat_shfile= open("hisat2-htseq.sh","w")
	hisat_shfile.write("#!/bin/bash\nunzip hisat2-2.2.1-Linux_x86_64.zip\n") #unpack hisat2
	hisat_shfile.write("chmod -R 755 hisat2-2.2.1/\n") #change permissions
	hisat_shfile.write("export PATH\n. /etc/profile.d/modules.sh\nmodule load GCC/8.3.0\n") # load GCC module
	hisat_shfile.write("hisat2-2.2.1/hisat2\n") # test program is working
	hisat_shfile.write("cp %s/$1 ./\n" % args.dir) # copy file from staging
	hisat_shfile.write("cp %s/$2 ./\n" % args.dir)
	hisat_shfile.write("gunzip $1\n") # unzip trimmed files
	hisat_shfile.write("gunzip $2\n")
	# build the index for the genome
	hisat_shfile.write("cp %s/%s ./\n" %(args.dir,args.genome_seq)) #cp genome from staging
	hisat_shfile.write('hisat2-2.2.1/hisat2-build %s %s\n'%(args.genome_seq,args.base))
	# map reads
	if args.layout == 'PE':
		# map the reads to the genome
		hisat_shfile.write("hisat2-2.2.1/hisat2 -x %s -1 $3 -2 $4 -S $5 -q --dta --new-summary -p %s --no-discordant\n" % (args.base,args.threads))
		#rm trim files
		hisat_shfile.write("rm *fastq.paired.trim*\n")
		# unzip samtools and install
		hisat_shfile.write("tar xjf samtools-1.15.1.tar.bz2\n")
		hisat_shfile.write("chmod -R 777 samtools-1.15.1/\n")
		hisat_shfile.write("cd samtools-1.15.1/\n")
		hisat_shfile.write("./configure --disable-lzma\n")
		hisat_shfile.write("make\n")
		hisat_shfile.write("make install\n")
		hisat_shfile.write("cd ../\n")
		# sort the sam file- htseq needs sam file sorted by name (-n)
		hisat_shfile.write("samtools-1.15.1/samtools sort -O sam -n $5 -o $6\n")
		# remove unsorted sam and index
		hisat_shfile.write("rm $5\n")
		hisat_shfile.write("rm %s*\n" % args.base)
		## untar your Python installation. Make sure you are using the right version!
		hisat_shfile.write("tar -xzf python37.tar.gz\n")
		## untar packages
		hisat_shfile.write("tar -xzf packages2.tar.gz\n")
		## export python and packages to working directory
		hisat_shfile.write("export PATH=$PWD/python/bin:$PATH\n")
		hisat_shfile.write("export PYTHONPATH=$PWD/packages2\n")
		hisat_shfile.write("export HOME=$PWD\n")
		# get uniquely mapped reads and read counts
		## run script to get read counts
		hisat_shfile.write("python3 -m HTSeq.scripts.count --format=sam -m union -s yes -t gene -i ID --nonunique=none -n %s $6 %s > $7\n" % (args.gff, args.threads))
		## remove sorted sam
		hisat_shfile.write("rm $6\n")
	elif args.layout == 'SE':
		# map the reads to the genome
		hisat_shfile.write('hisat2 -p 4 --dta -x %s -U $3 -S $5\n'%(args.base))
		# sort the sam file
		hisat_shfile.write('samtools sort -n -O sam, $5  -o $6\n')
		# get read counts
		hisat_shfile.write('htseq-count --format=sam -m union -s yes -t gene -i ID --nonunique=none -n %s $6 %s > $7\n'% (args.gff, args.threads))
	else:
		print("must indicate PE or SE")
	
	hisat_shfile.close()
	print("write hisat2-htseq.sub file")
	hisat_subfile= open("hisat2-htseq.sub","w")
	hisat_subfile.write("universe = vanilla\n")
	hisat_subfile.write("log = hisat2-htseq_$(Cluster).log\n")
	hisat_subfile.write("error = hisat2-htseq_$(Cluster).err\n")
	hisat_subfile.write("output = hisat2-htseq_$(cluster).out\n")
	hisat_subfile.write("executable = hisat2-htseq.sh\n")
	hisat_subfile.write("arguments = $(infilepath1) $(infilepath2) $(infilepath3) $(infilepath4) $(infilepath5) $(infilepath6) $(infilepath7)\n")
	hisat_subfile.write("should_transfer_files = YES\n")
	hisat_subfile.write("when_to_transfer_output = ON_EXIT\n")
	hisat_subfile.write("transfer_input_files = hisat2-2.2.1-Linux_x86_64.zip,samtools-1.15.1.tar.bz2,http://proxy.chtc.wisc.edu/SQUID/chtc/python37.tar.gz,packages2.tar.gz,%s\n" % (args.gff))
	hisat_subfile.write("requirements = (HasCHTCStaging == true) && (HasChtcSoftware == true) && (OpSysMajorVer =?= 7)\n")
	hisat_subfile.write("request_cpus = %s\n" % args.threads)
	hisat_subfile.write("request_memory = %sGB\n" % args.mem)
	hisat_subfile.write("request_disk = %sGB\n" % args.disk)
	hisat_subfile.write("queue infilepath1,infilepath2,infilepath3,infilepath4,infilepath5,infilepath6,infilepath7 from hisat2-htseq_filelist.csv\n")


if __name__ == '__main__':
	main()
	
	
