import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for fastqc and trimmomatic to write sh and sub files')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-fastq_list', help='tab-delimited file with list of fastq files to process, for PE R1 in first column, R2 in second column', required=True)
	req_group.add_argument('-layout', help='PE: paired end, SE: single end', required=True)
	req_group.add_argument('-dir', help='the path to your staging directory on chtc', required=True)
	#Optional
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-SRA', help='SRA sample ID', default='NA')
	inp_group.add_argument('-threads', help='how many threads to be used',default=4)
	inp_group.add_argument('-adapters', help='adapter sequences to be trimmed',default='TruSeq3-PE.fa')
	inp_group.add_argument('-seedMis', help='seed mismatches: specifies the maximum mismatch count which will still allow a full match to be performed',default=2)
	inp_group.add_argument('-pClipThres', help='palindromeClipThreshold: specifies how accurate the match between the two adapter ligated reads must be for PE palindrome read alignment.',default=30)
	inp_group.add_argument('-sClipThres', help='simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.',default=10)
	inp_group.add_argument('-LEADING', help='Specifies the minimum quality required to keep a base',default=3)
	inp_group.add_argument('-TRAILING', help='Specifies the minimum quality required to keep a base.',default=3)
	inp_group.add_argument('-windowSize', help='specifies the number of bases to average across',default=4)
	inp_group.add_argument('-requiredQuality', help='specifies the average quality required',default=20)
	inp_group.add_argument('-MINLEN', help='Specifies the minimum length of reads to be kept.',default=36)
	inp_group.add_argument('-HEADCROP', help='Specifies number of reads at beginning to be removed', default=0)
	inp_group.add_argument('-mem', help='mem asked for job running',default=100)
	inp_group.add_argument('-disk', help='disk space asked for job running',default=100)
	
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	
	
	# write fastq_list.csv for fastqc_trimm.sh file
	print("writing fastq list for fastqc and trimmomatic")
	fastq_file = args.fastq_list
	fastq_file_open= open(fastq_file,"r")
	fastq_out= open("fastq_list.csv","w")
	for line in fastq_file_open:
		L=line.strip().split("\t")
		R1z=L[0]
		R2z=L[1]
		R1= R1z.replace('.gz','')
		R2= R2z.replace('.gz','')
		R1p= R1+".paired.trim"
		R1u= R1+".unpaired.trim"
		R2p= R2+".paired.trim"
		R2u= R2+".unpaired.trim"
		fastq_out.write("%s,%s,%s,%s,%s,%s,%s,%s\n"%(R1z,R2z,R1,R2,R1p,R1u,R2p,R2u))
	fastq_file_open.close()
	fastq_out.close()
	
	#write fastqc_trimm.sh file
	print('write fastqc_trimm.sh file')
	fqct_file=open('fastqc_trimm.sh','w')
	fqct_file.write("#!/bin/bash\nunzip fastqc_v0.11.9.zip\nchmod 755 FastQC/fastqc\nFastQC/fastqc\n") #unzip fastqc, change permissions, test run
	if args.layout == 'PE' and args.SRA == 'NA':
		fqct_file.write("cp %s/$1 ./\n" % args.dir) #copy files from staging
		fqct_file.write("cp %s/$2 ./\n" % args.dir)
		fqct_file.write("gunzip $1\n") #unzip files
		fqct_file.write("gunzip $2\n")
		fqct_file.write("FastQC/fastqc -f fastq $3\n") # run fastqc
		fqct_file.write("FastQC/fastqc -f fastq $4\n")
		fqct_file.write("unzip Trimmomatic-0.39.zip\n") #unpack trimmomatic
		fqct_file.write("cp Trimmomatic-0.39/adapters/%s ./\n" % args.adapters) # copy adapters to working directory
		if args.HEADCROP == 0:
			fqct_file.write("java -jar Trimmomatic-0.39/trimmomatic-0.39.jar %s -phred33 -threads %s $3 $4 $5 $6 $7 $8 ILLUMINACLIP:%s:%s:%s:%s:2:True LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s\n" % (args.layout, args.threads, args.adapters, args.seedMis, args.pClipThres, args.sClipThres, args.LEADING,args.TRAILING, args.windowSize, args.requiredQuality, args.MINLEN)) #run trimmomatic
		else:
			fqct_file.write("java -jar Trimmomatic-0.39/trimmomatic-0.39.jar %s -phred33 -threads %s $3 $4 $5 $6 $7 $8 ILLUMINACLIP:%s:%s:%s:%s:2:True LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s HEADCROP:%s\n" % (args.layout, args.threads, args.adapters, args.seedMis, args.pClipThres, args.sClipThres, args.LEADING,args.TRAILING, args.windowSize, args.requiredQuality, args.MINLEN, args.HEADCROP))
		fqct_file.write("FastQC/fastqc -f fastq $5\n") # redo fastqc on new file
		fqct_file.write("FastQC/fastqc -f fastq $7\n")
		fqct_file.write("gzip *trim\n") #gzip files and copy to staging
		fqct_file.write("cp *trim.gz %s/\n" % args.dir)
		fqct_file.write("rm $3\n")# remove big files
		fqct_file.write("rm $4\n")
		fqct_file.write("rm *trim.gz\n")
		# done- QC files will be transferred to home directory
	if args.layout == 'SE' and args.SRA == 'NA':
		print("to add SE")
	# SRAToolkit needed for SRA files
	# to run PE: ('fastq-dump --split-files %s\n'%SRA)
	# to run SE: ('fastq-dump %s\n'%SRA)
		
	#write fastqc_trimm.sub file
	print('write fastqc_trimm_loop.sub file')
	fqct_sub=open('fastqc_trimm_loop.sub','w')
	fqct_sub.write("universe = vanilla\n")
	fqct_sub.write("log = fastqc_$(Cluster).$(Process).log\n")
	fqct_sub.write("error = fastqc_$(Cluster).$(Process).err\n")
	fqct_sub.write("output = fastqc_$(cluster).$(Process).out\n")
	fqct_sub.write("executable = fastqc_trimm.sh\n")
	fqct_sub.write("arguments = $(infilepath1) $(infilepath2) $(infilepath3) $(infilepath4) $(infilepath5) $(infilepath6) $(infilepath7) $(infilepath8)\n")
	fqct_sub.write("should_transfer_files = YES\n")
	fqct_sub.write("when_to_transfer_output = ON_EXIT\n")
	fqct_sub.write("transfer_input_files = fastqc_v0.11.9.zip,Trimmomatic-0.39.zip\n")
	fqct_sub.write("requirements = (Target.HasJava == true) && (Target.HasCHTCStaging == true)\n")
	fqct_sub.write("request_cpus = %s\n" % args.threads)
	fqct_sub.write("request_memory = %sGB\n" % args.mem)
	fqct_sub.write("request_disk = %sGB\n" % args.disk)
	fqct_sub.write("queue infilepath1,infilepath2,infilepath3,infilepath4,infilepath5,infilepath6,infilepath7,infilepath8 from fastq_list.csv\n")
	
	fqct_file.close()
	fqct_sub.close()

if __name__ == '__main__':
	main()
	
	
