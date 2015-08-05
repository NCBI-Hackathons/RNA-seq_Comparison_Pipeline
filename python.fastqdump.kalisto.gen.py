f=open('run.id.txt','r')
f1=f.readlines()
f.close()

f=open('fastqdump.kallisto.sh','w')

for row in f1:
	f.write('fastq-dump --split-3 --fasta '+row[:-1]+'\n')	
	f.write('kallisto quant -i ~/data/refseq.idx -o '+ 'output_'+row[:-1]+' -t 4 '+ row[:-1]+'_1.fasta '+row[:-1]+'_2.fasta \n')

