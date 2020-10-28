#_*_coding:utf8_*_
import os

bowtie_index = '/gss1/seqlib/tair10/TAIR10_GFF3_genes_transposons.gff' + ' ' + '/gss1/seqlib/tair10/TAIR10'
bowtie_index_fa = '/gss1/seqlib/tair10/TAIR10.fa'
input_fastq_dir = '/gss1/home/gaohao/RNA_seq_process'
fastq_name = ['SRR11308177_1.fastq','SRR11308177_2.fastq','SRR11308178_1.fastq','SRR11308178_2.fastq']
sample_label = ['stressed','unstressed']
current_dir = '/gss1/home/gaohao/RNA_seq_autoprocess' #程序运行所在的目录

gtf_path = open('assemblies.txt','w')
for x in range(1,int(len(fastq_name)/2)+1):
    tophat_command = 'bsub -n 16 -K -o tophat.log tophat -p 8 --solexa-quals -o %s -G %s %s%s%s %s%s%s'%\
                     (fastq_name[2*x-1][:11],bowtie_index,input_fastq_dir,os.sep,fastq_name[2*x-2],input_fastq_dir,os.sep,fastq_name[2*x-1])
    #-K 提交作业，并且等待作业完成。当提交作业后，终端打印“waiting for dispath”。当作业完成后，终端打印“job is finished”。作业没有完成，不能提交新的作业
    os.system(tophat_command)

    cufflinks_command = 'bsub -n 16 -K -o cufflinks.log cufflinks -p 8 -o %s/%s/cufflinks %s%s%s%saccepted_hits.bam'%\
                        (current_dir,fastq_name[2*x-1][:11],current_dir,os.sep,fastq_name[2*x-1][:11],os.sep)
    #print(cufflinks_command)
    os.system(cufflinks_command)

    gtf_path.write('%s%s%s/cufflinks/transcripts.gtf\n'%(current_dir,os.sep,fastq_name[2*x-1][:11]))
gtf_path.close()

cuffmerge_command = 'bsub -n 16 -K -o cuffmerge.log cuffmerge -g %s -s %s -p 8 assemblies.txt'%(bowtie_index.split(' ')[0],bowtie_index_fa)
#print(cuffmerge_command)
os.system('mkdir cuffmerge')
os.system('cp %s/assemblies.txt %s/cuffmerge/'%(current_dir,current_dir))
os.chdir('%s/cuffmerge'%(current_dir)) #切换目录要用os.chdir
os.system(cuffmerge_command)
os.chdir(current_dir)

cuffdiff_command = 'bsub -n 16 -K -o cuffdiff.log cuffdiff -o ./diff_out/ -b %s -p 8 -L %s,%s -u %s/cuffmerge/merged_asm/merged.gtf'%\
                   (bowtie_index_fa,sample_label[0],sample_label[1],current_dir)
for x in range(1,int(len(fastq_name)/2)+1):
    cuffdiff_command = cuffdiff_command + ' ' + '%s/%s/accepted_hits.bam'%(current_dir,fastq_name[2*x-1][:11])
#print(cuffdiff_command)
os.system('mkdir cuffdiff')
os.chdir('%s/cuffdiff'%(current_dir))
os.system(cuffdiff_command)
os.chdir('%s/cuffdiff/diff_out'%(current_dir))
os.system('''awk '{if(($10<-2)&&($11<0.001))print $3"\t"$8"\t"$9"\t"$10}' gene_exp.diff | grep -v 'inf' > down.txt''')
os.system('''awk '{if(($10>2)&&($11<0.001))print $3"\t"$8"\t"$9"\t"$10}' gene_exp.diff | grep -v 'inf' > up.txt''')



