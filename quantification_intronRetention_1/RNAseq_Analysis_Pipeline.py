import os
from subprocess import *

PARENT_DIR='/path/to/intron-retention-paper'
BASE_DIR='{}/quantification_intronRetention_1'.format(PARENT_DIR)

# the directory of storing the results:
# star mapping: TARGET_DIR/{tissue}/{sample}/
# kallisto: TARGET_DIR/{tissue}/{sample}/
# intron retention quantification: TARGET_DIR/{tissue}/{sample}/
TARGET_DIR = '{}/results/'.format(BASE_DIR) # (trailing slash)
# directory of fastq files: DATA_DIR/{tissue}/fastq/{sample}_{1,2}.fastq
DATA_DIR = '{}/data/'.format(BASE_DIR) # (trailing slash)
SIRI_cmd = os.path.dirname(os.path.abspath(__file__)) + '/SIRI/bin/SIRI'
STARmm10 = '{}/star_index/'.format(BASE_DIR) # the path for the STAR index (trailing slash)
KALLISTO_INDEX = '{}/Mus_musculus.GRCm38.cdna.all.fa.idx'.format(BASE_DIR) # the path for the Kallisto index
GTF_FILE = '{}/data/Mus_musculus.GRCm38.91.chr.gtf'.format(PARENT_DIR) # the path for the GTF file
tissue = 'mESC'
sample = 'Chr'

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

def get_bam_length(bam_file): # get read length given bam files
    cmd = '{0} view {1} | head -n 1 | cut -f 10'.format('samtools', bam_file)
    p = Popen(cmd, shell=True, stdout=PIPE).communicate()[0].strip().split('\n')
    return len(p[0])

def run_star(tissue, sample): # run STAR mapping
    os.chdir(TARGET_DIR + tissue + '/' + sample)
    fastq_1 = DATA_DIR + tissue + '/fastq/' + sample + '_1.fastq'
    fastq_2 = DATA_DIR + tissue + '/fastq/' + sample + '_2.fastq'
    cmd = 'STAR --genomeDir {} ' \
          '--runThreadN 8 --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 17179869184 ' \
          '--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 ' \
          '--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 ' \
          '--alignMatesGapMax 1000000 --seedSearchStartLmax 30 --alignEndsType EndToEnd ' \
          '--readFilesIn {} {}'.format(STARmm10, fastq_1, fastq_2)
    call(cmd, shell=True)

def run_kallisto(tissue, sample): # run Kallisto, gene expression quantification
    os.chdir(TARGET_DIR + tissue + '/' + sample)
    fastq_1 = DATA_DIR + tissue + '/fastq/' + sample + '_1.fastq'
    fastq_2 = DATA_DIR + tissue + '/fastq/' + sample + '_2.fastq'
    fastq_files = fastq_1 + ' ' + fastq_2
    makedirs('kallisto')
    os.chdir(TARGET_DIR + tissue + '/' + sample + '/kallisto')
    cmd = 'kallisto quant --index={1} --output-dir=./ --threads 4 -b 100 {0}'.format(fastq_files, KALLISTO_INDEX)
    call(cmd, shell=True)

def run_SIRI(tissue, sample): # run SIRI, intron retention quantification
    os.chdir(TARGET_DIR + tissue + '/' + sample)
    bam_file = TARGET_DIR + tissue + '/' + sample + '/Aligned.sortedByCoord.out.bam'
    out_dir = './siri_out'
    read_length = get_bam_length(bam_file)
    cmd = '{} ' \
          '--bam_files {} ' \
          '--gtf {} ' \
          '--anchor 8 --length {} --lib first --read P ' \
          '--output {}'.format(SIRI_cmd, bam_file, GTF_FILE, read_length,
                               out_dir)
    call(cmd, shell=True)

run_star(tissue, sample)
run_kallisto(tissue, sample)
run_SIRI(tissue, sample)
