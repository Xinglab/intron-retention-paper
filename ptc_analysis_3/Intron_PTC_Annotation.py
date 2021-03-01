from pysam import FastaFile
import re
from tqdm import tqdm
import string
import os
from collections import defaultdict

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'
DATA_DIR = os.path.join(CURRENT_DIR, '..', 'data/')

MM10_FASTQ = DATA_DIR + 'mm10.fa'
MOUSE_GTF = DATA_DIR + 'Mus_musculus.GRCm38.91.chr.gtf'
if not os.path.exists(MOUSE_GTF):
    print 'please unzip Mus_musculus.GRCm38.91.chr.gtf.gz first'
    sys.exit(0)
REVERSE_COMPLIMENT_TAB = string.maketrans("ACTG", "TGAC")

def fetch_seq(seq_fetch_obj, chrom, start, end, strand):
    seq = seq_fetch_obj.fetch(chrom, start, end).upper()
    if strand == '-':
        seq = seq.translate(REVERSE_COMPLIMENT_TAB)[::-1]
    return seq

def get_transcript_cds():
    fp = open(MOUSE_GTF)
    dict_transcript_cds = {}
    for index, line in enumerate(fp):
        if line.startswith('#'):
            continue
        sp = line.strip().split('\t')
        if sp[2] == 'gene':
            continue
        if sp[2] == 'transcript':
            transcript_id = re.sub('.*transcript_id "|\".*', '', sp[8])
            dict_transcript_cds[transcript_id] = {'CDS':[], 'start_codon':[], 'stop_codon':[], 'strand':sp[6], 'chrom':sp[0], 'exon':[]}
        if sp[2] == 'CDS':
            dict_transcript_cds[transcript_id]['CDS'].append((int(sp[3]) - 1, int(sp[4])))
        if sp[2] == 'start_codon':
            dict_transcript_cds[transcript_id]['start_codon'] = [int(sp[3]) - 1, int(sp[4])]
        if sp[2] == 'stop_codon':
            dict_transcript_cds[transcript_id]['stop_codon'] = [int(sp[3]) - 1, int(sp[4])]
        if sp[2] == 'exon':
            dict_transcript_cds[transcript_id]['exon'].append([int(sp[3]) - 1, int(sp[4])])
    fp.close()
    dict_transcript_cds_has_start_stop_codon = {}
    for transcript_id in dict_transcript_cds:
        ## require the CDS has clear start codon and stop codon annotation
        if len(dict_transcript_cds[transcript_id]['start_codon']) > 1 and len(dict_transcript_cds[transcript_id]['stop_codon']) > 1:
            dict_transcript_cds_has_start_stop_codon[transcript_id] = dict_transcript_cds[transcript_id]
    return dict_transcript_cds_has_start_stop_codon

def check_stop_codon(seq):
    seq = re.findall('...', seq)
    if 'TAG' in seq or 'TAA' in seq or 'TGA' in seq:
        return True
    return False

def get_stop_codon_pos(seq):
    seq = re.findall('...', seq)
    if 'TAG' in seq:
        return seq.index('TAG'), 'TAG'
    if 'TAA' in seq:
        return seq.index('TAA'), 'TAA'
    if 'TGA' in seq:
        return seq.index('TGA'), 'TGA'

def get_intron_PTC():
    ## Our annotation is based on whole transcript reading frame.
    ## We investigate the effect of intron per transcript. 
    U_Intron_List = []
    # get U introns from clustering results in step 2
    for cell in ['ESC', 'NPC', 'Ctx']:
        fp = open(CURRENT_DIR + '../xmeans_cluster_U_intron_2/{}.cluster.polyA.txt'.format(cell))
        fp.readline()
        for line in fp:
            sp = line.strip().split('\t')
            U_Intron_List.append(sp[-1])
        fp.close()
    U_Intron_List = list(set(U_Intron_List))

    CHROM_INDEX_LIST = ['chr' + str(x) for x in xrange(1, 20)] + ['chrX']
    seq_fetch_obj = FastaFile(MM10_FASTQ)
    dict_transcript_cds = get_transcript_cds()
    intron_id2transcripts = {}
    fp = open(DATA_DIR + 'Intron_transcript.txt')
    for line in fp:
        sp = line.strip().split('\t')
        intron_id2transcripts[sp[0]] = sp[1:]
    fp.close()

    fw = open(CURRENT_DIR + 'Intron_PTC_results.txt', 'w')
    fw.write('SIRI_ID\tTranscript_ID\tPTC\n')
    for intron_id in tqdm(U_Intron_List, total = len(U_Intron_List)):
        transcript_list = intron_id2transcripts[intron_id]
        sq = intron_id.split('_')
        chrom, start, end, strand = sq[0], int(sq[1]) -1, int(sq[2]), sq[3]
        if chrom not in CHROM_INDEX_LIST:
            continue
        for transcript in transcript_list:
            if transcript not in dict_transcript_cds:
                continue
            cds_exon_list = dict_transcript_cds[transcript]['CDS']
            exon_list = dict_transcript_cds[transcript]['exon']
            stop_codon_pos = dict_transcript_cds[transcript]['stop_codon']
            if strand == '+':
                last_intron = chrom + '_' + '{}_{}_{}'.format(exon_list[-2][1] + 1, exon_list[-1][0], strand)
                if len(exon_list) >= 3:
                    second_last_intron = chrom + '_' + '{}_{}_{}'.format(exon_list[-3][1] + 1, exon_list[-2][0], strand)
                else:
                    second_last_intron == ''
                if stop_codon_pos[1] <= start:
                    continue
                if cds_exon_list[0][0] >= end:
                    continue
            else:
                last_intron = chrom + '_' + '{}_{}_{}'.format(exon_list[-1][1] + 1, exon_list[-2][0], strand)
                if len(exon_list) >= 3:
                    second_last_intron = chrom + '_' + '{}_{}_{}'.format(exon_list[-2][1] + 1, exon_list[-3][0], strand)
                else:
                    second_last_intron = ''
                if stop_codon_pos[0] >= end:
                    continue
                if cds_exon_list[0][1] <= start:
                    continue
            cds_seq = ''
            cds_intron_seq = ''
            if strand == '+':
                for cds in cds_exon_list:
                    seq = fetch_seq(seq_fetch_obj, chrom, cds[0], cds[1], strand)
                    cds_seq += seq
                    cds_intron_seq += seq
                    if cds[1] == start:
                        cds_intron_seq += fetch_seq(seq_fetch_obj, chrom, start, end, strand)
            if strand == '-':
                for cds in cds_exon_list:
                    seq = fetch_seq(seq_fetch_obj, chrom, cds[0], cds[1], strand)
                    cds_seq += seq
                    cds_intron_seq += seq
                    if cds[0] == end:
                        cds_intron_seq += fetch_seq(seq_fetch_obj, chrom, start, end, strand)
            if check_stop_codon(cds_seq):
                continue
            ## if intron is last intron, it would not cause NMD. We skip intron if it is the last intron of the transcript.
            if intron_id == last_intron:
                continue
            ptc_flag = False
            if check_stop_codon(cds_intron_seq):
                intron_seq = fetch_seq(seq_fetch_obj, chrom, start, end, strand)
                pos, stop_codon = get_stop_codon_pos(cds_intron_seq)
                intron_junc_dis = cds_intron_seq.index(intron_seq) + end - start - pos * 3 - 3
                if intron_junc_dis < 0:
                    continue
                if second_last_intron == intron_id: ## 50nt PTC rule
                    if intron_junc_dis + exon_list[-2][1] - exon_list[-2][0] >= 50:
                        ptc_flag = True
                else:
                    ptc_flag = True
            if ptc_flag:
                fw.write('{}\t{}\t{}\n'.format(intron_id, transcript, 'True'))
            else:
                fw.write('{}\t{}\t{}\n'.format(intron_id, transcript, 'False'))
    fw.close()

def analyze_intron_cluster_ptc():
    dict_ptc_info = defaultdict(lambda:[])
    fp = open(CURRENT_DIR + 'Intron_PTC_results.txt')
    fp.readline()
    for line in fp:
        sp = line.strip().split('\t')
        dict_ptc_info[sp[0]].append((sp[1], sp[2]))
    fp.close()

    fw = open(CURRENT_DIR + 'ptc_response_xmeans.txt', 'w')
    fw.write('cell\tcategory\tproportion\ttype\n')
    for c_index, cell_name in enumerate(['ESC', 'NPC', 'Ctx']):
        fp = open(CURRENT_DIR + '../xmeans_cluster_U_intron_2/{}.cluster.polyA.txt'.format(cell_name))
        fp.readline()
        cluster = {}
        for line in fp:
            sp = line.strip().split('\t')
            cluster[sp[-1]] = sp[1].split(':')[0]
        fp.close()
        cluster_name_list = set(cluster.values())
        for cluster_name in sorted(cluster_name_list):
            ptc_count, no_ptc_count = 0, 0
            for intron in cluster:
                if cluster[intron] != cluster_name:
                    continue
                if intron not in dict_ptc_info:
                    continue
                ptc_flag = False
                for ptc_info in dict_ptc_info[intron]:
                    txID, ptc = ptc_info[0], ptc_info[1]
                    if ptc == 'True':
                        ptc_flag = True
                if ptc_flag:
                    ptc_count += 1
                else:
                    no_ptc_count += 1
            #print cell_name, cluster_name, ptc_count, no_ptc_count, ptc_count / (ptc_count + no_ptc_count + 0.0)
            fw.write('{}\t{}\t{}\t{}\n'.format(cell_name, cluster_name, ptc_count / (ptc_count + no_ptc_count + 0.0), '+ PTC'))
            fw.write('{}\t{}\t{}\t{}\n'.format(cell_name, cluster_name, no_ptc_count / (ptc_count + no_ptc_count + 0.0), '- PTC'))
    fw.close()

def main():
    get_intron_PTC()
    analyze_intron_cluster_ptc()

if __name__=="__main__":
    main()