import cluster.xmeans as xmeans
import cluster.kmeans as kmeans
from collections import defaultdict
import numpy as np
import os

chromIndexWithOutYChrom = ['chr' + str(x) for x in range(1, 20)] + ['chrX']
CELL_TYPE = ['ESC', 'NPC', 'Ctx']
## We pooled all samples together and the information is available in intron_FI_combined_results.txt. 
## The poly RNA-seq and total RNA-seq are indexed below. 
POLYA_SAMPLE_INDEX = [[[0, 1, 2], [3, 4, 5], [6, 7, 8]], [[18, 19, 20], [21, 22, 23], [24, 25, 26]], [[36, 37, 38], [39, 40, 41], [42, 43, 44]]]
TOTAL_SAMPLE_INDEX = [[[9, 10, 11], [12, 13, 14], [15, 16, 17]], [[27, 28, 29], [30, 31, 32], [33, 34, 35]], [[45, 46, 47], [48, 49, 50], [51, 52, 53]]]
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'
DATA_DIR = os.path.join(CURRENT_DIR, '..', 'data/')

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

def get_counts(counts, index):
    count = 0
    for x in index:
        if counts[x] == 'NA':
            continue
        else:
            count += int(counts[x])
    return count

def get_FI_value(inc_count, skp_count, inc_l, skp_l):
    _1 = float(inc_count) / float(inc_l)
    _2 = float(skp_count) / float(skp_l)
    if _1 + _2 == 0:
        return 'NA'
    else:
        return _1 / (_1 + _2)

def get_l(_list):
    for _l in _list:
        if _l == 'NA':
            continue
        return _l
    return 'NA'

## we require the intron with skipping counts >= 2 in at least one cell compartment (detailed information is available in our method section)
def generate_black_list(): 
    intron_info_list = read_file(DATA_DIR + 'intron_FI_combined_results.txt')
    black_list = []
    for info in intron_info_list[1:]:
        sp = info.split('\t')
        if sp[7] != 'U':
            continue
        sq = sp[0].split('_')
        length = abs(int(sq[2]) - int(sq[1]) + 1)
        if length < 60:
            continue
        inc_counts = sp[8].split(',')
        skp_counts = sp[9].split(',')
        black_flag = True
        for index, sample_index in enumerate(POLYA_SAMPLE_INDEX):
            Chr_index, Nuc_index, Cyto_index = sample_index
            inc_count_chr, skp_count_chr = get_counts(inc_counts, Chr_index), get_counts(skp_counts, Chr_index)
            inc_count_nuc, skp_count_nuc = get_counts(inc_counts, Nuc_index), get_counts(skp_counts, Nuc_index)
            inc_count_cyto, skp_count_cyto = get_counts(inc_counts, Cyto_index), get_counts(skp_counts, Cyto_index)
            if skp_count_chr >= 2 or skp_count_nuc >= 2 or skp_count_cyto >= 2:
                black_flag = False
        if black_flag:
            black_list.append(sp[0])
    fw = open(CURRENT_DIR + 'black_list_intron.txt', 'w')
    for intron_id in black_list:
        fw.write(intron_id + '\n')
    fw.close()

def get_intron_cluster(Chr_index, Nuc_index, Cyto_index, cell_type, experiment, READ_THRESHOLD = 20):
    print('parsing cell type {}, experiment {}'.format(cell_type, experiment))
    intron_info_list = read_file(DATA_DIR + 'intron_FI_combined_results.txt')
    black_list_intron = set(read_file(CURRENT_DIR + 'black_list_intron.txt'))
    X = []
    intron_id = []
    intron_id2gene = {}
    for info in intron_info_list[1:]:
        sp = info.split('\t')
        ## we only keep U introns with length no less than 60 and not in black list.
        if sp[7] != 'U':
            continue
        sq = sp[0].split('_')
        length = abs(int(sq[2]) - int(sq[1]) + 1)
        if length < 60:
            continue
        if sp[0] in black_list_intron:
            continue
        inc_counts = sp[8].split(',')
        skp_counts = sp[9].split(',')
        inc_l = get_l(sp[10].split(','))
        skp_l = get_l(sp[11].split(','))
        inc_count_chr, skp_count_chr = get_counts(inc_counts, Chr_index), get_counts(skp_counts, Chr_index)
        inc_count_nuc, skp_count_nuc = get_counts(inc_counts, Nuc_index), get_counts(skp_counts, Nuc_index)
        inc_count_cyto, skp_count_cyto = get_counts(inc_counts, Cyto_index), get_counts(skp_counts, Cyto_index)
        if inc_count_chr + skp_count_chr < READ_THRESHOLD or inc_count_nuc + skp_count_nuc < READ_THRESHOLD or inc_count_cyto + skp_count_cyto < READ_THRESHOLD:
            continue
        FI_chr = get_FI_value(inc_count_chr, skp_count_chr, inc_l, skp_l)
        FI_nuc = get_FI_value(inc_count_nuc, skp_count_nuc, inc_l, skp_l)
        FI_cyto = get_FI_value(inc_count_cyto, skp_count_cyto, inc_l, skp_l)
        X.append([FI_chr, FI_nuc, FI_cyto])
        intron_id.append(sp[0])
        intron_id2gene[sp[0]] = sp[1]
    X = np.array(X)

    print('...initial {} clusters'.format(2))
    initial_centers = xmeans.kmeans_plusplus_initializer(X, 2).initialize()
    xmeans_instance = xmeans.xmeans(X, initial_centers)
    print('...run xmeans')
    xmeans_instance.process()
    print('...get clusters')
    clusters = xmeans_instance.get_clusters()
    print('...cluster num {}'.format(len(clusters)))

    print('...saving to current directory')
    fw = open(CURRENT_DIR + '{}.cluster.{}.txt'.format(cell_type, experiment), 'w')
    fw.write('Cell_Compartment\tCluster\tFI\tSIRI_ID\n')
    values = [np.mean(X[clusters[x]][:,0]) for x in range(len(clusters))]
    cluster_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    cluster_name_list = []
    for index in np.argsort(values):
        cluster = clusters[index]
        for _c in cluster:
            fw.write('Chr\t{}\t{}\t{}\n'.format(cluster_name[index] + ': ' + str(len(cluster)), X[_c][0], intron_id[_c]))
            fw.write('Nuc\t{}\t{}\t{}\n'.format(cluster_name[index] + ': ' + str(len(cluster)), X[_c][1], intron_id[_c]))
            fw.write('Cyto\t{}\t{}\t{}\n'.format(cluster_name[index] + ': ' + str(len(cluster)), X[_c][2], intron_id[_c]))
        print('...{}'.format(cluster_name[index] + ' ' + str(len(cluster))))
        cluster_name_list.append(cluster_name[index] + '\\n ' + str(len(cluster)))
    fw.close()

def main():
    generate_black_list()
    for index, sample_index in enumerate(POLYA_SAMPLE_INDEX):
        Chr_index, Nuc_index, Cyto_index = sample_index
        # get intron clusters by using x-means
        get_intron_cluster(Chr_index, Nuc_index, Cyto_index, CELL_TYPE[index], 'polyA')

if __name__=="__main__":
    main()
