import cluster.xmeans as xmeans
import cluster.kmeans as kmeans
from collections import defaultdict

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

X = [] # array to store the FI values from Chromosome, Nucleus, Cytoplasm.
SIRI_ID2information = {} # store the information for each introns by SIRI ID
SIRI_ID_list = []
fp = open('example.txt')
header = fp.readline()
for line in fp:
    sp = line.strip().split('\t')
    X.append([float(sp[2]), float(sp[3]), float(sp[4])]) # FI values in Chromosome, Nucleus, Cytoplasm.
    SIRI_ID2information[sp[1]] = line.strip()
    SIRI_ID_list.append(sp[1])
fp.close()

print '...initial {} clusters'.format(2)
initial_centers = xmeans.kmeans_plusplus_initializer(X, 2).initialize()
xmeans_instance = xmeans.xmeans(X, initial_centers)

print '...run xmeans'
xmeans_instance.process()

print '...get clusters'
clusters = xmeans_instance.get_clusters()

print '...cluster num {}'.format(len(clusters))

cluster_num = len(clusters)
fw = open('results.txt', 'w')
fw.write(header.strip() + '\t' + 'cluster\n')
for index in xrange(0, cluster_num):
    cluster_name = 'Cluster-{}'.format(index + 1)
    for _c in clusters[index]:
        fw.write(SIRI_ID2information[SIRI_ID_list[_c]] + '\t' + cluster_name + '\n')
fw.close()

