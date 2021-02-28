## xmeans cluster of U introns for ESC, NPC and Ctx
import os
import xmeans_cluster

# We require the intron with skipping counts >= 2 in at least one cell compartment (detailed information is available in our method section)
# The next is followed by obtaining the introns that will not be used in our clustering analysis (black list)
xmeans_cluster.generate_black_list()

# run xmeans cluster for ESC
Chr_index, Nuc_index, Cyto_index = xmeans_cluster.POLYA_SAMPLE_INDEX[0]
cell_name = xmeans_cluster.CELL_TYPE[0]
xmeans_cluster.get_intron_cluster(Chr_index, Nuc_index, Cyto_index, cell_name, 'polyA')

# run xmeans cluster for NPC
Chr_index, Nuc_index, Cyto_index = xmeans_cluster.POLYA_SAMPLE_INDEX[1]
cell_name = xmeans_cluster.CELL_TYPE[1]
xmeans_cluster.get_intron_cluster(Chr_index, Nuc_index, Cyto_index, cell_name, 'polyA')

# run xmeans cluster for Ctx\n",
Chr_index, Nuc_index, Cyto_index = xmeans_cluster.POLYA_SAMPLE_INDEX[2]
cell_name = xmeans_cluster.CELL_TYPE[2]
xmeans_cluster.get_intron_cluster(Chr_index, Nuc_index, Cyto_index, cell_name, 'polyA')

os.system('Rscript draw.cluster.R')
