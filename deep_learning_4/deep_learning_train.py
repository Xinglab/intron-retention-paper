from keras.models import model_from_json
from keras import optimizers
from keras.layers import Input
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.models import Model
from keras import backend as K
from keras.layers.normalization import BatchNormalization
import numpy as np
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error, roc_curve
from pysam import FastaFile
import string
import os
import sys
from collections import defaultdict
from collections import Counter
import h5py
import random

REVERSE_COMPLIMENT_TAB = string.maketrans("ACTG", "TGAC")
chromIndexWithOutYChrom = ['chr' + str(x) for x in range(1, 20)] + ['chrX']
## please unzip mm10.fa.gz first
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'
DATA_DIR = os.path.join(CURRENT_DIR, '..', 'data/')
if not os.path.exists(os.path.join(DATA_DIR, 'mm10.fa')):
    print 'please unzip mm10.fa.gz first'
    sys.exit(0)

CLUSTER_DIR = os.path.join(CURRENT_DIR, '..', 'xmeans_cluster_U_intron_2/')
MM10SEQ_FASTA_OBJ = FastaFile(os.path.join(DATA_DIR, 'mm10.fa')) 

def seq_fetch(chrom, start, end, strand):
    start, end = int(start), int(end)
    seq = MM10SEQ_FASTA_OBJ.fetch(chrom, start, end).upper()
    if strand == '-':
        seq = seq.translate(REVERSE_COMPLIMENT_TAB)[::-1]
    return seq

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def split_by_fold(instances, fold_id, folds_num):
    samples_per_fold = len(instances) / folds_num
    folds = []
    for x in xrange(0, folds_num - 1):
        folds.append(instances[samples_per_fold * x:samples_per_fold * (x + 1)])
    folds.append(instances[samples_per_fold * (folds_num - 1):])
    test_list = folds[fold_id - 1]
    train_list = []
    for x in xrange(0, folds_num):
        if x == fold_id - 1:
            continue
        train_list.extend(folds[x])
    return set(train_list), set(test_list)

def masked_cross_entropy_loss(y_true, y_pred):
    mask = K.cast(K.not_equal(y_true, -1), K.floatx())
    return K.binary_crossentropy(y_true * mask, y_pred * mask)

def create_model(out_dim, FEATURE_NUM):
    layer1 = 1000
    layer2 = 600
    layer3 = 200
    drop1 = 0.5
    drop2 = 0.2
    drop3 = 0.2
    activation = 'relu'
    initi_function = 'lecun_uniform'
    input_seq = Input(shape=(FEATURE_NUM,))
    model = Dense(layer1, kernel_initializer=initi_function)(input_seq)
    model = BatchNormalization()(model)
    model = Activation(activation)(model)
    model = Dropout(drop1)(model)
    model = Dense(layer2, kernel_initializer=initi_function)(model)
    model = BatchNormalization()(model)
    model = Activation(activation)(model)
    model = Dropout(drop2)(model)
    model = Dense(layer3, kernel_initializer=initi_function)(model)
    model = BatchNormalization()(model)
    model = Activation(activation)(model)
    model = Dropout(drop3)(model)
    model = Dense(out_dim)(model)
    model = Activation('sigmoid')(model)
    model = Model(inputs = input_seq, outputs = model)
    return model

def make_train_data():
    dict_intron_id_value = defaultdict(lambda:[-1] * 12)
    for c_index, cell in enumerate(['ESC', 'NPC', 'Ctx']):
        ## cluster files are obtained from step 2 (xmeans_cluster_U_intron)
        fp = open(CLUSTER_DIR + '{}.cluster.polyA.txt'.format(cell))
        fp.readline()
        for line in fp:
            sp = line.strip().split('\t')
            intron_id = sp[-1]
            intron_split = intron_id.split('_')
            chrom, start, end, strand = intron_split[0], int(intron_split[1]) - 1, int(intron_split[2]), intron_split[3]
            if chrom not in chromIndexWithOutYChrom:
                continue
            _left_seq = seq_fetch(chrom, start, start + 2, strand)
            _right_seq = seq_fetch(chrom, end - 2, end, strand)
            if strand == '+':
                canonical_seq = _left_seq + '-' + _right_seq
            if strand == '-':
                canonical_seq = _right_seq + '-' + _left_seq
            if 'GT-AG' != canonical_seq: ## remove the introns without 'GT-AG' combination
                continue
            cluster = sp[1].split(':')[0]
            if cluster == 'A':
                dict_intron_id_value[sp[-1]][c_index * 4] = 1
                dict_intron_id_value[sp[-1]][c_index * 4 + 1] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 2] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 3] = 0
            elif cluster == 'B':
                dict_intron_id_value[sp[-1]][c_index * 4] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 1] = 1
                dict_intron_id_value[sp[-1]][c_index * 4 + 2] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 3] = 0
            elif cluster == 'C':
                dict_intron_id_value[sp[-1]][c_index * 4] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 1] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 2] = 1
                dict_intron_id_value[sp[-1]][c_index * 4 + 3] = 0
            elif cluster == 'D':
                dict_intron_id_value[sp[-1]][c_index * 4] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 1] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 2] = 0
                dict_intron_id_value[sp[-1]][c_index * 4 + 3] = 1

    intron_id_filtered_list = []
    for intron in dict_intron_id_value:
        y = dict_intron_id_value[intron]
        minus_1_count = y.count(-1)
        # at least one cell has the clustering result
        if minus_1_count > 8:
            continue
        intron_id_filtered_list.append(intron)
    intron_id_filtered_list = set(intron_id_filtered_list)

    intron_id_list = []
    # get feature annotation files
    fp = open(DATA_DIR + 'Intron.feature.annotation')
    header = fp.readline().strip()
    feature_start_index = header.split('\t').index('LogLen.E1')
    dict_intron_feature = {}
    for line in fp:
        sp = line.strip().split('\t')
        if sp[0] not in intron_id_filtered_list:
            continue
        intron_id_list.append(sp[0])
        dict_intron_feature[sp[0]] = sp[feature_start_index:]
    fp.close()

    feature_matrix = []
    for intron_id in intron_id_list:
        feature_matrix.append(dict_intron_feature[intron_id])
    feature_matrix = np.array(feature_matrix, dtype='float32')
    max_feature = np.array([max(1e-100, x) for x in np.max(np.abs(feature_matrix), axis=0)])
    feature_matrix /= max_feature

    intron_value_array = np.array([dict_intron_id_value[intron_id] for intron_id in intron_id_list])

    h5_data = h5py.File(CURRENT_DIR + 'train.data.hd5', 'w')
    h5_data.create_dataset('feature_matrix', data = feature_matrix)
    h5_data.create_dataset('intron_id_list', data = intron_id_list)
    h5_data.create_dataset('intron_value', data = intron_value_array)
    h5_data.close()

def split_train_data(random_seed, fold_id):
    random.seed(random_seed)
    h5_data = h5py.File(CURRENT_DIR + 'train.data.hd5', 'r')
    feature_matrix = h5_data.get('feature_matrix')[:]
    intron_id_list = h5_data.get('intron_id_list')[:]
    intron_value_list = h5_data.get('intron_value')[:]
    intorn_id2index = {intron_id:idx for idx, intron_id in enumerate(intron_id_list)}
    [random.shuffle(intron_id_list) for x in xrange(100)]
    train_list, test_list = split_by_fold(intron_id_list, fold_id, 5)
    X_train, y_train = [], []
    X_test, y_test = [], []

    for intron in train_list:
        X_train.append(feature_matrix[intorn_id2index[intron]])
        y_train.append(intron_value_list[intorn_id2index[intron]])

    for intron in test_list:
        X_test.append(feature_matrix[intorn_id2index[intron]])
        y_test.append(intron_value_list[intorn_id2index[intron]])

    X_train, y_train = np.array(X_train), np.array(y_train)
    X_test, y_test = np.array(X_test), np.array(y_test)

    train_data_dir = CURRENT_DIR + 'train_data/'
    makedirs(train_data_dir)
    np.save(train_data_dir + 'X_test.{}.npy'.format(fold_id), X_test)
    np.save(train_data_dir + 'y_test.{}.npy'.format(fold_id), y_test)
    np.save(train_data_dir + 'X_train.{}.npy'.format(fold_id), X_train)
    np.save(train_data_dir + 'y_train.{}.npy'.format(fold_id), y_train)

def deep_learning_train(random_seed, fold_id):
    np.random.seed(random_seed)
    train_data_dir = CURRENT_DIR + 'train_data/'
    X_train = np.load(train_data_dir + 'X_train.{}.npy'.format(fold_id))
    y_train = np.load(train_data_dir + 'y_train.{}.npy'.format(fold_id))
    FEATURE_NUM = X_train.shape[1]
    sgd = optimizers.SGD(lr=0.01, decay=1e-7, momentum=0.9, nesterov=True)
    model_dir = CURRENT_DIR + 'model/'
    makedirs(model_dir)
    model = create_model(12, FEATURE_NUM)
    model.compile(loss=masked_cross_entropy_loss, optimizer=sgd)
    model.fit(X_train, y_train, batch_size=2000, epochs=200, shuffle=True)
    model_json_file = open(model_dir + 'model.{}'.format(fold_id) + '.json', 'w')
    model_json_file.write(model.to_json())
    model_json_file.close()
    model.save_weights(model_dir + 'model.{}'.format(fold_id) + '.weight.hd5')

def deep_learning_pred(fold_id):
    train_data_dir = CURRENT_DIR + 'train_data/'
    model_dir = CURRENT_DIR + 'model/'
    X_test = np.load(train_data_dir + 'X_test.{}.npy'.format(fold_id))
    y_test = np.load(train_data_dir + 'y_test.{}.npy'.format(fold_id))
    model = model_from_json(open(model_dir + 'model.{}'.format(fold_id) + '.json').read())
    model.load_weights(model_dir + 'model.{}'.format(fold_id) + '.weight.hd5')
    y_pred = model.predict(X_test)
    np.save(train_data_dir + 'y_pred.{}.npy'.format(fold_id), y_pred)

def get_auc_calculation(cluster_idx1, cluster_idx2):
    y_true, y_scores = [], []
    train_data_dir = CURRENT_DIR + 'train_data/'
    for fold_id in range(1, 6):
        y_test = np.load(train_data_dir + 'y_test.{}.npy'.format(fold_id))
        y_pred = np.load(train_data_dir + 'y_pred.{}.npy'.format(fold_id))
        for _pred, _test in zip(y_pred, y_test):
            for cell_index in range(3):
                if _test[4 * cell_index] == -1:
                    continue
                if _test[4 * cell_index + cluster_idx2] == 1:
                    y_true.append(1)
                    y_scores.append(_pred[4 * cell_index + cluster_idx2] / (_pred[4 * cell_index + cluster_idx2] + _pred[4 * cell_index + cluster_idx1]))
                elif _test[4 * cell_index + cluster_idx1] == 1:
                    y_true.append(0)
                    y_scores.append(_pred[4 * cell_index + cluster_idx2] / (_pred[4 * cell_index + cluster_idx2] + _pred[4 * cell_index + cluster_idx1]))
    auc = roc_auc_score(y_true, y_scores)
    fpr, tpr, thresholds = roc_curve(y_true, y_scores, pos_label=1)
    return auc, fpr, tpr

def evaluate_performance():
    cluster_name_list = ['A', 'B', 'C', 'D']
    makedirs(CURRENT_DIR + 'results/')
    fw = open(CURRENT_DIR + 'results/' + 'evaluate_roc.performance.txt', 'w')
    fw.write('Category\tFPR\tTPR\n')
    for cluster_idx1 in range(0, 4):
        for cluster_idx2 in range(0, 4):
            if cluster_idx1 >= cluster_idx2:
                continue
            auc, fpr, tpr = get_auc_calculation(cluster_idx1, cluster_idx2)
            auc = round(auc, 2)
            cluster_name1 = cluster_name_list[cluster_idx1]
            cluster_name2 = cluster_name_list[cluster_idx2]
            print 'cluster{} vs. cluster{}, AUC = {}'.format(cluster_name1, cluster_name2, auc)
            comparison_info = '{} vs. {} ({})'.format(cluster_name1, cluster_name2, auc)
            for _fpr, _tpr in zip(fpr, tpr):
                fw.write('{}\t{}\t{}\n'.format(comparison_info, _fpr, _tpr))
    fw.close()

def main():
    print 'making training data matrix'
    make_train_data()
    random_seed = 1037
    # split the train data into five folds. Iteratively, four as training, the rest as testing.
    for fold_id in [1, 2, 3, 4, 5]:
        print '...getting data for fold{}'.format(fold_id)
        split_train_data(random_seed, fold_id)
        # train data
        print '...train fold{}'.format(fold_id)
        deep_learning_train(random_seed, fold_id)
        # test data
        print '...predict fold{}'.format(fold_id)
        deep_learning_pred(fold_id)
    print 'evaluate the prediction from deep learning'
    evaluate_performance()

if __name__=="__main__":
    main()
