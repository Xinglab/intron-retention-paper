from keras.models import model_from_json
from keras.preprocessing import sequence
from keras.optimizers import RMSprop, SGD
from keras.models import Sequential, Model
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers import Input, merge
from keras.layers.normalization import BatchNormalization
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras.utils import np_utils
from keras import backend as K
from keras.callbacks import Callback
from sklearn.decomposition import PCA, FastICA, LatentDirichletAllocation
from subprocess import *
from keras import optimizers
import numpy as np
import random
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error

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

def make_train_data(random_seed, fold_id):
    X = np.load('X.npy')
    y = np.load('y.npy')
    instances = range(0, X.shape[0]) # get sample index
    [random.shuffle(instances) for x in range(7)] # random shuffle instances 7 times
    folds_num = 5 # split the data into 5 folds
    train_list, test_list = split_by_fold(instances, fold_id, folds_num)
    X_train, y_train = [], []
    X_test, y_test = [], []
    for train_index in train_list:
        X_train.append(X[train_index])
        y_train.append(y[train_index])
    for test_index in test_list:
        X_test.append(X[test_index])
        y_test.append(y[test_index])
    return np.array(X_train), np.array(y_train), np.array(X_test), np.array(y_test)

def get_auc(y_pred, y_test, v1, v2):
    acc, count = 0.0, 0.0
    y_true, y_scores = [], []
    for _pred, _test in zip(y_pred, y_test):
        for x in xrange(3):
            if np.isnan(_test[3 * x]):
                continue
            argmax_test = np.argmax(_test[3 * x: 3 * x + 3])
            argmax_pred = np.argmax(_pred[3 * x: 3 * x + 3])
            _pred_list = _pred[3 * x : 3 * x + 3]
            if argmax_test == argmax_pred:
                acc += 1
            count += 1
            if _test[3 * x + v2] == 1:
                y_true.append(1)
                y_scores.append(_pred[3 * x + v2] / (_pred[3 * x + v2] + _pred[3 * x + v1]))
            elif _test[3 * x + v1] == 1:
                y_true.append(0)
                y_scores.append(_pred[3 * x + v2] / (_pred[3 * x + v2] + _pred[3 * x + v1]))
    auc = roc_auc_score(y_true, y_scores)
    return auc

random_seed = 1337
fold_id = 1
epoch = 200
batch_size = 2000

X_train, y_train, X_test, y_test = make_train_data(random_seed, fold_id)
FEATURE_NUM = X_train.shape[1]

sgd = optimizers.SGD(lr=0.01, decay=1e-7, momentum=0.9, nesterov=True)
model = create_model(12, FEATURE_NUM)
model.compile(loss=masked_cross_entropy_loss, optimizer=sgd)
model.fit(X_train, y_train, batch_size=batch_size, epochs=epoch, shuffle=True)
model_json_file = open('model.{}'.format(fold_id) + '.json', 'w')
model_json_file.write(model.to_json())
model_json_file.close()
model.save_weights('model.{}'.format(fold_id) + '.weight.hd5')

y_pred = model.predict(X_test)

for i in range(0, 4):
    for j in range(0, 4):
        if i >= j:
            continue
        auc = round(get_auc(y_pred, y_test, i, j), 2)
        print 'Cluster{} vs. Cluster{} (AUC = {})'.format(i, j, auc)