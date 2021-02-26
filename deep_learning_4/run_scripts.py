## Deep learning of intron clusters
import os
import deep_learning_train

# make training data matrix
deep_learning_train.make_train_data()

random_seed = 1037
# split the train data into five folds. Iteratively, four as training, the rest as testing.
for fold_id in [1, 2, 3, 4, 5]:
    deep_learning_train.split_train_data(random_seed, fold_id)

# train data fold 1
deep_learning_train.deep_learning_train(random_seed, 1)

# train data fold 2
deep_learning_train.deep_learning_train(random_seed, 2)

# train data fold 3
deep_learning_train.deep_learning_train(random_seed, 3)

# train data fold 4
deep_learning_train.deep_learning_train(random_seed, 4)

# train data fold 5
deep_learning_train.deep_learning_train(random_seed, 5)

# predict testing data
for fold_id in [1, 2, 3, 4, 5]:
    deep_learning_train.deep_learning_pred(fold_id)

# evaluate the performance of deep learning
deep_learning_train.evaluate_performance()

os.system('Rscript draw.roc.R')
