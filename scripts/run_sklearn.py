import numpy as np
from sklearn.datasets import load_svmlight_file
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import time
import sys

# Usage:
# - From command line: 
#     python run_sklearn.py path_to_input_data_file data_format k d N algorithm
#     example: python run_sklearn.py ../data/frey.dat binary 12 560 1965 kd_tree
# - In interpreter: 
#     knn_graph = run_sklearn.get_knn_graph(path_to_input, data_format, k, d, N, algorithm)
#     example: knn_graph = run_sklearn.get_knn_graph('../data/frey.dat', 'binary', 12, 560, 1965, 'kd_tree')

def get_knn_graph(data_file, data_format, k, d, N, alg):
  if data_format == "binary":
    a = np.fromfile(data_file, dtype=float).reshape((N,d))
  elif data_format == "libsvm":
    x, labels = load_svmlight_file(data_file)
    del labels
    a = x.todense()
    del x
  else:
    print "wrong data format!"
    return 0
  k_plus_1 = k+1
  t_start = time.time()
  nbrs = NearestNeighbors(n_neighbors=(k_plus_1), algorithm=alg, leaf_size=1).fit(a)
  t_tree = time.time()
  knn_graph = nbrs.kneighbors_graph(a)
  t_graph = time.time() - t_tree
  t = time.time() - t_start
  print 'overall time = ' + str(t) + " seconds"
  return knn_graph

if __name__ == '__main__':
  knn_graph = get_knn_graph(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])