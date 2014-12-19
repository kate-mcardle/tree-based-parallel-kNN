/********************************************************************
 * File: Galois_KD_tree.cpp
 * Author: Kate McArdle
 *
 * A file to find the k nearest neighbors of an input graph, using Galois.
 *
 * To use:
 * 	- Initial setup:
 * 		- Note: Ensure that your system has: GCC 4.7 or later; cmake; boost
 * 		- Create a directory (eg 'kdtree') in the 'apps' directory found in your Galois directory (eg Galois-2.2.1/apps/kdtree)
 * 		- Put this file in that new directory
 * 		- Create a file CMakeLists.txt in that new directory. If your app directory is 'kdtree', the CMakeLists.txt file would look like:
 * 				if(USE_EXP)
  				include_directories(../../exp/apps/kdtree .)
					endif()
					app(kdtree Galois_KD_tree.cpp)
 * 		- Add your app directory to the CMakeLists.txt file in the 'apps' directory. If your app directory is 'kdtree', add this line:
 * 				add_subdirectory(kdtree)
 * 		- In Galois/build/default, issue a make command. (In the future, if you modify this file and wish to rebuild it, issue: make -C apps/kdtree)
 * 	- Run on data:
 * 		Galois/build/default/apps/kdtree/kdtree <path to data> n_dimensions n_data_points k n_threads
 * 			- Note: this assumes the data file is in a binary double precision format.
 * 				If it's in another format, see the enum "FileFormat" and make the appropriate change to the read_data call in main.
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <algorithm>
#include <queue>
#include <map>
#include <cmath>
// for file reading:
#include <cstring>
#include <fstream>
#include <sstream>

#include "Galois/Galois.h"
#include "Galois/Graph/FirstGraph.h"
#include "Galois/Statistic.h"

using namespace std;

struct KDNode {
	double* pt;
	int idx;
	int dim;
	bool isLeftChild;

	KDNode(double* pt, int idx) {
		this->pt = pt;
		this->idx = idx;
		dim = -1;
		isLeftChild = true;
	}
};

struct CompareNodes {
	int dim;
	CompareNodes(int dim) : dim(dim) {}

	bool operator() (const KDNode& i, const KDNode& j) {
		return (i.pt[dim] < j.pt[dim]);
	}
};

typedef Galois::Graph::FirstGraph<KDNode,void,true> KDTree;
typedef KDTree::GraphNode TreeNode;
typedef pair<TreeNode, TreeNode> TreeEdge;

struct P_buildTree {
	int gnode_parent_idx;
	int kd_lo_idx;
	int kd_hi_idx;
	int dim;
	bool isLeftChild;

	P_buildTree(int gnode_parent_idx, int kd_lo_idx, int kd_hi_idx, int dim, bool isLeftChild) {
		this->gnode_parent_idx = gnode_parent_idx;
		this->kd_lo_idx = kd_lo_idx;
		this->kd_hi_idx = kd_hi_idx;
		this->dim = dim;
		this->isLeftChild = isLeftChild;
	}
};

int find_split_dim(vector<KDNode>& kdnodes, int lo_idx, int hi_idx, int n_dims) {
	int best_dim = 0;
	double max_var = 0;
	int n = hi_idx - lo_idx;

	for (int d = 0; d < n_dims; ++d) {
		double mean = 0;
		double sum_squares = 0;
		for (int i = lo_idx; i < hi_idx; ++i) {
			mean += kdnodes[i].pt[d];
			sum_squares += (kdnodes[i].pt[d])*(kdnodes[i].pt[d]);
		}
		mean /= n;
		double var = (sum_squares/n) - mean*mean;
		if (var > max_var) {
			max_var = var;
			best_dim = d;
		}
	}
	return best_dim;
}

struct P_addTreeEdge {
	vector<KDNode>& kdnodes;
	KDTree& tree;
	TreeNode* gnodes;
	int n_dims;

	P_addTreeEdge(vector<KDNode>& kdnodes, KDTree& tree, TreeNode* gnodes, int n_dims) : kdnodes(kdnodes), tree(tree) {
		this->gnodes = gnodes;
		this->n_dims = n_dims;
	}

	void operator() (P_buildTree& curr, Galois::UserContext<P_buildTree>& wl) {
		// base cases:
		if (curr.kd_hi_idx == curr.kd_lo_idx) { return; }
		if ((curr.kd_hi_idx - curr.kd_lo_idx) == 1) {
			int gnode_idx = kdnodes[curr.kd_lo_idx].idx;
			tree.addEdge(gnodes[curr.gnode_parent_idx], gnodes[gnode_idx]);
			KDNode& n = tree.getData(gnodes[gnode_idx]);
			n.isLeftChild = curr.isLeftChild;
			n.dim = curr.dim;
			return;
		}

		int split_dim = curr.dim;
		/* if splitting on max discrimination */
		split_dim = find_split_dim(kdnodes, curr.kd_lo_idx, curr.kd_hi_idx, n_dims);
		/* end if splitting on max discrimination */

		sort(kdnodes.begin()+curr.kd_lo_idx, kdnodes.begin()+curr.kd_hi_idx, CompareNodes(split_dim));
		int med_idx = curr.kd_lo_idx + (curr.kd_hi_idx-curr.kd_lo_idx)/2;
		int gnode_idx = kdnodes[med_idx].idx;
		wl.push(P_buildTree(gnode_idx, curr.kd_lo_idx, med_idx, (curr.dim+1)%n_dims, true));
		wl.push(P_buildTree(gnode_idx, med_idx+1, curr.kd_hi_idx, (curr.dim+1)%n_dims, false));
		tree.addEdge(gnodes[curr.gnode_parent_idx], gnodes[gnode_idx]);
		KDNode& n = tree.getData(gnodes[gnode_idx]);
		n.isLeftChild = curr.isLeftChild;
		n.dim = split_dim;
	}
};

class CompareDist {
public:
	bool operator() (const pair<int, double>& lhs, const pair<int, double>& rhs) const {
		if (lhs.second < rhs.second) { return true; }
		return false;
	}
};

typedef Galois::Graph::FirstGraph<KDNode,double,true> KNNGraph;
typedef KNNGraph::GraphNode KNNNode;

struct P_knn {
	int k;
	int n_dims;
	KNNGraph& knnGraph;
	KNNNode* knnNodes;
	KDTree& kdtree;
	TreeNode& kdroot;

	P_knn(int k, int D, KNNGraph& knnGraph, KNNNode* knnNodes, KDTree& kdtree, TreeNode& root) : knnGraph(knnGraph), knnNodes(knnNodes), kdtree(kdtree), kdroot(root) {
		this->k = k;
		this->n_dims = D;
	}

	double getDistance(const KDNode& key, KDNode& curr) {
		double dist = 0.0;
		for (int i = 0; i < n_dims; ++i) {
			dist += (key.pt[i] - curr.pt[i]) * (key.pt[i] - curr.pt[i]);
		}
		return sqrt(dist);
	}

	void search_subtree(priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist>& pq, TreeNode& curr, const KDNode& key_kd, int level) {
		KDNode& curr_kd = kdtree.getData(curr, Galois::MethodFlag::NONE);
		double dist = getDistance(key_kd, curr_kd);
		// if pq has less than k elems in it, push curr on:
		if (pq.size() < k) {
			pq.emplace((make_pair(curr_kd.idx, dist)));
		}
		// otherwise, only push on if distance of curr is less than distance of max elem, and pop max elem:
		else if (dist < pq.top().second) {
			pq.pop();
			pq.emplace((make_pair(curr_kd.idx, dist)));
		}
		bool left_child_searched = true;
		if (key_kd.pt[curr_kd.dim] < curr_kd.pt[curr_kd.dim]) {
			// check to see if curr has a left child. if it does, search_subtree on that child
			for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
				TreeNode dest = kdtree.getEdgeDst(edge);
				if (kdtree.getData(dest, Galois::MethodFlag::NONE).isLeftChild) {
					search_subtree(pq, dest, key_kd, (level+1)%n_dims);
					break;
				}
			}
		}
		else {
			left_child_searched = false;
			// check to see if curr has a right child. if it does, search_subtree on that child
			for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
				TreeNode dest = kdtree.getEdgeDst(edge);
				if (!kdtree.getData(dest, Galois::MethodFlag::NONE).isLeftChild) {
					search_subtree(pq, dest, key_kd, (level+1)%n_dims);
					break;
				}
			}
		}

		// as we walk back up the tree:
		if ( (pq.size() < k) || (fabs(key_kd.pt[curr_kd.dim] - curr_kd.pt[curr_kd.dim]) < pq.top().second) ) {
			if (left_child_searched) { // search right subtree
				for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
					TreeNode dest = kdtree.getEdgeDst(edge);
					if (!kdtree.getData(dest, Galois::MethodFlag::NONE).isLeftChild) {
						search_subtree(pq, dest, key_kd, (level+1)%n_dims);
						break;
					}
				}
			}
			else { // search left subtree
				for (KDTree::edge_iterator edge : kdtree.out_edges(curr)) {
					TreeNode dest = kdtree.getEdgeDst(edge);
					if (kdtree.getData(dest, Galois::MethodFlag::NONE).isLeftChild) {
						search_subtree(pq, dest, key_kd, (level+1)%n_dims);
						break;
					}
				}
			}
		}
	}

	void operator() (const KNNNode& node) {
		priority_queue<pair<int, double>, vector<pair<int, double> >, CompareDist> pq;
		const KDNode& key_kd = knnGraph.getData(node, Galois::MethodFlag::NONE);
		search_subtree(pq, kdroot, key_kd, 0);
		for (int i = 0; i < k; ++i) {
			int dest_idx = pq.top().first;
			KNNNode& knn_dest = knnNodes[dest_idx];
			double dist = pq.top().second;
			pq.pop();
			knnGraph.getEdgeData(knnGraph.addEdge(node, knn_dest, Galois::MethodFlag::NONE), Galois::MethodFlag::NONE) = dist;
//			knnGraph.getEdgeData(knnGraph.addEdge(node, knn_dest)) = dist;
		}
	}
};

/* data in file stored in binary double precision format. Source: Chen */
bool read_data_from_binary_double(double **data, char *filename, int N, int D) {
  FILE *fp = NULL;
  if (!(fp = fopen(filename, "rb"))) {
  	cout << filename << " didn't open" << endl;
  	return false;
  }

  int i;
  int num_in, num_total;
  for (i = 0; i < N; ++i)
    {
      num_total = 0;
      while (num_total < D)
        {
          num_in = fread(data[i]+num_total, sizeof(double), D, fp);
          num_total += num_in;
        }
    }

  fclose(fp);
  return true;
}

/* data in file stored in binary 1-byte uint8 precision format, eg for Tiny Images dataset. */
bool read_data_from_binary_uint(double **data, char *filename, int N, int D) {
  FILE *fp = NULL;
  if (!(fp = fopen(filename, "rb"))) {
  	cout << filename << " didn't open" << endl;
  	return false;
  }

	uint8_t** data_uint;
	data_uint = (uint8_t**) malloc(N*sizeof(uint8_t *));
	for (int i = 0; i < N; ++i) {
		data_uint[i] = (uint8_t*) malloc(D*sizeof(uint8_t));
	}

  int i;
  int num_in, num_total;
  for (i = 0; i < N; ++i) {
  	num_total = 0;
  	while (num_total < D) {
  		num_in = fread(data_uint[i]+num_total, sizeof(uint8_t), D, fp);
  		num_total += num_in;
  	}
  }
  fclose(fp);

  for (i = 0; i < N; ++i) {
  	for (int j = 0; j < D; ++j) {
  		data[i][j] = (double) data[i][j];
  	}
  }
  free(data_uint);
  return true;
}

/* data in file stored in libsvm format */
bool read_data_from_libsvm(double **data, char *filename, int N, int D) {
	for (int idx = 0; idx < N; ++idx) {
		for (int dim = 0; dim < D; ++dim) {
			data[idx][dim] = 0;
		}
	}

	ifstream f(filename);
	string line;
	int data_idx = 0;
	while(getline(f, line)) {
		istringstream iss(line);
		string label = "";
		iss >> label;
		while(iss) {
			string sub = "";
			iss >> sub;
			if(sub.length() == 0) { break; }
			string feature = sub.substr(0, sub.find(':'));
			int dim = atoi(feature.c_str()) - 1; // have to decrement to be 0-based
			string val = sub.substr(sub.find(':')+1, sub.length());
			data[data_idx][dim] = atof(val.c_str());
		}
		++data_idx;
	}
	f.close();
  return true;
}

struct Times {
	double build_tree;
	double get_knn;
	double total;

	Times() {
		build_tree = 0.0;
		get_knn = 0.0;
		total = 0.0;
	}
};

enum FileFormat { BINARY_DOUBLE, BINARY_UINT8, LIBSVM };

double** read_data(char* datafile, FileFormat fileformat, int D, int N) {
	double** data;
	data = (double**) malloc(N*sizeof(double *));
	for (int i = 0; i < N; ++i) {
		data[i] = (double*) malloc(D*sizeof(double));
	}

	bool read_result = false;
	switch(fileformat) {
	case BINARY_DOUBLE :
		read_result = read_data_from_binary_double(data, datafile, N, D);
		break;
	case BINARY_UINT8 :
		read_result = read_data_from_binary_uint(data, datafile, N, D);
		break;
	case LIBSVM :
		read_result = read_data_from_libsvm(data, datafile, N, D);
		break;
	default :
		read_result = false;
	}

	if (!read_result) {
		printf("error reading data!\n");
		exit(1);
	}
	return data;
}

void run_knn(double** data, int D, int N, int k, int n_threads, Times& times) {
	/* Read in data: */
	Galois::StatManager sm;

	Galois::setActiveThreads(1);
	Galois::StatTimer totalTime("totalTime");
	totalTime.start();

	/* Convert to vector of KDNodes and insert into graph: */
	Galois::StatTimer buildTree("buildTree");
	buildTree.start();
	KDTree kdtree;
	vector<KDNode> kdnodes; // vector of kdtree nodes, NOT in indexed order (starts in order, but will be resorted)
	kdnodes.reserve(N);
	TreeNode* gnodes; // array of graph nodes, IN indexed order: gnodes[0] corresponds to the first data point in the file
	gnodes = (TreeNode*) malloc(N*sizeof(TreeNode));
	for (int i = 0; i < N; ++i) {
		kdnodes.emplace_back(data[i], i);
		gnodes[i] = kdtree.createNode(kdnodes[i]);
		kdtree.addNode(gnodes[i], Galois::MethodFlag::NONE);
	}

	/* Build KDTree */
	int split_dim = 0;

	/* if splitting on max discrimination */
	split_dim = find_split_dim(kdnodes, 0, N, D);
	/* end if splitting on max discrimination */

	sort(kdnodes.begin(), kdnodes.end(), CompareNodes(split_dim));
	int median = N/2;
	int root_node_idx = kdnodes[median].idx; // corresponds to the root's index in gnodes: gnodes[root_node_idx]
	KDNode& root = kdtree.getData(gnodes[root_node_idx], Galois::MethodFlag::NONE);
	root.dim = split_dim;
	vector<P_buildTree> worklist;
	worklist.emplace_back(root_node_idx, 0, median, 1, true);
	worklist.emplace_back(root_node_idx, median+1, N, 1, false);
	Galois::for_each(worklist.begin(), worklist.end(), P_addTreeEdge(kdnodes, kdtree, gnodes, D));
	buildTree.stop();
	times.build_tree += buildTree.get();

	/* Build kNN Graph */
	Galois::setActiveThreads(n_threads);
	Galois::StatTimer knnTime("knnTime");
	knnTime.start();
	KNNGraph knnGraph;
	KNNNode* gnodes_knn;
	gnodes_knn = (KNNNode*) malloc(N*sizeof(KNNNode));
	for (int i = 0; i < N; ++i) {
		int idx = kdnodes[i].idx;
		gnodes_knn[idx] = knnGraph.createNode(kdnodes[i]);
		knnGraph.addNode(gnodes_knn[idx], Galois::MethodFlag::NONE);
	}
	Galois::do_all(knnGraph.begin(), knnGraph.end(), P_knn(k, D, knnGraph, gnodes_knn, kdtree, gnodes[root_node_idx]));
//	Galois::do_all_local(knnGraph, P_knn(k, D, knnGraph, gnodes_knn, kdtree, gnodes[root_node_idx]));
	knnTime.stop();
	times.get_knn += knnTime.get();
	totalTime.stop();
	times.total += totalTime.get();

	free(gnodes);
	free(gnodes_knn);
}

struct Data {
	char* datafile;
	FileFormat fileformat;
	int D;
	int N;
	int k;

	Data(char* f, FileFormat fileformat, int D, int N, int k) : datafile(f), fileformat(fileformat), D(D), N(N), k(k) { }
};

void measure_performance() {
	int n_reps = 1;
	vector<Data> datasets;
	int n_threads[] = { 16, 8, 4, 1 };

	Data mnist_test("../../../ScalableML/data/mnist-test.dat", BINARY_DOUBLE, 784, 10000, 8);
	datasets.push_back(mnist_test);
	Data mnist_train("../../../ScalableML/data/mnist-train.dat", BINARY_DOUBLE, 784, 60000, 8);
	datasets.push_back(mnist_train);
	Data mnist8m("/scratch/02234/kmcardle/data/mnist8m", LIBSVM, 784, 8100000, 8);
	datasets.push_back(mnist8m);
	Data tinyimgs("/scratch/02234/kmcardle/data/tiny_images.bin", BINARY_UINT8, 3072, 79302017, 8);
	datasets.push_back(tinyimgs);

	Data poker("/scratch/02234/kmcardle/data/poker.t", LIBSVM, 10, 1000000, 8);
	datasets.push_back(poker);
	Data rna("/scratch/02234/kmcardle/data/cod-rna.t", LIBSVM, 8, 271617, 8);
	datasets.push_back(rna);
	Data cadata("/scratch/02234/kmcardle/data/cadata", LIBSVM, 8, 20640, 8);
	datasets.push_back(cadata);
	Data covtype("/scratch/02234/kmcardle/data/covtype.libsvm.binary", LIBSVM, 54, 581012, 8);
	datasets.push_back(covtype);
	Data year("/scratch/02234/kmcardle/data/YearPredictionMSD", LIBSVM, 90, 463715, 8);
	datasets.push_back(year);
	Data aloi("/scratch/02234/kmcardle/data/aloi", LIBSVM, 128, 108000, 8);

	for (Data dataset : datasets) {
		cout << "---------- Results for " << dataset.datafile << " data: ----------" << endl;
		double** data = read_data(dataset.datafile, dataset.fileformat, dataset.D, dataset.N);
		for (int t : n_threads) {
			printf(" +++++ %d Threads +++++\n", t);
			Times times;
			for (int rep = 0; rep < n_reps; ++rep) {
				run_knn(data, dataset.D, dataset.N, dataset.k, t, times);
			}
			printf("\n\n +++++ %d Threads +++++\n", t);
			printf("average time to build kd tree = %f seconds\n", times.build_tree/(n_reps*1000));
			printf("average time to get knn graph = %f seconds\n", times.get_knn/(n_reps*1000));
			printf("average total time            = %f seconds\n\n\n", times.total/(n_reps*1000));
		}
		free(data);
	}

}

int main(int argc, char **argv) {
	if (argc == 6) {
		char* datafile = argv[1];
		cout << "Results for " << datafile << "data:" << endl;
		int D = atoi(argv[2]);
		int N = atoi(argv[3]);
		int k = atoi(argv[4]);
		int n_threads = atoi(argv[5]);
		double** data = read_data(datafile, BINARY_DOUBLE, D, N);
		Times times;
		run_knn(data, D, N, k, n_threads, times);
		free(data);
	}
	else {
		measure_performance();
	}
}