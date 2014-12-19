/*************************************************
 *
 * Author: Nazneen Rajani
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <vector>
#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <string>
#include <iostream>
#include <time.h>
#include <omp.h>
#include<queue>
#include<cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include<limits>

using namespace std;


int read_X_from_file(double **X, int n, int D, char *filename);

struct datapoint {
    double *dim; //list of dimensions
    int idx;  //its own index
    datapoint(double *dim, int idx){
        this -> dim = dim;
        this -> idx = idx;
    }
};

vector<datapoint> pts;

int D;

struct ballnode {
    vector<struct datapoint*> data; //list of nodes it owns
    double radius;
    double *pivot;
    struct ballnode* child1;
    struct ballnode* child2;
    ballnode(vector<struct datapoint*> data){
        this -> data = data;
        this->radius=0;
        this->pivot= (double*) calloc(D,sizeof(double));
        this->child1=NULL;
        this->child2=NULL;
    }
};

struct balltree {
    int num_points;
    struct ballnode *root;
    balltree(int n){
        this->num_points=n;
        this->root=NULL;
    }
};
struct balltree *btree;
pair<double,struct datapoint*> getRadius(double *target, vector<datapoint*>& pts, int D) {
    double radius =0.0;
    struct datapoint *child1;
    for (int k=0; k<pts.size(); k++) {
        double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target[i] - pts[k]->dim[i]) * (target[i] - pts[k]->dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child1 = pts[k];
        }
    }
    return make_pair(radius,child1);
}

struct datapoint* getMaxDist(const struct datapoint *target, vector<datapoint*>& pts, int D) {
    double radius =0.0;
    struct datapoint *child2;
    for (int k=0; k<pts.size(); k++) {
        double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target->dim[i] - pts[k]->dim[i]) * (target->dim[i] - pts[k]->dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child2 = pts[k];
        }
    }
    return child2;
}

double getDistance(struct datapoint* key, struct datapoint* curr, int D) {
    double dist = 0.0;
    for (int i = 0; i < D; ++i) {
        dist += (key->dim[i] - curr->dim[i]) * (key->dim[i] - curr->dim[i]);
    }
    return sqrt(dist);
}

double getDistancePivot(struct datapoint* key, double* curr, int D) {
    double dist = 0.0;
    for (int i = 0; i < D; ++i) {
        dist += (key->dim[i] - curr[i]) * (key->dim[i] - curr[i]);
        //cout<<curr[i]<<endl;
    }
    return sqrt(dist);
}

void recursive_insert(struct ballnode *root, int items, int D, int leaf_size){
    double dim[D];
    for (int d=0; d<D; d++) {
        dim[d]=0;
    }
    for (int k=0; k<D; k++) {
        for (int i=0; i<items; i++) {
            dim[k] += root->data.at(i)->dim[k];
        }
        dim[k] = dim[k]/items;
        root->pivot[k] = dim[k];
    }
    pair<double, struct datapoint*> answer = getRadius(root->pivot,root->data,D);
    root->radius = answer.first;
    
    struct datapoint *child1_point = answer.second;
    struct datapoint *child2_point = getMaxDist(child1_point,root->data,D);
    
    vector<datapoint*> child1_points;
    vector<datapoint*> child2_points;
    
    for (int k=0; k<items; k++) {
        datapoint *dt = root->data.at(k);
        if(getDistance(child1_point,dt,D)<=getDistance(child2_point,dt,D))
            child1_points.push_back(dt);
        else
            child2_points.push_back(dt);
    }
    
    struct ballnode *child1 = new ballnode(child1_points);
    struct ballnode *child2 = new ballnode(child2_points);
    
    root->child1 = child1;
    root->child2 = child2;
    
    if(child1_points.size()>leaf_size) //I might change
        recursive_insert(root->child1, child1_points.size(),D,leaf_size); //recursive insert on both childs
    else{
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<root->child1->data.size(); i++) {
                dim[k] += root->child1->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/root->child1->data.size();
            root->child1->pivot[k] = dim[k];
        }
        answer = getRadius(root->child1->pivot,root->child1->data,D);
        root->child1->radius = answer.first;
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<root->child2->data.size(); i++) {
                dim[k] += root->child2->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/root->child2->data.size();
            root->child2->pivot[k] = dim[k];
        }
        answer = getRadius(root->child2->pivot,root->child2->data,D);
        root->child2->radius = answer.first;
    }
    if(child2_points.size()>leaf_size) //I might change
    {
        recursive_insert(root->child2, child2_points.size(),D,leaf_size);
        
    }
    else{
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<root->child1->data.size(); i++) {
                dim[k] += root->child1->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/root->child1->data.size();
            root->child1->pivot[k] = dim[k];
        }
        answer = getRadius(root->child1->pivot,root->child1->data,D);
        root->child1->radius = answer.first;
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<root->child2->data.size(); i++) {
                dim[k] += root->child2->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/root->child2->data.size();
            root->child2->pivot[k] = dim[k];
        }
        answer = getRadius(root->child2->pivot,root->child2->data,D);
        root->child2->radius = answer.first;
    }
}

struct ballnode* balltree_insert_all(struct balltree* tree, double **dim_all, int n, int D, int leaf_size)
{
    int j;
    pts.reserve(n);
    for (j=0; j<n; j++) {
        pts.emplace_back(dim_all[j],j);
    }
    vector<datapoint*> points;
    points.reserve(n);
    for (j=0; j<n; j++) {
        points.push_back(&pts.at(j));
    }
    struct ballnode *root = new ballnode(points);
    tree->root = root;
    double dim[D];
    
    for (int d=0; d<D; d++) {
        dim[d]=0;
    }
    for (int k=0; k<D; k++) {
        for (int i=0; i<n; i++) {
            dim[k] += points.at(i)->dim[k];
        }
        dim[k] = dim[k]/n;
        tree->root->pivot[k] = dim[k];
    }
    
    pair<double, struct datapoint*> answer = getRadius(tree->root->pivot,tree->root->data,D);
    tree->root->radius = answer.first;
    
    struct datapoint* child1_point = answer.second;
    struct datapoint* child2_point = getMaxDist(child1_point,tree->root->data,D);
    
    vector<datapoint*> child1_points;
    vector<datapoint*> child2_points;
    
    for (int k=0; k<points.size(); k++) {
        datapoint *dt = points.at(k);
        if(getDistance(child1_point,dt,D)<=getDistance(child2_point,dt,D))
            child1_points.push_back(dt);
        else
            child2_points.push_back(dt);
    }
    
    struct ballnode *child1 = new ballnode(child1_points);
    struct ballnode *child2 = new ballnode(child2_points);
    
    tree->root->child1 = child1;
    tree->root->child2 = child2;
    
    if(child1_points.size()>leaf_size) // might change
        recursive_insert(tree->root->child1, child1_points.size(),D,leaf_size); //recursive insert on both childs

    else{
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<tree->root->child1->data.size(); i++) {
                dim[k] += tree->root->child1->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/tree->root->child1->data.size();
            tree->root->child1->pivot[k] = dim[k];
        }
        answer = getRadius(tree->root->child1->pivot,tree->root->child1->data,D);
        tree->root->child1->radius = answer.first;
        
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<tree->root->child2->data.size(); i++) {
                dim[k] += tree->root->child2->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/tree->root->child2->data.size();
            tree->root->child2->pivot[k] = dim[k];
        }
        answer = getRadius(tree->root->child2->pivot,tree->root->child2->data,D);
        tree->root->child2->radius = answer.first;
        
    }
    if(child2_points.size()>leaf_size) //might change
        recursive_insert(tree->root->child2, child2_points.size(),D,leaf_size);
    else{
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<tree->root->child1->data.size(); i++) {
                dim[k] += tree->root->child1->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/tree->root->child1->data.size();
            tree->root->child1->pivot[k] = dim[k];
        }
         answer = getRadius(tree->root->child1->pivot,tree->root->child1->data,D);
        tree->root->child1->radius = answer.first;
        
        dim[D];
        for (int d=0; d<D; d++) {
            dim[d]=0;
        }
        for (int k=0; k<D; k++) {
            for (int i=0; i<tree->root->child2->data.size(); i++) {
                dim[k] += tree->root->child2->data.at(i)->dim[k];
            }
            dim[k] = dim[k]/tree->root->child2->data.size();
            tree->root->child2->pivot[k] = dim[k];
        }
        answer = getRadius(tree->root->child2->pivot,tree->root->child2->data,D);
        tree->root->child2->radius = answer.first;
    }
    return tree->root;
}

class sortNodes {
public:
    bool operator() (const pair<int, double>& target, const pair<int, double>& query) {
        if (target.second < query.second)
            return true;
        else
            return false;
    }
};

void balltree_nearest_n(priority_queue<pair<int, double>, vector<pair<int, double> >, sortNodes>& pq, struct ballnode* node, struct datapoint *t, int k, int D,double d_parent){
    double d_node = max(getDistancePivot(t,node->pivot,D)-node->radius,d_parent);
    if(node->child1==NULL && node->child2 == NULL){
        for (int i=0; i<node->data.size(); i++) {
            double dist =getDistance(t, node->data.at(i),D);
            if (pq.size()==k){
                if(dist< pq.top().second){
                    pq.pop();
                    pq.emplace(make_pair(node->data.at(i)->idx,dist));
                }
            }
            else
                    pq.emplace(make_pair(node->data.at(i)->idx,dist));
        }
        
    }
    else if (pq.size()==k){
        if (d_node>=pq.top().second) {
         return;
        }
     }
    else{
        double dist_1 = getDistancePivot(t, node->child1->pivot,D);
        double dist_2 = getDistancePivot(t, node->child2->pivot,D);
        if(dist_1<=dist_2){
            balltree_nearest_n(pq,node->child1,t,k,D,d_node);
            balltree_nearest_n(pq,node->child2,t,k,D,d_node);
        }
        else{
            balltree_nearest_n(pq,node->child2,t,k,D,d_node);
            balltree_nearest_n(pq,node->child1,t,k,D,d_node);
        }
    }
}

int main(int argc, char **argv)
{   int i, j,leaf_size;
    unsigned int msec, start;
    double end, start1;
    int n = atoi(argv[1]);                 // # data points
    D = atoi(argv[2]);                  // data dimension
    char *fname = argv[3];
    double **X;
    X=(double**)malloc(sizeof(double*)*n);
    for (i=0; i<n; i++) {
        X[i] = (double*)malloc(sizeof(double)*D);
    }
    
    if (read_X_from_file(X, n, D, fname)==0)
    {
        printf("error open file.\n");
        exit(1);
    }
    omp_set_num_threads(atoi(argv[4]));
    leaf_size = 0.2*n;
    //atoi(argv[5]);
    double start_build = omp_get_wtime();
    btree = new balltree(n);
    struct ballnode* bnode = balltree_insert_all(btree,X,n,D,leaf_size);

    double end_build = omp_get_wtime()-start_build;
    cout<<"Time to create ball tree: "<<end_build<<endl;
    double start_query = omp_get_wtime();
    
    #pragma omp parallel for default(shared) private(i)
    for (i=0; i<n; i++) {
        priority_queue<pair<int, double>, vector<pair<int, double> >, sortNodes> pq;
        struct datapoint *target = new datapoint(X[i],i);
        balltree_nearest_n(pq, btree->root, target,9,D,0);
        /*while (!pq.empty()) {
            cout<<pq.top().first+1<<"\t"<<pq.top().second<<endl;
            pq.pop();
        }*/
    }
    double end_query = omp_get_wtime()-start_query;
    cout<<"Time to query ball tree: "<<end_query<<endl;
    return 0;
}

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
            int dim = atoi(feature.c_str()) - 1;
            string val = sub.substr(sub.find(':')+1, sub.length());
            data[data_idx][dim] = atof(val.c_str());
        }
        ++data_idx;
    }
    f.close();
    return true;
}

int read_X_from_file(double **X, int n, int D, char *filename)
{
    FILE *fp = NULL;
    if (!(fp = fopen(filename, "rb")))
        return 0;
    int i;
    int num_in, num_total;
    for (i = 0; i < n; i++)
    {
        num_total = 0;
        while (num_total < D)
        {
            num_in = fread(X[i]+num_total, sizeof(double), D, fp);
            num_total += num_in;
        }
    }
    
    fclose(fp);
    return 1;
}
