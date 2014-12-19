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
#include "kdtree.h"

struct kdnode {
	double *pos;
	int dir;
	int *data;
	struct kdnode *left, *right;	/* negative/positive side */
};

struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};

#define SIZE 100

struct rheap {
    struct res_node* node;
    //double* heaparray;
    int capacity;
    int size;
};

#define SQ(x)			((x) * (x))



struct kdtree *kd_create(int k)
{
	struct kdtree *tree;

	if(!(tree = malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = k;
	tree->root = 0;
	tree->destr = 0;
	tree->rect = 0;

	return tree;
}

static int insert_rec(struct kdnode **nptr, const double *pos, int data, int dir, int dim)
{
	int new_dir;
	struct kdnode *node;

	if(!*nptr) {
		if(!(node = malloc(sizeof *node))) {
			return -1;
		}
		if(!(node->pos = malloc(dim * sizeof *node->pos))) {
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->data = data;
		node->dir = dir;
		node->left = node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

int kd_insert(struct kdtree *tree, const double *pos, int data)
{
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) {
		return -1;
	}
	return 0;
}
bool sortNode(struct kdnode* i,struct kdnode* j, int dim){
    return i->pos[dim] < j->pos[dim];
}
static int insert_rec_all(struct kdnode *nptr, char* side, struct kdnode* nodes , int start, int end, int dim)
{
    if(!start == end){
        if ((end - start) == 1) {
            if (strcmp(side, "left")==0) {
                nptr->left=&nodes[start];
            }
            else if (strcmp(side, "right")==0) {
                nptr->right=&nodes[start];
            }
            return 0;
        }
        struct kdnode* new_node = &nodes[start];
        qsort_r(new_node, end-start,sizeof(nodes[0]),sortNode,dim);
        int median = start + (end-start)/2;
        if (strcmp(side, "left")==0) {
            nptr->left=&nodes[median];
            insert_rec_all(nptr->left, "left", nodes, start, median, dim+1);
            insert_rec_all(nptr->right, "right", nodes, median+1, end, dim+1);
        }
        else if (strcmp(side, "right")==0) {
           nptr->right=&nodes[median];
            insert_rec_all(nptr->left, "left", nodes, start, median, dim+1);
            insert_rec_all(nptr->right, "right", nodes, median+1, end, dim+1);
        }
    }
    return 0;
}



int kd_insert_all(struct kdtree *tree, const double **pos_all, int *data, int n, int D)
{   struct kdnode *nodes[n];
    int j;
    
    for (j=0; j<n; j++) {
        struct kdnode *node;
        
        if(!(node = malloc(sizeof *node))) {
            return -1;
        }
        if(!(node->pos = malloc(D * sizeof *node->pos))) {
            free(node);
            return -1;
        }
        
        memcpy(node->pos, &pos_all[j], D * sizeof *node->pos);
        node->data = data[j];
        
        node->left = node->right = 0;
        nodes[j] = node;
        free(node);
    }
       qsort_r(nodes, n, sizeof(nodes[0]),sortNode,0);
    tree->root = nodes[n/2];
    if (insert_rec_all(&tree->root,"left",&nodes, 0, n/2,1) && insert_rec_all(&tree->root, "right",&nodes,n/2+1 ,n, 1)) {
        return -1;
    }
    
    return 0;
}


struct rheap* initHeap() {
    struct rheap* h;
    h = (struct rheap*)(malloc(sizeof(struct rheap)));
    h->capacity = SIZE;
    //h->heaparray = (double*)malloc(sizeof(double)*(SIZE+1));
    h->size = 0;
    h->node=(struct res_node*)malloc(sizeof(struct res_node)*(SIZE+1));
    return h;
}
void swap(struct rheap *h, int index1, int index2) {
    struct res_node temp = h->node[index1];
    h->node[index1] = h->node[index2];
    h->node[index2] = temp;
}
void percolateUp(struct rheap *h, int index) {
    if (index > 1) {
        if (h->node[index/2].dist_sq < h->node[index].dist_sq) {
            swap(h, index, index/2);
            percolateUp(h, index/2);
        }
    }
}
int maximum(double a, int indexa, double b, int indexb) {
    if (a > b)
    return indexa;
    else
    return indexb;
}
void percolateDown(struct rheap *h, int index) {
    int max;
    if ((2*index+1) <= h->size) {
        max = maximum(h->node[2*index].dist_sq, 2*index, h->node[2*index+1].dist_sq, 2*index+1);
        if (h->node[index].dist_sq > h->node[max].dist_sq) {
            swap(h, index, max);
            percolateDown(h, max);
        }
    }
    else if (h->size == 2*index) {
        if (h->node[index].dist_sq > h->node[2*index].dist_sq)
        swap(h, index, 2*index);
    }
}

int rheap_insert(struct rheap *h, struct knode *node,double value) {
    struct res_node* tempnode;
    struct res_node* thrownode;
    int i;
    if (h->size == h->capacity) {
        h->capacity *= 2;
        tempnode=(struct res_node*)malloc(sizeof(struct res_node)*h->capacity+1);
        for (i=0; i<h->capacity; i++){
            tempnode[i] = h->node[i];
        }
        thrownode= h->node;
        h->node = thrownode;
        free(thrownode);
    }
    h->size++;
    h->node[h->size].item=node;
    h->node[h->size].dist_sq=value;
    percolateUp(h, h->size);
    return 1;
}

struct res_node* rheap_remove_max(struct rheap *h) {
    struct res_node* retval;
    if (h->size > 0) {
        retval = &h->node[1];
        h->node[1] = h->node[h->size];
        h->size--;
        percolateDown(h, 1);
        return &retval;
    }
}
void heapify(struct rheap *h) {
    int i;
    for (i=h->size/2; i>0; i--)
    percolateDown(h, i);
    
}

struct rheap * initHeapfromArray(struct res_node* values, int length) {
    int i;
    struct rheap* h;
    h = (struct rheap*)(malloc(sizeof(struct rheap)));
    h->node = (struct res_node*)malloc(sizeof(struct res_node)*(length+1));
    for (i=0; i<length; i++)
    h->node[i] = values[i];
    h->size = length;
    heapify(h);
    return h;
}

void rheap_heapsort(struct res_node* values[], int length) {
    struct rheap *h;
    int i;
    h =  initHeapfromArray(values, length);
    length = h->size;
    for (i=0; i<length; i++) {
        values[i] = rheap_remove_max(h);
    }
}

struct res_node* rheap_get_max(struct rheap *h){
    return &h->node[1];
}

static int find_nearest_n(struct kdnode *node, const double *pos, double range, int num, struct rheap *heap, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;
    double range_sq = SQ(range);
	if(!node) return 0;
    
	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= range_sq) {
		if(heap->size >= num) {
			/* get furthest element */
			struct res_node *maxelem = rheap_get_max(heap);

			/* and check if the new one is closer than that */
			if(maxelem->dist_sq > dist_sq) {
				rheap_remove_max(heap);

				if(rheap_insert(heap, node, dist_sq) == -1) {
					return -1;
				}
				added_res = 1;

				range_sq = dist_sq;
			}
		} else {
			if(rheap_insert(heap, node, dist_sq) == -1) {
				return -1;
			}
			added_res = 1;
		}
	}

	/* find signed distance from the splitting plane */
	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest_n(dx <= 0.0 ? node->left : node->right, pos, range, num, heap, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest_n(dx <= 0.0 ? node->right : node->left, pos, range, num, heap, dim);
	}

}

void *kd_nearest_n(struct kdtree *kd, const double *pos, int num)
{
	int ret;
    struct rheap *h;
    h = initHeap();
	if((ret = find_nearest_n(kd->root, pos, range, num, h, kd->dim)) == -1) {
		return 0;
	}
}


