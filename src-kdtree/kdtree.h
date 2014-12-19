/*************************************************
 *
 * Author: Nazneen Rajani
 *
 */

struct kdtree *kd_create(int k);

int kd_insert(struct kdtree *tree, const double *pos, int data);
struct kdres *kd_nearest_n(struct kdtree *tree, const double *pos, int num);
int *kd_res_item_data(struct kdres *set);

struct res_node* rheap_get_max(struct rheap *h);
void rheap_heapsort(struct res_node* values[], int length); 
struct rheap * initHeapfromArray(struct res_node* values, int length);
void heapify(struct rheap *h);
struct res_node* rheap_remove_max(struct rheap *h);
//void rheap_insert(struct rheap *h, struct knode *node,double value);
void percolateDown(struct rheap *h, int index);
int maximum(double a, int indexa, double b, int indexb);
void percolateUp(struct rheap *h, int index);
void swap(struct rheap *h, int index1, int index2);
struct rheap* initHeap();
