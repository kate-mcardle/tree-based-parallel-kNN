/*************************************************
 *
 * Author: Nazneen Rajani
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include "kdtree.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <omp.h>
#include<iostream>
using namespace std;

unsigned int get_msec(void)
{
	static struct timeval timeval, first_timeval;

	gettimeofday(&timeval, 0);

	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}

int read_data_from_libsvm(double **data, char *filename, int N, int D) {
	int idx; int dim;
	for (idx = 0; idx < N; ++idx) {
		for (dim = 0; dim < D; ++dim) {
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
  return 1;
}

int main(int argc, char **argv)
{
    int i, j;
    struct kdres *set;
    struct kdtree *kd;
    unsigned int msec, start;
    double end, start1;
    int n = atoi(argv[1]);                 // # data points
    int D = atoi(argv[2]);                  // data dimension
    char *fname = argv[3]; // data file
    double **X = NULL;
    X= (double**)malloc(sizeof(double*)*n);
    for (i=0; i<n; i++) {
        X[i] = NULL;
        X[i] = (double*)malloc(sizeof(double)*D);
    }
    read_data_from_libsvm(X,fname,n,D); 
	kd = kd_create(D);
    //printf("here now create\n");

    double pos[D];
	start = get_msec();
	for(i=0; i<n; i++) {
        	for (j=0; j<D; j++) {
            		pos[j] = X[i][j];
        	}
        	kd_insert(kd, pos, i);
	}
	msec = get_msec() - start;
     printf("%.3f sec for creating tree \n", (float)msec / 1000.0);
    int n_thread = atoi(argv[4]);	
    omp_set_num_threads(n_thread);
    //start = get_msec();
    start1 = omp_get_wtime();
    //printf("%f sec for querying\n", start1);
    #pragma omp parallel for default(shared) private(i)
    for (i=0; i<n; i++) {
        kd_nearest_n(kd, X[i],8);
    }
    end = omp_get_wtime()-start1;
    //msec = get_msec() - start;
    cout<<"Querying time: "<<end<<endl;
	kd_free(kd);
	return 0;
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
    
    //printf("done\n");
    
    return 1;
}




