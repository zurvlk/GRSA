#include "graph.h"
#include "bmp.h"

int function;

double theta(double n, double T);
double p(int *label, int height, int width);
double energy(Graph *G, int *label, int *I, double T);
double pairwise(double i, double j, double T);
double data(int i, int label);
void set_edge_for_grsa(Graph *G, int height, int width, int *ls, int *label, int *I, double T);

int list_isin_array(int *ls, int target);
void set_edge(Graph *G, int height, int width, int alpha, int beta, int label_size, int *label, int *label_index, int size, int *I, double T);
int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta);
void set_single_edges(Graph *G, int height, int width);
int isin_array(int *array, int target, int size);
int cmparray(int *array1, int *array2, int size);
void cpyarray(int *terget, int *source, int size);