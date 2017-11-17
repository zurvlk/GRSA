#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#define INF DBL_MAX

int nc2(int n) {
    return n * (n - 1) / 2;
}

int isin_array(int *array, int target, int size) {
    for (int i = 1; i <= size; i++) {
        if(array[i] > target) break;
        if(array[i] == target) return 1;
    }
    return 0;
}

int list_isin_array(int *ls, int target) {
    for(int i = 1; i <= ls[0]; i++) {
        if (ls[i] > target) break;
        if (ls[i] == target) return 1;
    }
    return 0;
}

double fmin(double i, double j) {
    return i < j ? i : j;
}

// 2つの配列の値がすべて同一 0
// 2つの配列の値が異なる 1
int cmparray(int *array1, int *array2, int size) {
    int i;
    for (i = 1; i <= size; i++) {
        if(array1[i] != array2[i]) return 1;
    }
    return 0;
}


double theta(double n, double T) {

    if(!function) return fmin(n > 0 ? n : -n, T);
    else return fmin(n * n, T);

}


void cpyarray(int *terget, int *source, int size) {
    int i;

    for (i = 1; i <= size; i++) {
        terget[i] = source[i];
    }
}

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}


double pairwise(double i, double j, double T) {
    return theta(i - j, T);
}

double data(int i, int label) {
    return 1.0 * dabs(label, i);
}

double d_p(int label, int i, int alpha) {
    return dabs(label, (i + alpha));
}

double energy(Graph *G, int *label, int *I, double T) {
    int i;
    double energy = 0;
    //* Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data(I[i], label[i]);
    }
    // */
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T);
    }
    return energy;
}

double r(int i, int k, int grids_edge, double T) {
    return grids_edge * (- 0.5) * (theta(k + 1 - i, T) + theta(i + 1, T));
}

double e_cost(int i, int j, double T) {
    double cost = 0.5 * (theta(i - j + 1, T)- 2 * theta(i - j, T) + theta(i - j - 1, T));
    return cost;
}



//枝を作る関数makeedge(グラフ,高さ,幅)
void set_single_edges(Graph *G, int height, int width) {
    int i, j, edge_count;
    int tail, head, source, sink;

    source = G->n - 1;
    sink = source + 1;
    setSource(G, source);
    setSink(G, sink);
    
    edge_count = 1;
    //点と点の間の枝（横）
    for (i = 1; i < height + 1; i++) {
        for (j = 1; j < width; j++) {
            tail =  (i - 1) * width + j;
            head =  tail + 1;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }
 
    //点と点の間の枝（縦）
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            tail = (i - 1) * width + j;
            head = tail + width;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }
 
    //sourceと点の間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, G->src, i, 0);
       edge_count++;
    }
    
    //点とsinkの間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, i, G->sink, 0);
       edge_count++;
    }
    return;
}

int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta) {
    int i, arraysize;
    for (i = 1; i <= G->n - 2; i++) label_index[i] = 0;
    arraysize = 1;
    for (i = 1; i <= G->n - 2; i++) {
        if (label[i] <= beta && label[i] >= alpha) {
            label_index[arraysize] = i;
            arraysize++;
        }
    }
    return arraysize;
}

double phi (int i, int j, int *ls, double T) {
    double p = 0;
    if(j > i) return 0;
    if(1 < j && j <= i) {
        p = theta(ls[i] - ls[j - 1], T) - theta(ls[i] - ls[j], T) - theta(ls[i - 1] - ls[j - 1], T) + theta(ls[i - 1] - ls[j], T);
        if(i == j) p *= 0.5;
    }
    return p;
}

double nn(int i, int label, int *ls, double T) {
    double p = 0;
    p = pairwise(ls[i], label, T);
    return p;
}

int nnp_4_grsa(int i, int height, int width, int *ls, int *label, double T) {
    int grids_node = height * width;

    double nnp_total = 0;


    if (i > width) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (!list_isin_array(ls, label[i - width])){
            // iの上の点がLs内に含まれない
            nnp_total += pairwise(ls[i], label[i - width], T) ;
        }
    }
    

    if (i <= grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (!list_isin_array(ls, label[i + width])){
            // iの下の点がLs内に含まれない
            nnp_total += pairwise(ls[i], label[i + width], T) ;
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (!list_isin_array(ls, label[i - 1])){
            // iの左の点がLs内に含まれない
            nnp_total += pairwise(ls[i], label[i - 1], T) ;
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (!list_isin_array(ls, label[i + 1])){
            // iの右の点がLs内に含まれない
            nnp_total += pairwise(ls[i], label[i + 1], T) ;
        }
    }
    return nnp_total;
}


// set_edge for rangeswap
void set_edge_for_grsa(Graph *G, int height, int width, int *ls, int *label, int *I, double T) {
    int i, j, k, l, node;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin, sink_label;
    double *min;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * ls[0] + 1;
    sink = source + 1;
    sink_label = ls[0] + 1;

    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        if(list_isin_array(ls, label[i])) {
            G->capa[edge_count] = INF;
        }
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }
    
    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < ls[0]; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);
            if(list_isin_array(ls , label[i])) {
                G->capa[edge_count] = data(I[i],ls[j]) + nnp_4_grsa(i, height, width, ls, label, T);
            }

            if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
            edge_count++;
            tail = head;
        }
    }
    // ik->sink
    i2t_begin = edge_count;
    // rk = r(beta , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (ls[0] - 1), sink, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(list_isin_array(ls , label[i])) {
            G->capa[edge_count] = data(I[i], ls[ls[0]]) + nnp_4_grsa(ls[ls[0]], height, width, ls, label, T);
        }

        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }


    // reverce edge
    // depth
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < ls[0]; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, head, tail, 0);
            if(list_isin_array(ls, label[i])) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }

    

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(list_isin_array(ls, label[t_base]) && list_isin_array(ls, label[h_base])) {
                        // head, tail in label_index
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(list_isin_array(ls, label[t_base]) && list_isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }




    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(list_isin_array(ls, label[t_base]) && list_isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(list_isin_array(ls, label[t_base]) && list_isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = 1; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < ls[0]; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    
    printf("total edge : %d\n", edge_count - 1);
    return;
}

int r_4_rangeswap(int i, int li, int height, int width, int label_size, int *label_index, int *label, int index_size, int alpha, int beta, double T) {
    int grids_node = height * width;

    double r_total = 0;
    double k = dabs(alpha, beta);

    if (i > width) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (isin_array(label_index, i - width, index_size)){
            // iの上の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1, T);
        } else {
            // iの上の点がalpha-beta間に含まれない-= r(li, k, 1);
            r_total += theta(li + alpha - label[i - width], T);
        }
    }
    

    if (i <= grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (isin_array(label_index, i + width, index_size)){
            // iの下の点がalpha-beta間に含まない
            r_total -= r(li, k, 1, T);
        } else {
            // iの下の点がalpha-beta間に含まれない
            r_total += theta(li + alpha -  label[i + width], T);
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (isin_array(label_index, i - 1, index_size)){
            // iの左の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1, T);
        } else {
            // iの左の点がalpha-beta間に含まれない
            r_total += theta(li + alpha -  label[i - 1], T);
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (isin_array(label_index, i + 1, index_size)){
            // iの右の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1, T);
        } else {
            // iの右の点がalpha-beta間に含まれない
            r_total += theta(li + alpha - label[i + 1], T);
        }
    }
    return r_total;
}

int near_node_counts(int i, int height, int width, int label_size, int *label_index, int index_size) {
    int nn = 0;
    int grids_node = height * width;

    if (i > width) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (isin_array(label_index, i - width, index_size)){
            // iの上の点がalpha-beta間に含まれる
            nn++;
        }
    }
    

    if (i <= grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (isin_array(label_index, i + width, index_size)){
            // iの下の点がalpha-beta間に含まない
            nn++;
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (isin_array(label_index, i - 1, index_size)){
            // iの左の点がalpha-beta間に含まれる
            nn++;
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (isin_array(label_index, i + 1, index_size)){
            // iの右の点がalpha-beta間に含まれる
            nn++;
        }
    }
    return nn;

}

// set_edge for rangeswap
void set_edge(Graph *G, int height, int width, int alpha, int beta, int label_size, int *label, int *label_index, int size, int *I, double T) {
    int i, j, k, l, node;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin, range_size, sink_label;
    double *min, r_total;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;
    range_size = beta - alpha;

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * range_size + 1;
    sink = source + 1;
    sink_label = range_size + 1;

    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        if(isin_array(label_index, i, size)) {
            r_total = r_4_rangeswap(i, 0, height, width, label_size, label_index, label, size, alpha, beta, T);
            // G->capa[edge_count] = data(I, i, alpha) - r0 + e_cost(alpha, alpha + 1);
            G->capa[edge_count] = d_p(label[i], 0, alpha) + r_total + near_node_counts(i, height, width, label_size, label_index, size) * e_cost(0, 1, T);
        }
        
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }
    // printf("%lf\n", G->capa[edge_count - 1]);
    // source->i2~ik
    for (i = grids_node + 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, source, i, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = near_node_counts(node, height, width, label_size, label_index, size) * e_cost(0, (i - 1) / grids_node + 1, T);
        }

        edge_count++;
    }
    // printf("%d %d \n",i,  i % grids_node == 0 ? grids_node : i % grids_node);
    // printf("::%d\n", range_size); 

    // i1~ik-1->sink
    for (i = 1; i <= (range_size - 1) * grids_node; i++) {
        setEdge(G, edge_count, i, sink, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = near_node_counts(node, height, width, label_size, label_index, size) * e_cost((i - 1) / grids_node + 1, sink_label, T);
        }

        edge_count++;
    }
    // printf("%lf\n", G->capa[edge_count - 1]);
    // ik->sink
    i2t_begin = edge_count;
    // rk = r(beta , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (range_size - 1), sink, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            r_total = r_4_rangeswap(i, beta, height, width, label_size, label_index, label, size, alpha, beta, T);
            G->capa[edge_count] = d_p(I[i], range_size , alpha) + r_total + near_node_counts(i, height, width, label_size, label_index, size) * e_cost(range_size, sink_label, T);
        }

        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }

    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < range_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);

            if(isin_array(label_index, i , size)) {
                r_total = r_4_rangeswap(i, j, height, width, label_size, label_index, label, size, alpha, beta, T);
                G->capa[edge_count] = d_p(I[i], j, alpha) + r_total;
            }

            if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        // head, tail in label_index
                        G->capa[edge_count] = e_cost(k + 1, l + 1, T);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1, T);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }


    // reverce edge
    // depth
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < range_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, head, tail, 0);
            if(isin_array(label_index, i, size)) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1, T);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1, T);
                    }
                    edge_count++;
                }
            }
        }
    }


    // each node 2 source 
    for (i = 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, i, source, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1, T);
        }
        edge_count++;
    }
     // each node 2 sink
    for (i = 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, sink, i, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, sink_label, T);
        }
        edge_count++;
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = 1; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < range_size; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    // printf("total edge : %d\n", edge_count - 1);
    return;
}