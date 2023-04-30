//
// Created by kaixin on 3/25/23.
//

#ifndef DESCOL_EBBKC_H
#define DESCOL_EBBKC_H
#include "def.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class EBBkC_Graph_t {

public:
    int v_size = 0;
    int e_size = 0;
    Edge_t* edges;

    int* sub_e_size;     // number of edges in the sub-branch
    int* sub_v_size;     // number of nodes in the sub-branch
    int** sub_v;  // nodes (id) in the sub-branch
    int** sub_e;  // edges (id) in the sub-branch

//    int *sub_v_size = nullptr;
//    int *sub_e_size = nullptr;
//    int **sub_v = nullptr;
//    int **sub_e = nullptr;

//    int color_num = 0;
    vector<int> new2old;
//    int **d = nullptr;
//    int **od = nullptr;
//    int *cd = nullptr;
//    int *adj = nullptr;

    vector<int> rank;
//    int *rank = nullptr;
    int* col;
//    int *color = nullptr;
    int* lab;
    int** out_deg;
    int** deg;
    int**  G_adj;       // G_adj[l][node (id)] = a set of nodes (id) -- neighbors
    int** DAG_adj;      // DAG_adj[l][node (id)] = a set of nodes (id) -- out-neighbors
    bool** used;

    // edge-oriented branch
    int truss_num = 0;
    HashMap_t edge2id;

    int** T;   // T[edge (id)] = a set of nodes (id)
    int* T_size;
    int** C;   // C[edge (id)] = a set of edges (id)
    int* C_size;

    // early-stop
    int plex_sub = 0;
    int plex_size = 0;
    int plex_p_size = 0;
    int plex_f_size = 0;
    int *plex_new2old = nullptr;
    int *plex_p = nullptr;
    int *plex_f = nullptr;
    int *plex_d_num = nullptr;
    int **plex_d = nullptr;
    int *plex_lev = nullptr;
    int *plex_loc = nullptr;

    EBBkC_Graph_t();
    ~EBBkC_Graph_t();

    void read_edges_from_file(const char* file_name);
    void truss_decompose(const char* w_file_name);
    void read_ordered_edges_from_file(const char* file_name);
    void build_from_G();
    void color_DAG();
    void build_from_DAG(int height);

    void basic_clique_list(int l, unsigned long long *cliques);
    void basic_clique_list_plus(int l, unsigned long long *cliques);
    void improved_clique_list(int, unsigned long long *cliques);
    void improved_clique_list_plus(int l, unsigned long long *cliques);
    void build();

    bool list_clique_in_plex(int l, unsigned long long* cliques);
    void list_in_plex(int p, int q, unsigned long long* cliques);
};


class EBBkC_t {
public:
    EBBkC_t();

    double truss_order(const char* r_file_name, const char* w_file_name);
    double list_k_clique(const char* file_name, int type);
};


#endif //DESCOL_EBBKC_H
