//
// Created by kaixin on 3/25/23.
//

#ifndef DESCOL_VBBkC_H
#define DESCOL_VBBkC_H
#include "def.h"

class VBBkC_Graph_t {
public:
    int n = 0;
    int m = 0;
    int depth = 0;
    int core_num = 0;
    int color_num = 0;
    int *new2old = nullptr;
    int *ns = nullptr;
    int **d = nullptr;
    int **od = nullptr;
    int *cd = nullptr;
    int *adj = nullptr;
    int *rank = nullptr;
    int *color = nullptr;
    int *lab = nullptr;
    int **sub = nullptr;
    Edge_t *edges = nullptr;

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

    VBBkC_Graph_t();
    ~VBBkC_Graph_t();

    void read_edges_from_file(const char* file_name);
    void core_DAG();
    void color_DAG();
    void build_from_DAG(int height);
    void vertex_oriented_branch(VBBkC_Graph_t* sg, int node);
    void vertex_oriented_twice_branch(VBBkC_Graph_t* sg, int edge);
    void list_clique(int l, unsigned long long* cliques);
    void list_clique_plus(int l, unsigned long long* cliques); // early-stop
    bool list_clique_in_plex(int l, unsigned long long* cliques);
    void list_in_plex(int p, int q, unsigned long long* cliques);
};

class VBBkC_t {
public:
    VBBkC_t();

    double list_k_clique(const char* file_name, int type);
};

#endif //DESCOL_VBBkC_H
