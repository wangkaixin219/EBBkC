#include <set>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <cassert>
#include <chrono>

#include "set_operation.h"
#include "edge_oriented.h"
#include "truss/util/graph/graph.h"
#include "truss/util/log/log.h"
#include "truss/util/timer.h"
#include "truss/decompose/parallel_all_edge_cnc.h"
#include "truss/util/reordering/reorder_utils.h"
#include "truss/decompose/iter_helper.h"

extern const int K, L;
extern unsigned long long N;

EBBkC_Graph_t::EBBkC_Graph_t() = default;

EBBkC_Graph_t::~EBBkC_Graph_t() {
    int i, k = K, node_size = v_size, link_size = e_size;

    if (is_sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = truss_num * (truss_num - 1) / 2;
    }

    delete [] edges;

    if (T) {
        for (i = 0; i < link_size; i++) delete [] T[i];
        delete [] T;
    }

    if (C) {
        for (i = 0; i < link_size; i++) delete [] C[i];
        delete [] C;
    }

    delete [] T_size;

    delete [] C_size;

    if (sub_v) {
        for (i = 0; i <= k; i++) delete [] sub_v[i];
        delete [] sub_v;
    }

    if (sub_e) {
        for (i = 0; i <= k; i++) delete [] sub_e[i];
        delete [] sub_e;
    }

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    if (DAG_deg) {
        for (i = 0; i <= k; i++) delete [] DAG_deg[i];
        delete [] DAG_deg;

    }

    if (G_deg) {
        for (i = 0; i <= k; i++) delete [] G_deg[i];
        delete [] G_deg;
    }


    delete [] col;

    if (DAG_adj) {
        for (i = 0; i < node_size; i++) delete [] DAG_adj[i];
        delete [] DAG_adj;
    }

    if (G_adj) {
        for (i = 0; i < node_size; i++) delete [] G_adj[i];
        delete [] G_adj;
    }

    if (used) {
        for (i = 0; i <= k; i++) delete [] used[i];
        delete [] used;
    }


    delete [] v_lab;

    if (out_v_size) {
        for (i = 0; i <= k; i++) delete [] out_v_size[i];
        delete [] out_v_size;
    }

    if (out_e_size) {
        for (i = 0; i <= k; i++) delete [] out_e_size[i];
        delete [] out_e_size;
    }

    delete [] F;

    delete [] P;

    delete [] lack_size;

    if (lack) {
        for (i = 0; i < node_size; i++) delete [] lack[i];
        delete [] lack;
    }

    delete [] lev;

    delete [] loc;
}

//void EBBkC_Graph_t::read_edges_from_file(const char *file_name) {
//    FILE *fp;
//    if ((fp = fopen(file_name, "r")) == nullptr) {
//        printf("Cannot open file %s.\n", file_name);
//        exit(0);
//    }
//
//    Edge_t e;
//    int u, v, i;
//    int *old2new = new int [N_NODES];
//    for (i = 0; i < N_NODES; i++) old2new[i] = -1;
//
//    e_size = 0;
//    edges = new Edge_t [N_EDGES];
//
//    while (fscanf(fp, "%d %d%*[^\n]%*c", &u, &v) == 2) {
//
//        if (u > N_NODES || v > N_NODES) {
//            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
//            exit(0);
//        }
//        if (old2new[u] == -1) {
//            old2new[u] = (int) new2old.size();
//            new2old.push_back(u);
//        }
//        if (old2new[v] == -1) {
//            old2new[v] = (int) new2old.size();
//            new2old.push_back(v);
//        }
//
//        e = Edge_t(old2new[u], old2new[v], false);
//        edges[e_size++] = e;
//    }
//
//    v_size = (int) new2old.size();
//
//    fclose(fp);
//
//    delete [] old2new;
//}
//
//void EBBkC_Graph_t::read_ordered_edges_from_file(const char *file_name) {
//    FILE *fp;
//    if ((fp = fopen(file_name, "r")) == nullptr) {
//        printf("Cannot open file %s.\n", file_name);
//        exit(0);
//    }
//
//    Edge_t e, e_;
//    int u, v, w, i, j, k, idx, edge_rank, edge_sub_size;
//    int *old2new = new int [N_NODES];
//    for (i = 0; i < N_NODES; i++) old2new[i] = -1;
//    vector<int> t_;
//    vector<vector<int>> T_;
//
//    e_size = 0;
//    edges = new Edge_t [N_EDGES];
//
//    while (fscanf(fp, "%d %d %d %d %d", &u, &v, &k, &edge_rank, &edge_sub_size) == 5) {
//
//        if (k <= K) {
//            fscanf(fp, "%*[^\n]%*c");
//            continue;
//        }
//        truss_num = truss_num > k ? truss_num : k;
//
//        if (u > N_NODES || v > N_NODES) {
//            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
//            exit(0);
//        }
//        if (old2new[u] == -1) {
//            old2new[u] = (int) new2old.size();
//            new2old.push_back(u);
//        }
//        if (old2new[v] == -1) {
//            old2new[v] = (int) new2old.size();
//            new2old.push_back(v);
//        }
//
//        e = Edge_t(old2new[u], old2new[v], false);
//        rank.push_back(edge_rank);
//        edge2id.insert(e, e_size);
//        edges[e_size++] = e;
//
//        for (j = 0; j < edge_sub_size; j++) {
//            fscanf(fp, "%d", &w);
//            if (old2new[w] == -1) {
//                old2new[w] = (int) new2old.size();
//                new2old.push_back(w);
//            }
//            t_.push_back(old2new[w]);
//        }
//        T_.push_back(t_);
//        t_.clear();
//    }
//
//    v_size = (int) new2old.size();
//
//    T = new int* [e_size];
//    for (i = 0; i < e_size; i++) T[i] = new int [truss_num + 1];
//    T_size = new int [e_size];
//
//    C = new int* [e_size];
//    for (i = 0; i < e_size; i++) C[i] = new int [T_[i].size() * (T_[i].size() - 1) / 2];
//    C_size = new int [e_size];
//
//    for (i = 0; i < e_size; i++) {
//        T_size[i] = C_size[i] = 0;
//
//        for (j = 0; j < T_[i].size(); j++) T[i][T_size[i]++] = T_[i][j];
//
//        for (j = 0; j < T_size[i]; j++) {
//            for (k = j + 1; k < T_size[i]; k++) {
//
//                e = Edge_t(T[i][j], T[i][k], false);
//
//                if ((idx = edge2id.exist(e)) != -1  && rank[idx] > rank[i]) {
//                    C[i][C_size[i]++] = idx;
//                }
//            }
//        }
//    }
//
//    printf("|V| = %d, |E| = %d\n", v_size, e_size);
//    printf("Truss number = %d\n", truss_num - 2);
//
//    fclose(fp);
//    delete [] old2new;
//}
//
//void EBBkC_Graph_t::truss_decompose(const char* w_file_name) {
//
//    int i, j, k, s, t, w, end, edge_seq = 1, sw, wt;
//
//    auto *_d = new int [v_size]();
//    auto *_cd = new int [v_size + 1];
//    auto *_adj = new int [2 * e_size];
//
//    for (i = 0; i < e_size; i++) {
//        _d[edges[i].s]++;
//        _d[edges[i].t]++;
//    }
//    _cd[0] = 0;
//    for (i = 1; i < v_size + 1; i++) {
//        _cd[i] = _cd[i - 1] + _d[i - 1];
//        _d[i - 1] = 0;
//    }
//    for (i = 0; i < e_size; i++) {
//        _adj[_cd[edges[i].s] + _d[edges[i].s]++] = edges[i].t;
//        _adj[_cd[edges[i].t] + _d[edges[i].t]++] = edges[i].s;
//    }
//
//    KeyVal_t kv;
//    Heap_t heap;
//    Edge_t  e;
//
//    auto *sup = new int [e_size]();
//    auto *h_table = new bool [v_size]();
//    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_id;
//
//    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_truss;
//    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_rank;
//    unordered_map<Edge_t, vector<int>, Edge_t::Hash_Edge_t> edge_sub;
//
//    for (i = 0; i < e_size; i++) {
//
//        edge_id[edges[i]] = i;
//
//        s = _d[edges[i].s] < _d[edges[i].t] ? edges[i].s : edges[i].t;
//        t = _d[edges[i].s] < _d[edges[i].t] ? edges[i].t : edges[i].s;
//
//        for (j = _cd[s]; j < _cd[s + 1]; j++) {
//            w = _adj[j];
//            h_table[w] = true;
//        }
//
//        for (j = _cd[t]; j < _cd[t + 1]; j++) {
//            w = _adj[j];
//            if (h_table[w])
//                sup[i]++;
//        }
//
//        for (j = _cd[s]; j < _cd[s + 1]; j++) {
//            w = _adj[j];
//            h_table[w] = false;
//        }
//    }
//
//    heap.make_heap(sup, e_size);
//
//    for (k = 3; !heap.empty(); k++) {
//
//        while (!heap.empty() && heap.min_element().val < k - 2) {
//            kv = heap.pop();
//            e = edges[kv.key];
//            edge_truss[e] = k;
//            edge_rank[e] = edge_seq++;
//
//            s = _d[e.s] < _d[e.t] ? e.s : e.t;
//            t = _d[e.s] < _d[e.t] ? e.t : e.s;
//
//            end = _cd[s] + _d[s];
//            for (j = _cd[s]; j < end; j++) {
//                w = _adj[j];
//                h_table[w] = true;
//            }
//
//            end = _cd[t] + _d[t];
//            for (j = _cd[t]; j < end; j++) {
//                w = _adj[j];
//
//                if (w == s) {
//                    _adj[j--] = _adj[--end];
//                    _d[t]--;
//                }
//
//                if (h_table[w]) {
//                    edge_sub[e].push_back(w);
//                    heap.update(edge_id[Edge_t(s, w, false)]);
//                    heap.update(edge_id[Edge_t(w, t, false)]);
//                }
//            }
//
//            end = _cd[s] + _d[s];
//            for (j = _cd[s]; j < end; j++) {
//                w = _adj[j];
//                h_table[w] = false;
//
//                if (w == t) {
//                    _adj[j--] = _adj[--end];
//                    _d[s]--;
//                }
//            }
//        }
//    }
//
//    heap.release_heap();
//
//    FILE *fp = fopen(w_file_name, "w");
//
//    for (i = 0; i < e_size; i++) {
//        e = edges[i];
//        fprintf(fp, "%d %d %d %d %ld ", new2old[e.s], new2old[e.t], edge_truss[e], edge_rank[e], edge_sub[e].size());
//        for (j = 0; j < edge_sub[e].size(); j++) fprintf(fp, "%d ", new2old[edge_sub[e][j]]);
//        fprintf(fp, "\n");
//    }
//
//    fclose(fp);
//
//    delete [] _d;
//    delete [] _cd;
//    delete [] _adj;
//    delete [] h_table;
//}


void EBBkC_Graph_t::build(bool sub) {
    int i, k = K, node_size = v_size, link_size = e_size;

    is_sub = sub;

    if (sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = (truss_num) * (truss_num - 1) / 2;
    }

    sub_v = new int* [k + 1];

    sub_e = new int* [k + 1];

    sub_e_size = new int [k + 1];

    sub_v_size = new int [k + 1];

    for (i = 0; i < k; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[k] = new int [node_size];

    for (i = 0; i < k; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[k] = new int [link_size];

    sub_v_size[k] = 0;

    sub_e_size[k] = 0;

    lab = new int [node_size];
    for (i = 0; i < node_size; i++) lab[i] = k;

    DAG_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) DAG_deg[i] = new int [node_size];

    G_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) G_deg[i] = new int [node_size];

    col = new int [node_size];

    DAG_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [k + 1];
    for (i = 0; i <= k; i++) used[i] = new bool [node_size + 1]();

    v_lab = new int [node_size];
    for (i = 0; i < node_size; i++) v_lab[i] = k;

    out_v_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_v_size[i] = new int [link_size];

    out_e_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_e_size[i] = new int [link_size];

    F = new int [truss_num + 1];

    P = new int [truss_num + 1];

    lack_size = new int [node_size];

    lack = new int* [node_size];
    for (i = 0; i < node_size; i++) lack[i] = new int [L + 1];

    lev = new int [node_size]();

    loc = new int [node_size];
}

//
//void EBBkC_Graph_t::EBBkC(int l, unsigned long long *cliques) {
//    int i, j, k, u, e, e_, _e, end;
//
//    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;
//
//    if (l == K) {
//        if (K == 3) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                for (j = 0; j < T_size[e]; j++) {
//                    (*cliques)++;
//                }
//            }
//        }
//        else if (K == 4) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                for (j = 0; j < C_size[e]; j++) {
//                    (*cliques)++;
//                }
//            }
//        }
//        else {
//
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                sort(T[e], T[e] + T_size[e]);
//                sort(C[e], C[e] + C_size[e]);
//            }
//
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;
//
//                sub_v_size[l - 2] = 0;
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
//                }
//
//                sub_e_size[l - 2] = 0;
//
//                for (j = 0; j < C_size[e]; j++) {
//                    e_ = C[e][j];
//                    sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                }
//
//                EBBkC(l - 2, cliques);
//
//            }
//        }
//
//        return;
//    }
//
//    for (i = 0; i < sub_e_size[l]; i++) {
//        e = sub_e[l][i];
//
//        if (l == 3) {
//            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
//            for (j = 0; j < sub_v_size[l - 2]; j++) {
//                (*cliques)++;
//            }
//        }
//
//        else if (l == 4) {
//             sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
//             for (j = 0; j < sub_e_size[l - 2]; j++) {
//                (*cliques)++;
//             }
//        }
//
//        else {
//            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
//            sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
//            EBBkC(l - 2, cliques);
//        }
//    }
//}
//
//
//void EBBkC_Graph_t::EBBkC_plus(int l, unsigned long long *cliques) {
//    int c, i, j, k, e, e_, u, v, w, s, t, end, dist;
//
//    if (sub_v_size[l] < l) return;
//
//    if (l == K) {
//        if (K == 3) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                 e = sub_e[l][i];
//
//                 for (j = 0; j < T_size[e]; j++) {
//                    (*cliques)++;
//                 }
//            }
//        }
//        else if (K == 4) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                for (j = 0; j < C_size[e]; j++) {
//                    (*cliques)++;
//                }
//            }
//        }
//        else {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    col[u] = 0;
//                    DAG_deg[0][u] = 0;
//                    G_deg[l - 2][u] = 0;
//                }
//
//                for (j = 0; j < C_size[e]; j++) {
//                    e_ = C[e][j];
//                    s = edges[e_].s;
//                    t = edges[e_].t;
//                    G_adj[s][G_deg[l - 2][s]++] = t;
//                    G_adj[t][G_deg[l - 2][t]++] = s;
//                }
//
//                auto *list = new KeyVal_t [truss_num + 1];
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    list[j].key = u;
//                    list[j].val = G_deg[l - 2][u];
//                }
//                sort(list, list + T_size[e]);
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = list[j].key;
//                    for (k = 0; k < G_deg[l - 2][u]; k++) {
//                        v = G_adj[u][k];
//                        used[K][col[v]] = true;
//                    }
//                    for (c = 1; used[K][c]; c++) ;
//                    col[u] = c;
//                    for (k = 0; k < G_deg[l - 2][u]; k++) {
//                        v = G_adj[u][k];
//                        used[K][col[v]] = false;
//                    }
//                }
//                delete [] list;
//
//                sub_v_size[l - 2] = 0;
//                dist = 0;
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
//                    if (!used[K][col[u]]) {
//                        used[K][col[u]] = true;
//                        dist++;
//                    }
//                }
//
//                if (dist >= l - 2) {
//                    sort(sub_v[l - 2], sub_v[l - 2] + sub_v_size[l - 2]);
//
//                    sub_e_size[l - 2] = 0;
//                    for (j = 0; j < C_size[e]; j++) {
//                        e_ = C[e][j];
//                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                        s = edges[e_].s;
//                        t = edges[e_].t;
//                        edges[e_].s = (col[s] > col[t]) ? s : t;
//                        edges[e_].t = (col[s] > col[t]) ? t : s;
//                        s = edges[e_].s;
//                        t = edges[e_].t;
//                        DAG_adj[s][DAG_deg[0][s]++] = t;
//                    }
//
//                    for (j = 0; j < T_size[i]; j++) {
//                        u = T[e][j];
//                        // sorted array for SIMD usage, DAG_adj can be considered as const.
//                        sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]);
//                    }
//
//                    EBBkC_plus(l - 2, cliques);
//                }
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    used[K][col[u]] = false;
//                }
//            }
//        }
//
//        return;
//    }
//
//    if (l == 1) {
//        for (i = 0; i < sub_v_size[l]; i++) {
//            (*cliques)++;
//        }
//
//        return;
//    }
//
//    for (i = 0; i < sub_v_size[l]; i++) {
//        u = sub_v[l][i];
//
//        if (col[u] < l) continue;
//
//        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);
//
//        if (l == 2) {
//            for (j = 0; j < sub_v_size[l - 1]; j++) {
//                (*cliques)++;
//            }
//        }
//
//        else {
//            if (sub_v_size[l - 1] >= l - 1) {
//                EBBkC_plus(l - 1, cliques);
//            }
//        }
//    }
//}

//
//void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
//    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;
//
//    if (l == K) {
//
//        for (i = 0; i < v_size; i++) sub_v[K][sub_v_size[K]++] = i;
//
//        for (i = 0; i < e_size; i++) sub_e[K][sub_e_size[K]++] = i;
//
//        if (K == 3) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//                for (j = 0; j < T_size[e]; j++) {
//                    (*cliques)++;
//                }
//            }
//        }
//        else if (K == 4) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//                for (j = 0; j < C_size[e]; j++) {
//                    (*cliques)++;
//                }
//            }
//        }
//        else {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    col[u] = 0;
//                    DAG_deg[l - 2][u] = 0;
//                    G_deg[l - 2][u] = 0;
//                }
//
//                for (j = 0; j < C_size[e]; j++) {
//                    e_ = C[e][j];
//                    s = edges[e_].s;
//                    t = edges[e_].t;
//                    G_adj[s][G_deg[l - 2][s]++] = t;
//                    G_adj[t][G_deg[l - 2][t]++] = s;
//                }
//
//                auto *list = new KeyVal_t [truss_num + 1];
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    list[j].key = u;
//                    list[j].val = G_deg[l - 2][u];
//                }
//                sort(list, list + T_size[e]);
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = list[j].key;
//                    for (k = 0; k < G_deg[l - 2][u]; k++) {
//                        v = G_adj[u][k];
//                        used[K][col[v]] = true;
//                    }
//                    for (c = 1; used[K][c]; c++) ;
//                    col[u] = c;
//                    for (k = 0; k < G_deg[l - 2][u]; k++) {
//                        v = G_adj[u][k];
//                        used[K][col[v]] = false;
//                    }
//                }
//                delete [] list;
//
//                sub_v_size[l - 2] = 0;
//                dist = 0;
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
//                    if (!used[K][col[u]]) {
//                        used[K][col[u]] = true;
//                        dist++;
//                    }
//                }
//
//                if (dist >= l - 2) {
//                    sub_e_size[l - 2] = 0;
//                    for (j = 0; j < C_size[e]; j++) {
//                        e_ = C[e][j];
//                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                        s = edges[e_].s;
//                        t = edges[e_].t;
//                        edges[e_].s = (col[s] > col[t]) ? s : t;
//                        edges[e_].t = (col[s] > col[t]) ? t : s;
//                        s = edges[e_].s;
//                        t = edges[e_].t;
//
//                        DAG_adj[s][DAG_deg[l - 2][s]++] = t;
//                    }
//
//                    EBBkC_plus_plus(l - 2, cliques);
//                }
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    used[K][col[u]] = false;
//                }
//            }
//        }
//
//        return;
//    }
//
//    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;
//
//    if (l == 2) {
//        for (i = 0; i < sub_v_size[l]; i++) {
//            u = sub_v[l][i];
//
//            for (j = 0; j < DAG_deg[l][u]; j++) {
//                (*cliques)++;
//            }
//        }
//
//        return;
//    }
//
//    if (l == 3) {
//
//        for (i = 0; i < sub_v_size[l]; i++) {
//            u = sub_v[l][i];
//
//            if (col[u] < l) continue;
//
//            for (j = 0; j < DAG_deg[l][u]; j++) {
//                v = DAG_adj[u][j];
//                lab[v] = l - 1;
//            }
//
//            for (j = 0; j < DAG_deg[l][u]; j++) {
//                v = DAG_adj[u][j];
//
//                if (col[v] < l - 1) continue;
//
//                for (k = 0; k < DAG_deg[l][v]; k++) {
//                    w = DAG_adj[v][k];
//                    if (lab[w] == l - 1) (*cliques)++;
//                }
//            }
//
//            for (j = 0; j < DAG_deg[l][u]; j++) {
//                v = DAG_adj[u][j];
//                lab[v] = l;
//            }
//        }
//
//        return;
//    }
//
//    if (can_terminate(l, cliques)) {
//        return;
//    }
//
//    for (i = 0; i < sub_v_size[l]; i++) {
//        u = sub_v[l][i];
//
//        if (col[u] < l) continue;
//
//        sub_v_size[l - 1] = 0;
//        dist = 0;
//
//        for (j = 0; j < DAG_deg[l][u]; j++) {
//            v = DAG_adj[u][j];
//            lab[v] = l - 1;
//            sub_v[l - 1][sub_v_size[l - 1]++] = v;
//            DAG_deg[l - 1][v] = 0;
//            G_deg[l - 1][v] = 0;
//
//            if (!used[l][col[v]]) {
//                used[l][col[v]] = true;
//                dist++;
//            }
//        }
//
//        if (dist >= l - 1) {
//
//            sub_e_size[l - 1] = 0;
//            for (j = 0; j < sub_v_size[l - 1]; j++) {
//                v = sub_v[l - 1][j];
//
//                end = DAG_deg[l][v];
//                for (k = 0; k < end; k++) {
//                    w = DAG_adj[v][k];
//                    if (lab[w] == l - 1) {
//                        DAG_deg[l - 1][v]++;
//                        sub_e_size[l - 1]++;
//
//                        // just for early-termination
//                        G_deg[l - 1][v]++;
//                        G_deg[l - 1][w]++;
//
//                    } else {
//                        DAG_adj[v][k--] = DAG_adj[v][--end];
//                        DAG_adj[v][end] = w;
//                    }
//                }
//            }
//
//            EBBkC_plus_plus(l - 1, cliques);
//        }
//
//        for (j = 0; j < sub_v_size[l - 1]; j++) {
//            v = sub_v[l - 1][j];
//            lab[v] = l;
//            used[l][col[v]] = false;
//        }
//    }
//}

void EBBkC_Graph_t::branch(int e, EBBkC_Graph_t* g) {
    int c, i, j, k, p, e_, u, v, w, s, t, end, dist, l = K;
    int *old2new = new int[v_size];

    g->v_size = 0;
    g->e_size = 0;
    g->edges = new Edge_t[C_size[e] + 1];
    g->new2old = vector<int>(T_size[e]);

    for (i = 0; i < T_size[e]; i++) {
        u = T[e][i];
        old2new[u] = g->v_size;
        g->new2old[g->v_size++] = u;
    }

    for (i = 0; i < C_size[e]; i++) {
        e_ = C[e][i];
        s = edges[e_].s;
        t = edges[e_].t;
        g->edges[g->e_size].s = old2new[s];
        g->edges[g->e_size++].t = old2new[t];
    }

    delete [] old2new;

    if (l == 3 || l == 4) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    if (g->v_size < l - 2 || g->e_size < (l - 2) * (l - 3) / 2) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    for (j = 0; j < g->v_size; j++) {
        g->col[j] = 0;
        g->DAG_deg[l - 2][j] = 0;
        g->G_deg[l - 2][j] = 0;
    }

    for (j = 0; j < g->e_size; j++) {
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->G_adj[s][g->G_deg[l - 2][s]++] = t;
        g->G_adj[t][g->G_deg[l - 2][t]++] = s;
    }

    auto *list = new KeyVal_t[truss_num + 1];
    for (j = 0; j < g->v_size; j++) {
        list[j].key = j;
        list[j].val = g->G_deg[l - 2][j];
    }
    sort(list, list + g->v_size);

    for (j = 0; j < g->v_size; j++) {
        u = list[j].key;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = true;
        }
        for (c = 1; g->used[l - 2][c]; c++);
        g->col[u] = c;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = false;
        }
    }
    delete[] list;

    g->sub_v_size[l - 2] = 0;

    for (j = 0; j < g->v_size; j++) {
        g->sub_v[l - 2][g->sub_v_size[l - 2]++] = j;
    }

    g->sub_e_size[l - 2] = 0;
    for (j = 0; j < g->e_size; j++) {
        g->sub_e[l - 2][g->sub_e_size[l - 2]++] = j;
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->edges[j].s = (g->col[s] > g->col[t]) ? s : t;
        g->edges[j].t = (g->col[s] > g->col[t]) ? t : s;
        s = g->edges[j].s;
        t = g->edges[j].t;

        g->DAG_adj[s][g->DAG_deg[l - 2][s]++] = t;
    }

    return;
}

void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (K == 3) {
        (*cliques) += sub_v_size[l];
        return;
    }

    if (K == 4) {
        (*cliques) += sub_e_size[l];
        return;
    }

    if (l == 2) {
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];

            for (j = 0; j < DAG_deg[l][u]; j++) {
                (*cliques)++;
            }
        }

        return;
    }

    if (l == 3) {

        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];

            if (col[u] < l) continue;

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l - 1;
            }

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];

                if (col[v] < l - 1) continue;

                for (k = 0; k < DAG_deg[l][v]; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) (*cliques)++;
                }
            }

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l;
            }
        }

        return;
    }

    if (can_terminate(l, cliques)) {
        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = 0;
        dist = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }

        if (dist >= l - 1) {

            sub_e_size[l - 1] = 0;
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                v = sub_v[l - 1][j];

                end = DAG_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        DAG_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;

                        // just for early-termination
                        G_deg[l - 1][v]++;
                        G_deg[l - 1][w]++;

                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }

            EBBkC_plus_plus(l - 1, cliques);
        }

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }
    }
}

void EBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long *cliques) {
    if (picked == k) {
        (*cliques)++;
        return;
    }

    for (int i = start; i < list_size; i++) {
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques);
    }
}

void EBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (F_size < q) return;

    if (p == 0) {
        if (q > F_size - q)
            EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        else
            EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
        return;
    }

    int i, j, u, v, vis = 0;

    for (i = start; i < P_size && P_act >= p; i++) {
        u = P[i];

        if (lev[u]) continue;

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == 0) {
                lev[v] = p;
                P_act--;
            }
        }

        list_in_plex(i + 1, p - 1, q, cliques);

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == p) {
                lev[v] = 0;
                P_act++;
            }
        }

        P_act--;
        vis++;
    }
    P_act += vis;
}

bool EBBkC_Graph_t::can_terminate(int l, unsigned long long *cliques) {
    int i, j, k, u, v, end, p_;

    if (sub_e_size[l] < sub_v_size[l] * (sub_v_size[l] - L) / 2) return false;

    if (sub_e_size[l] == sub_v_size[l] * (sub_v_size[l] - 1) / 2) {
        if (l > sub_v_size[l] - l)
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, sub_v_size[l] - l, cliques);
        else
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, l, cliques);

        return true;
    }

    if (L == 1) return false;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (sub_v_size[l] - G_deg[l][u] > L) {
            return false;
        }
    }

    F_size = 0;
    P_size = 0;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (G_deg[l][u] == sub_v_size[l] - 1) {
            loc[u] = -1;
            F[F_size++] = u;
            continue;
        }

        loc[u] = P_size;
        P[P_size++] = u;
    }

    int* e = new int [P_size * P_size]();

    for (i = 0; i < P_size; i++) {
        u = P[i];
        lack_size[u] = 0;

        end = DAG_deg[l][u];
        for (j = 0; j < end; j++) {
            v = DAG_adj[u][j];

            if (loc[v] != -1) {
                e[loc[u] * P_size + loc[v]] = 1;
                e[loc[v] * P_size + loc[u]] = 1;
            }
        }
    }

    for (i = 0; i < P_size * P_size; i++) {
        if (!e[i]) {
            j = i / P_size, k = i % P_size;
            u = P[j], v = P[k];
            lack[u][lack_size[u]++] = v;
        }
    }

    delete [] e;

    for(i = 0; i <= l; i++) {
        P_act = P_size;
        list_in_plex(0, i, l - i, cliques);
    }

    return true;
}

void EBBkC_Graph_t::truss_decompose(const char *dir) {
    graph_t g;

    //load the graph from file
    Graph G(dir);
    g.adj = G.edge_dst;
    g.num_edges = G.node_off;
    g.n = G.nodemax;
    g.m = G.edgemax;

    string reorder_method("core");

    vector <int32_t> new_vid_dict;
    vector <int32_t> old_vid_dict;
    ReorderWrapper(g, dir, reorder_method, new_vid_dict, old_vid_dict);

    //edge list array
    Timer get_eid_timer;

    Edge *edgeIdToEdge = (Edge *) malloc((g.m / 2) * sizeof(Edge));
    assert(edgeIdToEdge != nullptr);
    log_info("Malloc Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    //Populate the edge list array
    getEidAndEdgeList(&g, edgeIdToEdge);
    log_info("Init Eid Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    int *EdgeSupport = (int *) malloc(g.m / 2 * sizeof(int));
    assert(EdgeSupport != nullptr);

    auto max_omp_threads = omp_get_max_threads();
    omp_set_num_threads(max_omp_threads);
    log_info("Max Threads: %d", max_omp_threads);
#pragma omp parallel for
    for (auto i = 0; i < max_omp_threads; i++) {
        auto avg = g.m / 2 / max_omp_threads;
        auto iter_beg = avg * i;
        auto iter_end = (i == max_omp_threads - 1) ? g.m / 2 : avg * (i + 1);
        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
    }
    log_info("Init EdgeSupport Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    Timer global_timer;
    truss_num = PKT_intersection(&g, EdgeSupport, edgeIdToEdge);

#pragma omp single
    {
        int u, v, w;
        int *old2new = new int[N_NODES];
        for (int i = 0; i < N_NODES; i++) old2new[i] = -1;

        e_size = 0;
        edges = new Edge_t[g.m / 2];

        T = new int *[g.m / 2];
        T_size = new int[g.m / 2];

        for (int i = 0; i < g.m / 2; i++) {
            if (g.edge_truss[i] <= K) continue;

            Edge e = edgeIdToEdge[i];
            u = e.u;
            v = e.v;

            if (old2new[u] == -1) {
                old2new[u] = (int) new2old.size();
                new2old.push_back(u);
            }

            if (old2new[v] == -1) {
                old2new[v] = (int) new2old.size();
                new2old.push_back(v);
            }

            edges[e_size] = Edge_t(old2new[u], old2new[v], false);
            edge2id.insert(edges[e_size], e_size);
            rank.push_back(g.edge_rank[i]);

            int sz = g.v_set[i].size();
            T_size[e_size] = sz;

            T[e_size] = new int[truss_num + 1];
            for (int j = 0; j < sz; j++) {
                int w = g.v_set[i][j];

                if (old2new[w] == -1) {
                    old2new[w] = (int) new2old.size();
                    new2old.push_back(w);
                }

                T[e_size][j] = old2new[w];
            }

            e_size++;
        }

        C = new int *[e_size];
        C_size = new int[e_size];

        for (int i = 0; i < e_size; i++) {

            C_size[i] = 0;
            int sz = T_size[i];
            C[i] = new int[sz * (sz - 1) / 2];

            for (int j = 0; j < T_size[i]; j++) {
                for (int k = j + 1; k < T_size[i]; k++) {

                    Edge_t e = Edge_t(T[i][j], T[i][k], false);
                    int idx = edge2id.exist(e);

                    if (idx != -1 && rank[idx] > rank[i]) {
                        C[i][C_size[i]++] = idx;
                    }
                }
            }
        }
    };

    v_size = new2old.size();

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Truss number = %d\n", truss_num - 2);

    //Free memory
    free_graph(&g);
    free(edgeIdToEdge);
    free(EdgeSupport);

}


double EBBkC_t::list_k_clique(const char *file_name) {
    double runtime;
    struct rusage start, end;
    EBBkC_Graph_t G, g;

    printf("Reading edges from %s ...\n", file_name);
    GetCurTime(&start);
    G.truss_decompose(file_name);

    printf("Building necessary data structure ...\n");
    G.build(false);

    printf("Iterate over all cliques\n");

    g.truss_num = G.truss_num;
    g.build(true);
    for (int i = 0; i < G.e_size; i++) {
        G.branch(i, &g);
        g.EBBkC_plus_plus(K - 2, &N);
    }

    GetCurTime(&end);
    runtime = GetTime(&start, &end);

    return runtime;
}

double EBBkC_t::list_k_clique_parallel(const char *file_name) {
    double runtime = 0, runtime1 = 0;
    int n_edges, i;
//    struct rusage start, end;
    double start, end;
    EBBkC_Graph_t G, g;

    printf("Reading edges from %s ...\n", file_name);
    start = omp_get_wtime();
    G.truss_decompose(file_name);

    printf("Building necessary data structure ...\n");
    G.build(false);

    printf("Iterate over all cliques\n");

#pragma omp parallel private(g, start, end, n_edges, i) reduction(+:N)
    {
        n_edges = 0;
        g.truss_num = G.truss_num;
        g.build(true);
#pragma omp for schedule(dynamic, 1) nowait
        for (i = 0; i < G.e_size; i++) {
            G.branch(i, &g);
            g.EBBkC_plus_plus(K - 2, &N);
            n_edges++;
        }
        printf("Thread: %d, handled %d edges\n", omp_get_thread_num(), n_edges);
    }

    return (double) (omp_get_wtime() - start) * 1e3;
}
