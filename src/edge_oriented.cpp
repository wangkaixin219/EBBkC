#include "edge_oriented.h"
#include <set>
#include <algorithm>
#include <unordered_map>
#include "set_operation.h"

extern const int K, L;
extern unsigned long long N;

EBBkC_Graph_t::EBBkC_Graph_t() = default;

EBBkC_Graph_t::~EBBkC_Graph_t() {
    int i;

    delete [] edges;

    for (i = 0; i < e_size; i++) delete [] T[i];
    delete [] T;

    for (i = 0; i < e_size; i++) delete [] C[i];
    delete [] C;

    delete [] T_size;

    delete [] C_size;

    for (i = 0; i <= K; i++) delete [] sub_v[i];
    delete [] sub_v;

    for (i = 0; i <= K; i++) delete [] sub_e[i];
    delete [] sub_e;

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    for (i = 0; i <= K; i++) delete [] DAG_deg[i];
    delete [] DAG_deg;

    for (i = 0; i <= K; i++) delete [] G_deg[i];
    delete [] G_deg;

    delete [] col;

    for (i = 0; i < v_size; i++) delete [] DAG_adj[i];
    delete [] DAG_adj;

    for (i = 0; i < v_size; i++) delete [] G_adj[i];
    delete [] G_adj;

    for (i = 0; i <= K; i++) delete [] used[i];
    delete [] used;

    delete [] v_lab;

    delete [] e_lab;

    for (i = 0; i <= K; i++) delete [] out_v_size[i];
    delete [] out_v_size;

    for (i = 0; i <= K; i++) delete [] out_e_size[i];
    delete [] out_e_size;

    delete [] F;

    delete [] P;

    delete [] lack_size;

    for (i = 0; i < v_size; i++) delete [] lack[i];
    delete [] lack;

    delete [] lev;

    delete [] loc;
}

void EBBkC_Graph_t::read_edges_from_file(const char *file_name) {
    FILE *fp;
    if ((fp = fopen(file_name, "r")) == nullptr) {
        printf("Cannot open file %s.\n", file_name);
        exit(0);
    }

    Edge_t e;
    int u, v, i;
    int *old2new = new int [N_NODES];
    for (i = 0; i < N_NODES; i++) old2new[i] = -1;

    e_size = 0;
    edges = new Edge_t [N_EDGES];

    while (fscanf(fp, "%d %d%*[^\n]%*c", &u, &v) == 2) {

        if (u > N_NODES || v > N_NODES) {
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
            exit(0);
        }
        if (old2new[u] == -1) {
            old2new[u] = (int) new2old.size();
            new2old.push_back(u);
        }
        if (old2new[v] == -1) {
            old2new[v] = (int) new2old.size();
            new2old.push_back(v);
        }

        e = Edge_t(old2new[u], old2new[v], false);
        edges[e_size++] = e;
    }

    v_size = (int) new2old.size();

    fclose(fp);

    delete [] old2new;
}

void EBBkC_Graph_t::read_ordered_edges_from_file(const char *file_name) {
    FILE *fp;
    if ((fp = fopen(file_name, "r")) == nullptr) {
        printf("Cannot open file %s.\n", file_name);
        exit(0);
    }

    Edge_t e, e_;
    int u, v, w, i, j, k, idx, edge_rank, edge_sub_size;
    int *old2new = new int [N_NODES];
    for (i = 0; i < N_NODES; i++) old2new[i] = -1;
    vector<int> t_;
    vector<vector<int>> T_;

    e_size = 0;
    edges = new Edge_t [N_EDGES];

    while (fscanf(fp, "%d %d %d %d %d", &u, &v, &k, &edge_rank, &edge_sub_size) == 5) {

        if (k <= K) {
            fscanf(fp, "%*[^\n]%*c");
            continue;
        }
        truss_num = truss_num > k ? truss_num : k;

        if (u > N_NODES || v > N_NODES) {
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
            exit(0);
        }
        if (old2new[u] == -1) {
            old2new[u] = (int) new2old.size();
            new2old.push_back(u);
        }
        if (old2new[v] == -1) {
            old2new[v] = (int) new2old.size();
            new2old.push_back(v);
        }

        e = Edge_t(old2new[u], old2new[v], false);
        rank.push_back(edge_rank);
        edge2id.insert(e, e_size);
        edges[e_size++] = e;

        for (j = 0; j < edge_sub_size; j++) {
            fscanf(fp, "%d", &w);
            if (old2new[w] == -1) {
                old2new[w] = (int) new2old.size();
                new2old.push_back(w);
            }
            t_.push_back(old2new[w]);
        }
        T_.push_back(t_);
        t_.clear();
    }

    v_size = (int) new2old.size();

    T = new int* [e_size];
    for (i = 0; i < e_size; i++) T[i] = new int [truss_num + 1];
    T_size = new int [e_size];

    C = new int* [e_size];
    for (i = 0; i < e_size; i++) C[i] = new int [T_[i].size() * (T_[i].size() - 1) / 2];
    C_size = new int [e_size];

    for (i = 0; i < e_size; i++) {
        T_size[i] = C_size[i] = 0;

        for (j = 0; j < T_[i].size(); j++) T[i][T_size[i]++] = T_[i][j];

        for (j = 0; j < T_size[i]; j++) {
            for (k = j + 1; k < T_size[i]; k++) {

                e = Edge_t(T[i][j], T[i][k], false);

                if ((idx = edge2id.exist(e)) != -1  && rank[idx] > rank[i]) {
                    C[i][C_size[i]++] = idx;
                }
            }
        }
    }

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Truss number = %d\n", truss_num - 2);

    fclose(fp);
    delete [] old2new;
}

void EBBkC_Graph_t::truss_decompose(const char* w_file_name) {

    int i, j, k, s, t, w, end, edge_seq = 1;

    auto *_d = new int [v_size]();
    auto *_cd = new int [v_size + 1];
    auto *_adj = new int [2 * e_size];

    for (i = 0; i < e_size; i++) {
        _d[edges[i].s]++;
        _d[edges[i].t]++;
    }
    _cd[0] = 0;
    for (i = 1; i < v_size + 1; i++) {
        _cd[i] = _cd[i - 1] + _d[i - 1];
        _d[i - 1] = 0;
    }
    for (i = 0; i < e_size; i++) {
        _adj[_cd[edges[i].s] + _d[edges[i].s]++] = edges[i].t;
        _adj[_cd[edges[i].t] + _d[edges[i].t]++] = edges[i].s;
    }

    KeyVal_t kv;
    Heap_t heap;
    Edge_t  e;

    auto *sup = new int [e_size]();
    auto *h_table = new bool [v_size]();
    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_id;
    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_truss;
    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_rank;
    unordered_map<Edge_t, vector<int>, Edge_t::Hash_Edge_t> edge_sub;

    for (i = 0; i < e_size; i++) {

        edge_id[edges[i]] = i;

        s = _d[edges[i].s] < _d[edges[i].t] ? edges[i].s : edges[i].t;
        t = _d[edges[i].s] < _d[edges[i].t] ? edges[i].t : edges[i].s;

        for (j = _cd[s]; j < _cd[s + 1]; j++) {
            w = _adj[j];
            h_table[w] = true;
        }

        for (j = _cd[t]; j < _cd[t + 1]; j++) {
            w = _adj[j];
            if (h_table[w])
                sup[i]++;
        }

        for (j = _cd[s]; j < _cd[s + 1]; j++) {
            w = _adj[j];
            h_table[w] = false;
        }
    }

    heap.make_heap(sup, e_size);

    for (k = 3; !heap.empty(); k++) {

        while (!heap.empty() && heap.min_element().val < k - 2) {
            kv = heap.pop();
            e = edges[kv.key];
            edge_truss[e] = k;
            edge_rank[e] = edge_seq++;

            s = _d[e.s] < _d[e.t] ? e.s : e.t;
            t = _d[e.s] < _d[e.t] ? e.t : e.s;

            end = _cd[s] + _d[s];
            for (j = _cd[s]; j < end; j++) {
                w = _adj[j];
                h_table[w] = true;
            }

            end = _cd[t] + _d[t];
            for (j = _cd[t]; j < end; j++) {
                w = _adj[j];

                if (w == s) {
                    _adj[j--] = _adj[--end];
                    _d[t]--;
                }

                if (h_table[w]) {
                    edge_sub[e].push_back(w);
                    heap.update(edge_id[Edge_t(s, w, false)]);
                    heap.update(edge_id[Edge_t(w, t, false)]);
                }
            }

            end = _cd[s] + _d[s];
            for (j = _cd[s]; j < end; j++) {
                w = _adj[j];
                h_table[w] = false;

                if (w == t) {
                    _adj[j--] = _adj[--end];
                    _d[s]--;
                }
            }
        }
    }

    heap.release_heap();

    FILE *fp = fopen(w_file_name, "w");

    for (i = 0; i < e_size; i++) {
        e = edges[i];
        fprintf(fp, "%d %d %d %d %ld ", new2old[e.s], new2old[e.t], edge_truss[e], edge_rank[e], edge_sub[e].size());
        for (j = 0; j < edge_sub[e].size(); j++) fprintf(fp, "%d ", new2old[edge_sub[e][j]]);
        fprintf(fp, "\n");
    }

    fclose(fp);

    delete [] _d;
    delete [] _cd;
    delete [] _adj;
    delete [] h_table;
}


void EBBkC_Graph_t::build_from_G() {
    int i;

    sub_v = new int* [K + 1];

    sub_e = new int* [K + 1];

    sub_e_size = new int [K + 1];

    sub_v_size = new int [K + 1];

    for (i = 0; i < K; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[K] = new int [v_size];

    for (i = 0; i < K; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[K] = new int [e_size];

    sub_v_size[K] = 0;
    for (i = 0; i < v_size; i++) sub_v[K][sub_v_size[K]++] = i;

    sub_e_size[K] = 0;
    for (i = 0; i < e_size; i++) sub_e[K][sub_e_size[K]++] = i;

    lab = new int [v_size];
    for (i = 0; i < v_size; i++) lab[i] = K;

    DAG_deg = new int* [K + 1];
    for (i = 0; i <= K; i++) DAG_deg[i] = new int [v_size];

    G_deg = new int* [K + 1];
    for (i = 0; i <= K; i++) G_deg[i] = new int [v_size];

    col = new int [v_size];

    DAG_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [K + 1];
    for (i = 0; i <= K; i++) used[i] = new bool [v_size + 1]();

    v_lab = new int [v_size];
    for (i = 0; i < v_size; i++) v_lab[i] = K;

    e_lab = new int [e_size];
    for (i = 0; i < e_size; i++) e_lab[i] = K;

    out_v_size = new int* [K + 1];
    for (i = 0; i <= K; i++) out_v_size[i] = new int [e_size];

    out_e_size = new int* [K + 1];
    for (i = 0; i <= K; i++) out_e_size[i] = new int [e_size];

    F = new int [truss_num + 1];

    P = new int [truss_num + 1];

    lack_size = new int [v_size];

    lack = new int* [v_size];
    for (i = 0; i < v_size; i++) lack[i] = new int [L + 1];

    lev = new int [v_size]();

    loc = new int [v_size];
}


void EBBkC_Graph_t::EBBkC(int l, unsigned long long *cliques) {
    int i, j, k, u, e, e_, _e, end;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {

            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                sort(T[e], T[e] + T_size[e]);
                sort(C[e], C[e] + C_size[e]);
            }
            
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;


                sub_v_size[l - 2] = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
//                    v_lab[u] = l - 2;
                }


                sub_e_size[l - 2] = 0;

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                    e_lab[e_] = l - 2;
                }


//                sub_e_size[l - 2] = 0;

//                sort(C[e], C[e] + C_size[e]);

//                for (j = 0; j < C_size[e]; j++) {
//                    e_ = C[e][j];
//                    sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                    e_lab[e_] = l - 2;
//                }

//                sort(sub_e[l - 2], sub_e[l - 2] + sub_e_size[l - 2]);


//                sort(T[e], T[e] + T_size[e]);

//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
//                    v_lab[u] = l - 2;
//                }

//                sort(sub_v[l - 2], sub_v[l - 2] + sub_v_size[l - 2]);

//                for (j = 0; j < sub_e_size[l - 2]; j++) {
//                    e_ = sub_e[l - 2][j];
//                    out_e_size[l - 2][e_] = 0;
//
//                    end = C_size[e_];
//                    for (k = 0; k < end; k++) {
//                        _e = C[e_][k];
//
//                        if (e_lab[_e] == l - 2) {
//                            out_e_size[l - 2][e_]++;
//                        }
//                        else {
//                            C[e_][k--] = C[e_][--end];
//                            C[e_][end] = _e;
//                        }
//                    }
//
//                    out_v_size[l - 2][e_] = 0;
//
//                    end = T_size[e_];
//                    for (k = 0; k < end; k++) {
//                        u = T[e_][k];
//                        if (v_lab[u] == l - 2) {
//                            out_v_size[l - 2][e_]++;
//                        }
//                        else {
//                            T[e_][k--] = T[e_][--end];
//                            T[e_][end] = u;
//                        }
//                    }
//                }

                EBBkC(l - 2, cliques);

//                for (j = 0; j < C_size[e]; j++) {
//                    e_ = C[e][j];
//                    e_lab[e_] = l;
//                }
//
//                for (j = 0; j < T_size[e]; j++) {
//                    u = T[e][j];
//                    v_lab[u] = l;
//                }
            }
        }

        return;
    }

//    if (l == 2) {
//        for (i = 0; i < sub_e_size[l]; i++) {
//            (*cliques)++;
//        }
//
//        return;
//    }
//
//    if (l == 3) {
//        for (i = 0; i < sub_e_size[l]; i++) {
//            e = sub_e[l][i];
//
//            for (j = 0; j < out_v_size[l][e]; j++) {
//                u = T[e][j];
//                v_lab[u] = l - 2;
//            }
//
//            for (j = 0; j < sub_v_size[l]; j++) {
//                u = sub_v[l][j];
//                if (v_lab[u] == l - 2) (*cliques)++;
//            }
//
//            for (j = 0; j < out_v_size[l][e]; j++) {
//                u = T[e][j];
//                v_lab[u] = l;
//            }
//        }
//        return;
//    }

    for (i = 0; i < sub_e_size[l]; i++) {
        e = sub_e[l][i];

        if (l == 3) {
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            for (j = 0; j < sub_v_size[l - 2]; j++) {
                (*cliques)++;
            }
        }

        else if (l == 4) {
             sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
             for (j = 0; j < sub_e_size[l - 2]; j++) {
                (*cliques)++;
             }
        }

        else {
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
            EBBkC(l - 2, cliques);
        }

//        if (out_v_size[l][e] < l - 2 || out_e_size[l][e] < (l - 2) * (l - 3) / 2) continue;
//
//        sub_v_size[l - 2] = 0;
//        sub_e_size[l - 2] = 0;
//
//        for (j = 0; j < out_e_size[l][e]; j++) {
//            e_ = C[e][j];
//            sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//            e_lab[e_] = l - 2;
//        }
//
//        for (j = 0; j < out_v_size[l][e]; j++) {
//            u = T[e][j];
//            sub_v[l - 2][sub_v_size[l - 2]++] = u;
//            v_lab[u] = l - 2;
//        }
//
//        for (j = 0; j < sub_e_size[l - 2]; j++) {
//            e_ = sub_e[l - 2][j];
//            out_e_size[l - 2][e_] = 0;
//
//            end = out_e_size[l][e_];
//            for (k = 0; k < end; k++) {
//                _e = C[e_][k];
//
//                if (e_lab[_e] == l - 2) {
//                    out_e_size[l - 2][e_]++;
//                }
//                else {
//                    C[e_][k--] = C[e_][--end];
//                    C[e_][end] = _e;
//                }
//            }
//
//            out_v_size[l - 2][e_] = 0;
//
//            end = out_v_size[l][e_];
//            for (k = 0; k < end; k++) {
//                u = T[e_][k];
//
//                if (v_lab[u] == l - 2) {
//                    out_v_size[l - 2][e_]++;
//                }
//                else {
//                    T[e_][k--] = T[e_][--end];
//                    T[e_][end] = u;
//                }
//            }
//        }
//
//        EBBkC(l - 2, cliques);
//
//        for (j = 0; j < out_e_size[l][e]; j++) {
//            e_ = C[e][j];
//            e_lab[e_] = l;
//        }
//
//        for (j = 0; j < out_v_size[l][e]; j++) {
//            u = T[e][j];
//            v_lab[u] = l;
//        }
    }
}


void EBBkC_Graph_t::EBBkC_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                 e = sub_e[l][i];

                 for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                 }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    col[u] = 0;
                    DAG_deg[0][u] = 0;
                    G_deg[l - 2][u] = 0;
                }

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    s = edges[e_].s;
                    t = edges[e_].t;
                    G_adj[s][G_deg[l - 2][s]++] = t;
                    G_adj[t][G_deg[l - 2][t]++] = s;
                }

                auto *list = new KeyVal_t [truss_num + 1];
                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    list[j].key = u;
                    list[j].val = G_deg[l - 2][u];
                }
                sort(list, list + T_size[e]);

                for (j = 0; j < T_size[e]; j++) {
                    u = list[j].key;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = true;
                    }
                    for (c = 1; used[K][c]; c++) ;
                    col[u] = c;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = false;
                    }
                }
                delete [] list;

                sub_v_size[l - 2] = 0;
                dist = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
                    if (!used[K][col[u]]) {
                        used[K][col[u]] = true;
                        dist++;
                    }
                }

                if (dist >= l - 2) {
                    sort(sub_v[l - 2], sub_v[l - 2] + sub_v_size[l - 2]);

                    sub_e_size[l - 2] = 0;
                    for (j = 0; j < C_size[e]; j++) {
                        e_ = C[e][j];
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        edges[e_].s = (col[s] > col[t]) ? s : t;
                        edges[e_].t = (col[s] > col[t]) ? t : s;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        DAG_adj[s][DAG_deg[0][s]++] = t;
                    }

                    for (j = 0; j < T_size[i]; j++) {
                        u = T[e][j];
                        // sorted array for SIMD usage, DAG_adj can be considered as const.
                        sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]);
                    }

                    EBBkC_plus(l - 2, cliques);
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    used[K][col[u]] = false;
                }
            }
        }

        return;
    }

    if (l == 1) {
        for (i = 0; i < sub_v_size[l]; i++) {
            (*cliques)++;
        }

        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);

        if (l == 2) {
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                (*cliques)++;
            }
        }

//        else {
//            if (sub_v_size[l - 1] >= l - 1) {
//
//                for (j = 0; j < sub_v_size[l - 1]; j++) {
//                    v = sub_v[l - 1][j];
//
//                    if (col[v] < l - 1) continue;
//
//                    sub_v_size[l - 2] = intersect_simd4x(sub_v[l - 1], sub_v_size[l - 1], DAG_adj[v], DAG_deg[0][v], sub_v[l - 2]);
//
//                    if (sub_v_size[l - 2] >= l - 2) {
//                        EBBkC_plus(l - 2, cliques);
//                    }
//                }
//            }
//        }

        else {
            if (sub_v_size[l - 1] >= l - 1) {
                EBBkC_plus(l - 1, cliques);
            }
        }
    }



//    if (l == 2) {
//        for (i = 0; i < sub_v_size[l]; i++) {
//            u = sub_v[l][i];
//
//            for (j = 0; j < DAG_deg[l][u]; j++) {
//                (*cliques)++;
//            }
//        }
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
//                        // for efficiency issues, we don't save in sub_e[l]
//                    } else {
//                        DAG_adj[v][k--] = DAG_adj[v][--end];
//                        DAG_adj[v][end] = w;
//                    }
//                }
//            }
//
//            EBBkC_plus(l - 1, cliques);
//        }
//
//        for (j = 0; j < sub_v_size[l - 1]; j++) {
//            v = sub_v[l - 1][j];
//            lab[v] = l;
//            used[l][col[v]] = false;
//        }
//    }
}


void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];
                for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];
                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    col[u] = 0;
                    DAG_deg[l - 2][u] = 0;
                    G_deg[l - 2][u] = 0;
                }

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    s = edges[e_].s;
                    t = edges[e_].t;
                    G_adj[s][G_deg[l - 2][s]++] = t;
                    G_adj[t][G_deg[l - 2][t]++] = s;
                }

                auto *list = new KeyVal_t [truss_num + 1];
                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    list[j].key = u;
                    list[j].val = G_deg[l - 2][u];
                }
                sort(list, list + T_size[e]);

                for (j = 0; j < T_size[e]; j++) {
                    u = list[j].key;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = true;
                    }
                    for (c = 1; used[K][c]; c++) ;
                    col[u] = c;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = false;
                    }
                }
                delete [] list;

                sub_v_size[l - 2] = 0;
                dist = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
                    if (!used[K][col[u]]) {
                        used[K][col[u]] = true;
                        dist++;
                    }
                }

                if (dist >= l - 2) {
                    sub_e_size[l - 2] = 0;
                    for (j = 0; j < C_size[e]; j++) {
                        e_ = C[e][j];
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        edges[e_].s = (col[s] > col[t]) ? s : t;
                        edges[e_].t = (col[s] > col[t]) ? t : s;
                        s = edges[e_].s;
                        t = edges[e_].t;

                        DAG_adj[s][DAG_deg[l - 2][s]++] = t;
                    }

                    EBBkC_plus_plus(l - 2, cliques);
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    used[K][col[u]] = false;
                }
            }
        }

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
    if (p == 0) {
        if (q > F_size - q)
            EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        else
            EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
        return;
    }

    int i, j, u, v, vis = 0;

    for (i = start; i < P_size && (P_act >= p && F_size >= q); i++) {
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
        EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, min(l, sub_v_size[l] - l), cliques);
        return true;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (sub_v_size[l] - G_deg[l][u] > L) return false;
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


double EBBkC_t::truss_order(const char *r_file_name, const char *w_file_name) {
    double runtime;
    struct rusage start, end;
    EBBkC_Graph_t G;

    GetCurTime(&start);
    G.read_edges_from_file(r_file_name);
    G.truss_decompose(w_file_name);
    GetCurTime(&end);
    runtime = GetTime(&start, &end);

    return runtime;
}

double EBBkC_t::list_k_clique(const char *file_name, int type) {
    double runtime;
    struct rusage start, end;
    EBBkC_Graph_t G;

    printf("Reading edges from %s ...\n", file_name);
    G.read_ordered_edges_from_file(file_name);

    printf("Building necessary data structure ...\n");
    G.build_from_G();

    printf("Iterate over all cliques\n");
    if (type == 0) {
        GetCurTime(&start);
        G.EBBkC(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }
    else if (type == 1) {
        GetCurTime(&start);
        G.EBBkC_plus(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }
    else if (type == 2) {
        GetCurTime(&start);
        G.EBBkC_plus_plus(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }
    else {
        runtime = 0;
    }

    return runtime;
}