
#include "vertex_oriented.h"
#include "set_operation.h"
#include <cassert>
#include <algorithm>
#include <set>

using namespace std;

extern const int K, L;
extern unsigned long long N;

VBBkC_Graph_t::VBBkC_Graph_t() = default;

VBBkC_Graph_t::~VBBkC_Graph_t() {
    int i;

    delete [] edges;

    for (i = 0; i < v_size; i++) delete [] T[i];
    delete [] T;

    for (i = 0; i < v_size; i++) delete [] C[i];
    delete [] C;

    delete [] T_size;

    delete [] C_size;

//    for (i = 0; i <= K; i++) delete [] sub_v[i];
//    delete [] sub_v;

    for (i = 0; i <= K; i++) delete [] sub_e[i];
    delete [] sub_e;

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    for (i = 0; i <= K; i++) delete [] DAG_deg[i];
    delete [] DAG_deg;

    for (i = 0; i <= K; i++) delete [] G_deg[i];
    delete [] G_deg;

    delete [] rank;

    delete [] col;

    for (i = 0; i < v_size; i++) delete [] DAG_adj[i];
    delete [] DAG_adj;

    for (i = 0; i < v_size; i++) delete [] G_adj[i];
    delete [] G_adj;

    for (i = 0; i <= K; i++) delete [] used[i];
    delete [] used;

    delete [] F;

    delete [] P;

    delete [] lack_size;

    for (i = 0; i < v_size; i++) delete [] lack[i];
    delete [] lack;

    delete [] lev;

    delete [] loc;
}

void VBBkC_Graph_t::read_edges_from_file(const char *file_name) {      // read from .clean file
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
        edge2id.insert(e, e_size);
        edges[e_size++] = e;
    }

    v_size = (int) new2old.size();

    core_decompose();

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Core number = %d\n", core_num);

    fclose(fp);
    delete [] old2new;
}

void VBBkC_Graph_t::core_decompose() {
    int i, j, k, s, t, idx;
    KeyVal_t kv;
    Heap_t heap;

    rank = new int [v_size];
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

    heap.make_heap(_d, v_size);
    for (i = 0; i < v_size; i++) {
        kv = heap.pop();
        core_num = (kv.val > core_num) ? kv.val : core_num;
        rank[kv.key] = v_size - i - 1;
        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++)
            heap.update(_adj[j]);
    }
    heap.release_heap();

    delete [] _d;
    delete [] _cd;
    delete [] _adj;

    T = new int* [v_size];
    for (i = 0; i < v_size; i++) T[i] = new int [core_num + 1];
    T_size = new int [v_size]();

    C = new int* [v_size];
    for (i = 0; i < v_size; i++) C[i] = new int [core_num * core_num];
    C_size = new int [v_size]();

    auto * _lab = new int [v_size];
    for (i = 0; i < v_size; i++) _lab[i] = K;

    for (i = 0; i < e_size; i++) {
        s = edges[i].s;
        t = edges[i].t;
        edges[i].s = (rank[s] > rank[t]) ? s : t;
        edges[i].t = (rank[s] > rank[t]) ? t : s;
        s = edges[i].s;
        t = edges[i].t;

        T[s][T_size[s]++] = t;
    }

    for (i = 0; i < v_size; i++) {

        for (j = 0; j < T_size[i]; j++) {
            s = T[i][j];
            _lab[s] = K - 1;
        }

        for (j = 0; j < T_size[i]; j++) {
            s = T[i][j];

            for (k = 0; k < T_size[s]; k++) {
                t = T[s][k];

                if (_lab[t] == K - 1) {
                    idx = edge2id.exist(Edge_t(s, t, false));
                    assert(idx != -1);
                    C[i][C_size[i]++] = idx;
                }
            }
        }

        for (j = 0; j < T_size[i]; j++) {
            s = T[i][j];
            _lab[s] = K;
        }
    }

    delete [] _lab;
}

void VBBkC_Graph_t::build_from_G() {
    int i;

    sub_v = new int* [K + 1];

    sub_e = new int* [K + 1];

    sub_e_size = new int [K + 1];

    sub_v_size = new int [K + 1];

    for (i = 0; i < K; i++) sub_v[i] = new int [core_num + 1];
    sub_v[K] = new int [v_size];

    for (i = 0; i < K; i++) sub_e[i] = new int [core_num * core_num];
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
    for (i = 0; i < v_size; i++) DAG_adj[i] = new int [core_num + 1];

    G_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) G_adj[i] = new int  [core_num + 1];

    used = new bool* [K + 1];
    for (i = 0; i <= K; i++) used[i] = new bool [v_size + 1]();

    F = new int [core_num + 1];

    P = new int [core_num + 1];

    lack_size = new int [v_size];

    lack = new int* [v_size];
    for (i = 0; i < v_size; i++) lack[i] = new int [L + 1];

    lev = new int [v_size]();

    loc = new int [v_size];
}

void VBBkC_Graph_t::DDegCol(int l, unsigned long long *cliques) {
    int c, i, j, k, e, e_, u, v, w, s, t, end;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        for (i = 0; i < v_size; i++) {

            if (T_size[i] < l - 1 || C_size[i] < (l - 1) * (l - 2) / 2) continue;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                col[u] = 0;
                DAG_deg[l - 1][u] = 0;
                G_deg[l - 1][u] = 0;
            }

            for (j = 0; j < C_size[i]; j++) {
                e = C[i][j];
                s = edges[e].s;
                t = edges[e].t;
                G_adj[s][G_deg[l - 1][s]++] = t;
                G_adj[t][G_deg[l - 1][t]++] = s;
            }

            auto *list = new KeyVal_t [core_num + 1];
            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                list[j].key = u;
                list[j].val = G_deg[l - 1][u];
            }
            sort(list, list + T_size[i]);

            for (j = 0; j < T_size[i]; j++) {
                u = list[j].key;
                for (k = 0; k < G_deg[l - 1][u]; k++) {
                    v = G_adj[u][k];
                    used[l][col[v]] = true;
                }
                for (c = 1; used[l][c]; c++) ;
                col[u] = c;
                for (k = 0; k < G_deg[l - 1][u]; k++) {
                    v = G_adj[u][k];
                    used[l][col[v]] = false;
                }
            }
            delete [] list;

            sub_v_size[l - 1] = 0;
            sub_e_size[l - 1] = 0;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                sub_v[l - 1][sub_v_size[l - 1]++] = u;
            }

            for (j = 0; j < C_size[i]; j++) {
                e_ = C[i][j];
                sub_e_size[l - 1]++;
                s = edges[e_].s;
                t = edges[e_].t;
                edges[e_].s = (col[s] > col[t]) ? s : t;
                edges[e_].t = (col[s] > col[t]) ? t : s;
                s = edges[e_].s;
                t = edges[e_].t;
                DAG_adj[s][DAG_deg[l - 1][s]++] = t;
            }

            DDegCol(l - 1, cliques);

        }

        return;
    }

    if (l == 1) {
        for (i = 0; i < sub_v_size[l]; i++) {
            (*cliques)++;
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

//    if (can_terminate(l, cliques)) {
//        return;
//    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = 0;
        sub_e_size[l - 1] = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;
        }

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

        DDegCol(l - 1, cliques);

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
        }
    }
}


void VBBkC_Graph_t::DDegree(int l, unsigned long long *cliques) {
    int i, j, k, e, e_, u, v, w, s, t, end;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        for (i = 0; i < v_size; i++) {

            if (T_size[i] < l - 1 || C_size[i] < (l - 1) * (l - 2) / 2) continue;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                DAG_deg[l - 1][u] = 0;
                G_deg[l - 1][u] = 0;
            }

            for (j = 0; j < C_size[i]; j++) {
                e = C[i][j];
                s = edges[e].s;
                t = edges[e].t;
                G_adj[s][G_deg[l - 1][s]++] = t;
                G_adj[t][G_deg[l - 1][t]++] = s;
            }

            auto *list = new KeyVal_t [core_num + 1];
            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                list[j].key = u;
                list[j].val = G_deg[l - 1][u];
            }
            sort(list, list + T_size[i]);

            for (j = 0; j < T_size[i]; j++) {
                u = list[j].key;
                rank[u] = j;            // DAG: large degree --> small degree
            }
            delete [] list;

            sub_v_size[l - 1] = 0;
            sub_e_size[l - 1] = 0;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                sub_v[l - 1][sub_v_size[l - 1]++] = u;
            }

            for (j = 0; j < C_size[i]; j++) {
                e_ = C[i][j];
                sub_e_size[l - 1]++;
                s = edges[e_].s;
                t = edges[e_].t;
                edges[e_].s = (rank[s] > rank[t]) ? s : t;
                edges[e_].t = (rank[s] > rank[t]) ? t : s;
                s = edges[e_].s;
                t = edges[e_].t;
                DAG_adj[s][DAG_deg[l - 1][s]++] = t;
            }

            DDegree(l - 1, cliques);
        }

        return;
    }

    if (l == 1) {
        for (i = 0; i < sub_v_size[l]; i++) {
            (*cliques)++;
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

//    if (can_terminate(l, cliques)) {
//        return;
//    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        sub_v_size[l - 1] = 0;
        sub_e_size[l - 1] = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;
        }

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

        DDegree(l - 1, cliques);

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
        }
    }
}

void VBBkC_Graph_t::SDegree(int l, unsigned long long *cliques) {
    int i, j, k, e, e_, u, v, w, s, t, end;

    if (sub_v_size[l] < l) return;

    if (l == K) {
        for (i = 0; i < v_size; i++) {

            if (T_size[i] < l - 1 || C_size[i] < (l - 1) * (l - 2) / 2) continue;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                DAG_deg[0][u] = 0;
                G_deg[l - 1][u] = 0;
            }

            for (j = 0; j < C_size[i]; j++) {
                e = C[i][j];
                s = edges[e].s;
                t = edges[e].t;
                G_adj[s][G_deg[l - 1][s]++] = t;
                G_adj[t][G_deg[l - 1][t]++] = s;
            }

            auto *list = new KeyVal_t [core_num + 1];
            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                list[j].key = u;
                list[j].val = G_deg[l - 1][u];
            }
            sort(list, list + T_size[i]);

            for (j = 0; j < T_size[i]; j++) {
                u = list[j].key;
                rank[u] = j;            // DAG: large degree --> small degree
            }
            delete [] list;

            sub_v_size[l - 1] = 0;
            sub_e_size[l - 1] = 0;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                sub_v[l - 1][sub_v_size[l - 1]++] = u;
            }

            sort(sub_v[l - 1], sub_v[l - 1] + sub_v_size[l - 1]);       // sorted array for SIMD usage

            for (j = 0; j < C_size[i]; j++) {
                e_ = C[i][j];
                sub_e_size[l - 1]++;
                s = edges[e_].s;
                t = edges[e_].t;
                edges[e_].s = (rank[s] > rank[t]) ? s : t;
                edges[e_].t = (rank[s] > rank[t]) ? t : s;
                s = edges[e_].s;
                t = edges[e_].t;
                DAG_adj[s][DAG_deg[0][s]++] = t;
            }

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                // sorted array for SIMD usage, DAG_adj can be considered as const.
                sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]);
            }

            SDegree(l - 1, cliques);
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

        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);

        if (l == 2) {
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                (*cliques)++;
            }
        }

        else {
            if (sub_v_size[l - 1] >= l - 1) {
                SDegree(l - 1, cliques);
            }
        }
    }
}

void VBBkC_Graph_t::BitCol(int l, unsigned long long *cliques) {
    int c, i, j, k, e, e_, u, v, w, s, t, end;

    if (sub_v_size[l] < l) return;

    if (l == K) {
        for (i = 0; i < v_size; i++) {

            if (T_size[i] < l - 1 || C_size[i] < (l - 1) * (l - 2) / 2) continue;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                col[u] = 0;
                DAG_deg[0][u] = 0;
                G_deg[l - 1][u] = 0;
            }

            for (j = 0; j < C_size[i]; j++) {
                e = C[i][j];
                s = edges[e].s;
                t = edges[e].t;
                G_adj[s][G_deg[l - 1][s]++] = t;
                G_adj[t][G_deg[l - 1][t]++] = s;
            }

            auto *list = new KeyVal_t [core_num + 1];
            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                list[j].key = u;
                list[j].val = G_deg[l - 1][u];
            }
            sort(list, list + T_size[i]);

            for (j = 0; j < T_size[i]; j++) {
                u = list[j].key;
                for (k = 0; k < G_deg[l - 1][u]; k++) {
                    v = G_adj[u][k];
                    used[l][col[v]] = true;
                }
                for (c = 1; used[l][c]; c++) ;
                col[u] = c;
                for (k = 0; k < G_deg[l - 1][u]; k++) {
                    v = G_adj[u][k];
                    used[l][col[v]] = false;
                }
            }
            delete [] list;

            sub_v_size[l - 1] = 0;
            sub_e_size[l - 1] = 0;

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                sub_v[l - 1][sub_v_size[l - 1]++] = u;
            }

            sort(sub_v[l - 1], sub_v[l - 1] + sub_v_size[l - 1]);       // sorted array for SIMD usage

            for (j = 0; j < C_size[i]; j++) {
                e_ = C[i][j];
                sub_e_size[l - 1]++;
                s = edges[e_].s;
                t = edges[e_].t;
                edges[e_].s = (col[s] > col[t]) ? s : t;
                edges[e_].t = (col[s] > col[t]) ? t : s;
                s = edges[e_].s;
                t = edges[e_].t;
                DAG_adj[s][DAG_deg[0][s]++] = t;
            }

            for (j = 0; j < T_size[i]; j++) {
                u = T[i][j];
                // sorted array for SIMD usage, DAG_adj can be considered as const.
                sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]);
            }

            BitCol(l - 1, cliques);

        }

        return;
    }

    if (l == 1) {
        for (i = 0; i < sub_v_size[l]; i++) {
            (*cliques)++;
        }

        return;
    }

//    if (can_terminate(l, cliques)) {
//        return;
//    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);

        if (l == 2) {
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                (*cliques)++;
            }
        }

        else {
            if (sub_v_size[l - 1] >= l - 1) {
                BitCol(l - 1, cliques);
            }
        }
    }
}

void VBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long *cliques) {
    if (picked == k) {
        (*cliques)++;
        return;
    }

    for (int i = start; i < list_size; i++) {
        VBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques);
    }
}

void VBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (p == 0) {
        if (q > F_size - q)
            VBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        else
            VBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
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

bool VBBkC_Graph_t::can_terminate(int l, unsigned long long *cliques) {
    int i, j, k, u, v, end;

    if (sub_e_size[l] < sub_v_size[l] * (sub_v_size[l] - L) / 2) return false;

    if (sub_e_size[l] == sub_v_size[l] * (sub_v_size[l] - 1) / 2) {
        VBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, min(l, sub_v_size[l] - l), cliques);
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

double VBBkC_t::list_k_clique(const char* file_name, int type) {
    double runtime;
    struct rusage start, end;
    VBBkC_Graph_t G;

    printf("Reading edges from %s ...\n", file_name);
    G.read_edges_from_file(file_name);

    printf("Building necessary data structure ...\n");
    G.build_from_G();

    printf("Iterate over all cliques\n");
    if (type == 0) {
        GetCurTime(&start);
        G.DDegCol(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 1) {
        GetCurTime(&start);
        G.DDegree(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 2) {
        GetCurTime(&start);
        G.SDegree(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 3) {
        GetCurTime(&start);
        G.BitCol(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else {
        runtime = 0;
//#pragma omp parallel default(none) private(g, start, end) shared(G, K) reduction(max:runtime) reduction(+:N)
//        {
//            GetCurTime(&start);
//            G.core_DAG();
//            G.build_from_DAG(K);
//#pragma omp for schedule(dynamic, 1) nowait
//            for (int i = 0; i < G.m; i++) {
//                G.vertex_oriented_twice_branch(&g, i);
//                g.color_DAG();
//                g.build_from_DAG(K - 2);
//                g.list_clique_plus(K - 2, &N);
//            }
//            GetCurTime(&end);
//            runtime = GetTime(&start, &end);
//        }
    }

    return runtime;
}
