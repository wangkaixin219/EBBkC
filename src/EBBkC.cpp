#include "EBBkC.h"
#include <set>
#include <cassert>
#include <algorithm>
#include <map>

extern int K, P;
extern unsigned long long N;

EBBkC_Graph_t::EBBkC_Graph_t() = default;

EBBkC_Graph_t::~EBBkC_Graph_t() {
//    delete [] v_sizes;
//    delete [] cd;
//    delete [] adj;
//    delete [] rank;
//    delete [] new2old;
//    delete [] lab;
//    delete [] color;
//
//    if (d) {
//        for (int i = 0; i < depth + 1; ++i) delete [] d[i];
//        delete [] d;
//    }
//
////    if (od) {
////        for (int i = 0; i < depth + 1; ++i) delete [] od[i];
////        delete [] od;
////    }
//
//    if (sub) {
//        for (int i = 0; i < depth + 1; ++i) delete [] sub[i];
//        delete [] sub;
//
//    }
//
//    if (branch_edges) {
//        for (int i = 0; i < depth + 1; ++i) delete [] branch_edges[i];
//        delete [] branch_edges;
//    }
//    else {
//        delete [] edges;
//    }
//
//    delete [] plex_p;
//    delete [] plex_f;
//    for (int i = 0; i < v_size; i++) if(plex_d) delete [] plex_d[i];
//    delete [] plex_d;
//    delete [] plex_d_num;
//    delete [] plex_lev;
//    delete [] plex_new2old;
//    delete [] plex_loc;
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
    for (i = 0; i < e_size; i++) C[i] = new int [truss_num * truss_num];
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

//void EBBkC_Graph_t::color_DAG() {
//    int i, j, s, t, u, c, max_deg = 0;
//    Key_Val_t kv;
//
//    color_num = 1;
//    color = new int [v_size]();
//    auto *_d = new int [v_size]();
//    auto *_cd = new int [v_size + 1];
//    auto *_adj = new int [2 * edges.size()];
//    auto *used = new bool [v_size + 1]();
//    auto *list = new Key_Val_t [v_size];
//
//    for (i = 0; i < edges.size(); i++) {
//        assert(edges[i].s < v_size);
//        assert(edges[i].t < v_size);
//        _d[edges[i].s]++;
//        _d[edges[i].t]++;
//    }
//    _cd[0] = 0;
//    for (i = 1; i < v_size + 1; i++) {
//        _cd[i] = _cd[i - 1] + _d[i - 1];
//        max_deg = (max_deg > _d[i - 1]) ? max_deg : _d[i - 1];
//        list[i - 1].key = i - 1;
//        list[i - 1].val = _d[i - 1];
//        _d[i - 1] = 0;
//    }
//    for (i = 0; i < edges.size(); i++) {
//        _adj[_cd[edges[i].s] + _d[edges[i].s]++] = edges[i].t;
//        _adj[_cd[edges[i].t] + _d[edges[i].t]++] = edges[i].s;
//    }
//
//    sort(list, list + v_size);
//
//    for (i = 0; i < v_size; i++) {
//        kv = list[i];
//
//        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++) {
//            u = _adj[j];
//            used[color[u]] = true;
//        }
//
//        for (c = 1; used[c]; c++) ;
//        color[kv.key] = c;
//        color_num = (c > color_num) ? c : color_num;
//
//        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++) {
//            u = _adj[j];
//            used[color[u]] = false;
//        }
//    }
//
//    delete [] _d;
//    delete [] _cd;
//    delete [] _adj;
//    delete [] used;
//    delete [] list;
//
//    for (i = 0; i < edges.size(); i++) {
//        s = edges[i].s;
//        t = edges[i].t;
//
//        edges[i].s = (color[s] > color[t]) ? s : t;
//        edges[i].t = (color[s] > color[t]) ? t : s;
//    }
//}

//void EBBkC_Graph_t::build_from_DAG(int height) {
//    int s, t, i, max_deg = 0;
//    auto *_od = new int [v_size]();
//    auto *_d = new int [v_size]();
//    auto *_sub = new int [v_size];
//
//    cd = new int [v_size + 1];
////    adj = new int [edges.size()];
//    lab = new int [v_size];
//
//    for (i = 0; i < edges.size(); i++) {
//        _od[edges[i].s]++;
//        _d[edges[i].s]++;
//        _d[edges[i].t]++;
//    }
//    cd[0] = 0;
//    for (i = 1; i < v_size + 1; i++) {
//        cd[i] = cd[i - 1] + _od[i - 1];
//        max_deg = (_od[i - 1] > max_deg) ? _od[i - 1] : max_deg;
//        _od[i - 1] = 0;
//        _sub[i - 1] = i - 1;
//        lab[i - 1] = height;
//    }
//    for (i = 0; i < edges.size(); i++) {
//        s = edges[i].s, t = edges[i].t;
////        adj[cd[s] + _od[s]++] = t;
//    }
//
//    sub_v = vector<unordered_set<int>>(K);
//    sub_e_size = new int []; (K);
//
//    d = new int* [height + 1];
//    od = new int* [height + 1];
//
//    for (i = 0; i < height; i++) {
//        d[i] = new int [v_size];
//        od[i] = new int [v_size];
//    }
//    d[height] = _d;
//    od[height] = _od;
//
//
//    // for early-termination
//    plex_new2old = new int [v_size];
//    plex_f = new int [v_size];
//    plex_p = new int [v_size];
//    plex_d_num = new int [v_size];
//    plex_d = new int* [v_size];
//    for (i = 0; i < v_size; i++) plex_d[i] = new int [P + 1];
//    plex_lev = new int [v_size];
//    plex_loc = new int [v_size];
//}

void EBBkC_Graph_t::build_from_G() {
    int i, j, k, u, v, w;

    sub_v = new int* [K + 1];
    sub_e = new int* [K + 1];
    sub_e_size = new int [K + 1];
    sub_v_size = new int [K + 1];

    for (i = 0; i < K; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[K] = new int [v_size];

    for (i = 0; i < K; i++) sub_e[i] = new int [truss_num * truss_num];
    sub_e[K] = new int [e_size];

    sub_v_size[K] = 0;
    for (i = 0; i < v_size; i++) sub_v[K][sub_v_size[K]++] = i;

    sub_e_size[K] = 0;
    for (i = 0; i < e_size; i++) sub_e[K][sub_e_size[K]++] = i;

    lab = new int [v_size];
    for (i = 0; i < v_size; i++) lab[i] = K;

    out_deg = new int* [K + 1];
    for (i = 0; i <= K; i++) out_deg[i] = new int [v_size];

    deg = new int* [K + 1];
    for (i = 0; i <= K; i++) deg[i] = new int [v_size];

    col = new int [v_size];

    DAG_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [K + 1];
    for (i = 0; i <= K; i++) used[i] = new bool [truss_num + 1]();

//    auto *_d = new int [v_size]();
//    auto *_sub = new int [v_size];
//
//    cd = new int [v_size + 1];
//    adj = new int [2 * edges.size()];
//    lab = new int [v_size];
//
//    for (i = 0; i < edges.size(); i++) {
//        _d[edges[i].s]++;
//        _d[edges[i].t]++;
//    }
//
//    cd[0] = 0;
//    for (i = 1; i < v_size + 1; i++) {
//        cd[i] = cd[i - 1] + _d[i - 1];
//        _d[i - 1] = 0;
//        _sub[i - 1] = i - 1;
//        lab[i - 1] = K;
//    }
//    for (i = 0; i < edges.size(); i++) {
//        u = edges[i].s, v = edges[i].t;
//        adj[cd[u] + _d[u]++] = v;
//        adj[cd[v] + _d[v]++] = u;
//    }
//
//    d = new int* [K + 1];
//
//    for (i = 0; i < K; i++) {
//        d[i] = new int [v_size];
//    }
//
//    d[K] = _d;
//
//    // for early-termination
//    plex_new2old = new int [v_size];
//    plex_f = new int [v_size];
//    plex_p = new int [v_size];
//    plex_d_num = new int [v_size];
//    plex_d = new int* [v_size];
//    for (i = 0; i < v_size; i++) plex_d[i] = new int [P + 1];
//    plex_lev = new int [v_size];
//    plex_loc = new int [v_size];
}

//void EBBkC_Graph_t::edge_oriented_branch(EBBkC_Graph_t *sg, int edge_id) {
//    Edge_t e;
//    int i, j, k, u, v, w, end;
//    static int *old2new = nullptr;
//
//#pragma omp threadprivate(old2new)
//    if (!old2new) {
//        old2new = new int [v_size];
//    }
//
//    sg->n = 0;
//    sg->new2old = new int [truss_num];
//    e = edges[edge_id];
//
//    for (i = 0; i < edge2sub[e].size(); i++) {
//        w = edge2sub[e][i];
//        old2new[w] = sg->n;
//        sg->new2old[sg->n++] = w;
//    }
//
//    assert(sg->n <= truss_num);
//
//    sg->m = 0;
//    sg->edges = new Edge_t [sg->n * sg->n];
//
//    for (i = 0; i < sg->n; i++) {
//        for (j = i + 1; j < sg->n; j++) {
//            u = sg->new2old[i], v = sg->new2old[j];
//            e = Edge_t(u, v);
//            if (edge2id.count(e) && rank[edge2id[e]] > rank[edge_id]) {
//                sg->edges[sg->m++] = Edge_t(i, j);
//            }
//        }
//    }
//}


void EBBkC_Graph_t::basic_clique_list(int l, unsigned long long *cliques) {
//    int i, j, e;
//
//    if (sub_v[l].size() < l || sub_e_size[l] < l * (l - 1) / 2) return;
//
//    if (l == K) {
//        if (K == 3) {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//                (*cliques) += T[e].size();
//            }
//        }
//        else {
//            for (i = 0; i < sub_e_size[l]; i++) {
//                e = sub_e[l][i];
//
//                if (T[e].size() < l - 2 || C[e].size() < (l - 2) * (l - 3) / 2) continue;
//
//                sub_v[l - 2].clear();
//                sub_e_size[l - 2] = 0;
//
//                for (auto v : T[e]) sub_v[l - 2][sub_v_size[l - 2]++] = v;
//                for (auto e_ : C[e]) sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//
//                basic_clique_list(l - 2, cliques);
//            }
//        }
//        printf("r3 = %.3lf, r4 = %.3lf ms\n", r3, r4);
//        return;
//    }
//
//    if (l == 2) {
//        (*cliques) += sub_e_size[l];
//        return;
//    }
//
//    if (l == 3) {
//        for (i = 0; i < sub_e_size[l]; i++) {
//            e = sub_e[l][i];
//            for (auto v : sub_v[l])
//                if (T[e].count(v))
//                    (*cliques)++;
//        }
//        return;
//    }
//
//    for (i = 0; i < sub_e_size[l]; i++) {
//        e = sub_e[l][i];
//        GetCurTime(&start1);
//        sub_v[l - 2].clear();
//        GetCurTime(&end1);
//        r3 += GetTime(&start1, &end1);
//        sub_e_size[l - 2] = 0;
//
//        for (j = 0; j < sub_e_size[l]; j++) {
//            int e_ = sub_e[l][j];
//
//            GetCurTime(&start1);
//            bool sss = C[e].find(e_) != C[e].end();
//            GetCurTime(&end1);
//            r4 += GetTime(&start1, &end1);
//
//            if (sss) {
//                sub_e[l - 2][sub_e_size[l - 2]++] = e_;
//                sub_v[l - 2].insert(edges[e_].s);
//                sub_v[l - 2].insert(edges[e_].t);
//            }
//        }
//
//        basic_clique_list(l - 2, cliques);
//    }
}

void EBBkC_Graph_t::basic_clique_list_plus(int l, unsigned long long *cliques) {

}

int c1 = 0, c2 = 0, c3 = 0;
rusage start1, end1;
double r;

void EBBkC_Graph_t::improved_clique_list(int l, unsigned long long *cliques) {

    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                 e = sub_e[l][i];
                (*cliques) += T_size[e];
            }
        }
        else {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) {
                    c1++;
                    continue;
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    col[u] = 0;
                    out_deg[l - 2][u] = 0;
                    deg[l - 2][u] = 0;
                }

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    s = edges[e_].s;
                    t = edges[e_].t;
                    G_adj[s][deg[l - 2][s]++] = t;
                    G_adj[t][deg[l - 2][t]++] = s;
                }

                auto *list = new KeyVal_t [truss_num + 1];
                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    list[j].key = u;
                    list[j].val = deg[l - 2][u];
                }
                sort(list, list + T_size[e]);

                for (j = 0; j < T_size[e]; j++) {
                    u = list[j].key;
                    for (k = 0; k < deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = true;
                    }
                    for (c = 1; used[K][c]; c++) ;
                    col[u] = c;
                    for (k = 0; k < deg[l - 2][u]; k++) {
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
                        DAG_adj[s][out_deg[l - 2][s]++] = t;
                    }

//                    GetCurTime(&start1);
                    improved_clique_list(l - 2, cliques);
//                    GetCurTime(&end1);
//                    r += GetTime(&start1, &end1);
                } else {
                    c2++;
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    used[K][col[u]] = false;
                }
            }
        }
        printf("c1 = %d, c2 = %d, runtime = %.3lf ms\n", c1, c2, r);
        return;
    }

    if (l == 2) {
        (*cliques) += sub_e_size[l];
        return;
    }

//    auto *list = new KeyVal_t [truss_num + 1];
//    for (i = 0; i < sub_v_size[l]; i++) {
//        u = sub_v[l][i];
//        list[i].key = u;
//        list[i].val = col[u];
//    }
//    sort(list, list + sub_v_size[l]);
//    for (i = 0; i < sub_v_size[l]; i++) {
//        sub_v[l][i] = list[i].key;
//    }
//    delete [] list;

    if (l == 3) {

        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];

            if (col[u] < l) continue;

            for (j = 0; j < out_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l - 1;
            }

            for (j = 0; j < out_deg[l][u]; j++) {
                v = DAG_adj[u][j];

                if (col[v] < l - 1) continue;

                for (k = 0; k < out_deg[l][v]; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) (*cliques)++;
                }
            }

            for (j = 0; j < out_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l;
            }
        }

//        for (i = 0; i < sub_e_size[l]; i++) {
//            e = sub_e[l][i];
//            s = edges[e].s, t = edges[e].t;
//
//            for (j = 0; j < out_deg[l][t]; j++) {
//                w = DAG_adj[t][j];
//                lab[w] = l - 2;
//            }
//
//            for (j = 0; j < out_deg[l][s]; j++) {
//                w = DAG_adj[s][j];
//                if (lab[w] == l - 2) (*cliques)++;
//            }
//
//            for (j = 0; j < out_deg[l][t]; j++) {
//                w = DAG_adj[t][j];
//                lab[w] = l;
//            }
//        }
        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = 0;
        dist = 0;

        for (j = 0; j < out_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            out_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }

        if (dist >= l - 1) {

            sub_e_size[l - 1] = 0;
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                v = sub_v[l - 1][j];

                end = out_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        out_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;
//                        sub_e[l - 1][sub_e_size[l - 1]++] = edge2id.exist(Edge_t(v, w, true));
                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }

            improved_clique_list(l - 1, cliques);
        }

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }

//        for (j = 0; j < out_deg[l][u]; j++) {
//            v = DAG_adj[u][j];
//            lab[v] = l - 1;
//        }
//
//        for (j = 0; j < out_deg[l][u]; j++) {
//            v = DAG_adj[u][j];
//
//            if (col[v] < l - 1) continue;
//
//            sub_v_size[l - 2] = 0;
//            dist = 0;
//
//            for (k = 0; k < out_deg[l][v]; k++) {
//                w = DAG_adj[v][k];
//                if (lab[w] == l - 1) {
//                    lab[w] = l - 2;
//                    sub_v[l - 2][sub_v_size[l - 2]++] = w;
//                    out_deg[l - 2][w] = 0;
//
//                    if (!used[l][col[w]]) {
//                        used[l][col[w]] = true;
//                        dist++;
//                    }
//                }
//            }
//
//            if (dist >= l - 2) {
//
//                sub_e_size[l - 2] = 0;
//                for (k = 0; k < sub_v_size[l - 2]; k++) {
//                    s = sub_v[l - 2][k];
//
//                    end = out_deg[l][s];
//                    for (p = 0; p < end; p++) {
//                        t = DAG_adj[s][p];
//                        if (lab[t] == l - 2) {
//                            out_deg[l - 2][s]++;
//                            sub_e_size[l - 2]++;
////                            sub_e[l - 2][sub_e_size[l - 2]++] = edge2id[Edge_t(s, t)];
//                        } else {
//                            DAG_adj[s][p--] = DAG_adj[s][--end];
//                            DAG_adj[s][end] = t;
//                        }
//                    }
//                }
//
//                improved_clique_list(l - 2, cliques);
//            }
//
//            for (k = 0; k < sub_v_size[l - 2]; k++) {
//                w = sub_v[l - 2][k];
//                lab[w] = l - 1;
//                used[l][col[w]] = false;
//            }
//        }
//
//        for (k = 0; k < out_deg[l][u]; k++) {
//            w = DAG_adj[u][k];
//            lab[w] = l;
//        }
    }
}

//
//void EBBkC_Graph_t::improved_clique_list_plus(int l, unsigned long long *cliques) {
//    int u, v, w, i, j, k, end;
//
//    if (v_sizes[l] < l)    return;
//
//    if (l == 2) {
//        for (i = 0; i < v_sizes[l]; i++) {
//            u = sub[l][i];
//            end = cd[u] + od[l][u];
//            for (j = cd[u]; j < end; j++) {
//                (*cliques)++;
//            }
//        }
//        return;
//    }
//
//    if (l == 1) {
//        (*cliques) += v_sizes[l];
//        return;
//    }
//
//    if (list_clique_in_plex(l, cliques)) return;
//
//    for (i = 0; i < v_sizes[l]; i++) {
//        u = sub[l][i];
//
//        if (color[u] < l) continue;
//
//        v_sizes[l - 1] = 0;
//        end = cd[u] + od[l][u];
//        for (j = cd[u]; j < end; j++) {
//            v = adj[j];
//            lab[v] = l - 1;
//            sub[l - 1][v_sizes[l - 1]++] = v;
//            od[l - 1][v] = 0;
//            d[l - 1][v] = 0;
//
//        }
//
//        e_size = 0;
//
//        for (j = 0; j < v_sizes[l - 1]; j++) {
//            v = sub[l - 1][j];
//            end = cd[v] + od[l][v];
//            for (k = cd[v]; k < end; k++) {
//                w = adj[k];
//                if (lab[w] == l - 1) {
//                    od[l - 1][v]++;
//                    d[l - 1][v]++;
//                    d[l - 1][w]++;
//                    e_size++;
//                }
//                else {
//                    adj[k--] = adj[--end];
//                    adj[end] = w;
//                }
//            }
//        }
//
//        improved_clique_list_plus(l - 1, cliques);
//
//        for (j = 0; j < v_sizes[l - 1]; j++) {
//            v = sub[l - 1][j];
//            lab[v] = l;
//        }
//    }
//}
//
//
//void EBBkC_Comb_list(int *list, int list_size, int start, int l, int k, unsigned long long *cliques) {
//    if (l == k) {
//        (*cliques)++;
//        return;
//    }
//
//    for (int i = start; i < list_size; i++) {
//        EBBkC_Comb_list(list, list_size, i + 1, l + 1, k, cliques);
//    }
//}
//
//void EBBkC_Graph_t::list_in_plex(int p, int q, unsigned long long *cliques) {
//    if (p == 0) {
//        EBBkC_Comb_list(plex_f, plex_f_size, 0, 0, min(q, plex_f_size - q), cliques);
//        return;
//    }
//
//    int i, j, u, v, prev, vis = 0;
//
//    for (i = plex_sub; i < plex_size && (plex_p_size >= p && plex_f_size>= q); i++) {
//        u = plex_p[i];
//        if (plex_lev[u]) continue;
//
//        for (j = 0; j < plex_d_num[u]; j++) {
//            v = plex_d[u][j];
//            if (plex_loc[v] >= i && plex_lev[v] == 0) {
//                plex_lev[v] = p;
//                plex_p_size--;
//                assert(plex_p_size >= 0);
//            }
//        }
//
//        plex_sub = i + 1;
//        list_in_plex(p - 1, q, cliques);
//
//        for (j = 0; j < plex_d_num[u]; j++) {
//            v = plex_d[u][j];
//            if (plex_loc[v] >= i && plex_lev[v] == p) {
//                plex_lev[v] = 0;
//                plex_p_size++;
//            }
//        }
//
//        plex_p_size--;
//        vis++;
//    }
//    plex_p_size += vis;
//}
//
//bool EBBkC_Graph_t::list_clique_in_plex(int l, unsigned long long* cliques) {
//    int i, j, u, v, end, dis_m = 0, dis_num;
//
//    if (e_size < v_sizes[l] * (v_sizes[l] - P) / 2) return false;
//
//    if (e_size == v_sizes[l] * (v_sizes[l] - 1) / 2) {
//        EBBkC_Comb_list(sub[l], v_sizes[l], 0, 0, min(l, v_sizes[l] - l), cliques);
//        return true;
//    }
//
//    for (i = 0; i < v_sizes[l]; i++) {
//        u = sub[l][i];
//        dis_num = v_sizes[l] - d[l][u];
//        dis_m = dis_m > dis_num ? dis_m : dis_num;
//        if (dis_m > P) return false;
//    }
//
//    if (dis_m >= 2) {
//        int *old2new = new int [v_size];
//        int *e = new int [v_sizes[l] * v_sizes[l]]();
//
//        for (i = 0; i < v_sizes[l]; i++) {
//            u = sub[l][i];
//            plex_new2old[i] = u;
//            plex_d_num[i] = 0;
//            plex_lev[i] = 0;
//            old2new[u] = i;
//        }
//
//        plex_f_size = 0;
//        plex_p_size = 0;
//
//        for (i = 0; i < v_sizes[l]; i++) {
//            u = sub[l][i];
//            end = cd[u] + od[l][u];
//            for (j = cd[u]; j < end; j++) {
//                v = adj[j];
//                e[old2new[u] * v_sizes[l] + old2new[v]] = 1;
//                e[old2new[v] * v_sizes[l] + old2new[u]] = 1;
//            }
//            if (d[l][u] == v_sizes[l] - 1) plex_f[plex_f_size++] = i;
//            else {
//                plex_loc[i] = plex_p_size;
//                plex_p[plex_p_size++] = i;
//            }
//        }
//
//        plex_size = plex_p_size;
//
//        for (i = 0; i < v_sizes[l]; i++) {
//            u = sub[l][i];
//            if (d[l][u] == v_sizes[l] - 1) continue;
//
//            for (j = 0; j < v_sizes[l]; j++) {
//                if (e[i * v_sizes[l] + j]) continue;
//
//                plex_d[i][plex_d_num[i]++] = j;
//                assert(plex_d_num[i] <= dis_m);
//            }
//        }
//
//        delete [] old2new;
//        delete [] e;
//
//        for(i = 0; i <= l; i++) {
//            plex_sub = 0;
//            list_in_plex(i, l - i, cliques);
//        }
//
//        return true;
//    }
//
//    return false;
//}

EBBkC_t::EBBkC_t() = default;

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

    G.read_ordered_edges_from_file(file_name);
    G.build_from_G();

    if (type == 0) {
        GetCurTime(&start);
        G.basic_clique_list(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 1) {
        GetCurTime(&start);
        G.improved_clique_list(K, &N);
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }
    else {
        runtime = 0;
    }
//
//    else if (type == 2) {
//        GetCurTime(&start);
//        EBBkC_Graph_t *g = new EBBkC_Graph_t [G.m];
//        for (int i = 0; i < G.m; i++) {
//            G.edge_oriented_branch(g+i, i);
//            g[i].color_DAG();
//            g[i].build_from_DAG(K - 2);
//        }
//        GetCurTime(&end);
//        runtime = GetTime(&start, &end);
//        printf("%.2lf ms\n", runtime);
//        GetCurTime(&start);
//        for (int i = 0; i < G.m; i++) {
//            g[i].improved_clique_list_plus(K - 2, &N);
//        }
//        GetCurTime(&end);
//        runtime = GetTime(&start, &end);
//    }
//    else {
//        runtime = 0;
//    }

    return runtime;
}