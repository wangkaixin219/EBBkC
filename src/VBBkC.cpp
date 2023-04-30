
#include "VBBkC.h"
#include <cassert>
#include <algorithm>
#include <set>

using namespace std;

extern int K, P;
extern unsigned long long N;

VBBkC_Graph_t::VBBkC_Graph_t() = default;

void VBBkC_Graph_t::color_DAG() {
    int i, j, s, t, u, c, max_deg = 0;
    KeyVal_t kv;

    this->color_num = 1;
    this->color = new int [this->n]();
    auto *_d = new int [this->n]();
    auto *_cd = new int [this->n + 1];
    auto *_adj = new int [2 * this->m];
    auto *used = new bool [this->n + 1]();
    auto *list = new KeyVal_t [this->n];

    for (i = 0; i < this->m; i++) {
        assert(this->edges[i].s < this->n);
        assert(this->edges[i].t < this->n);
        _d[this->edges[i].s]++;
        _d[this->edges[i].t]++;
    }
    _cd[0] = 0;
    for (i = 1; i < this->n + 1; i++) {
        _cd[i] = _cd[i - 1] + _d[i - 1];
        max_deg = (max_deg > _d[i - 1]) ? max_deg : _d[i - 1];
        list[i - 1].key = i - 1;
        list[i - 1].val = _d[i - 1];
        _d[i - 1] = 0;
    }
    for (i = 0; i < this->m; i++) {
        _adj[_cd[this->edges[i].s] + _d[this->edges[i].s]++] = this->edges[i].t;
        _adj[_cd[this->edges[i].t] + _d[this->edges[i].t]++] = this->edges[i].s;
    }

    sort(list, list + this->n);

    for (i = 0; i < this->n; i++) {
        kv = list[i];

        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++) {
            u = _adj[j];
            used[this->color[u]] = true;
        }

        for (c = 1; used[c]; c++) ;
        this->color[kv.key] = c;
        this->color_num = (c > this->color_num) ? c : this->color_num;

        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++) {
            u = _adj[j];
            used[this->color[u]] = false;
        }

    }

    delete [] _d;
    delete [] _cd;
    delete [] _adj;
    delete [] used;
    delete [] list;

    for (i = 0; i < this->m; i++) {
        s = this->edges[i].s;
        t = this->edges[i].t;

        this->edges[i].s = (this->color[s] > this->color[t]) ? s : t;
        this->edges[i].t = (this->color[s] > this->color[t]) ? t : s;
    }
}


void VBBkC_Graph_t::build_from_DAG(int height) {
    int i, max_deg = 0;
    auto *_od = new int [this->n]();
    auto *_d = new int [this->n]();
    auto *_sub = new int [this->n];

    this->depth = height;
    this->cd = new int [this->n + 1];
    this->adj = new int [this->m];
    this->lab = new int [this->n];

    for (i = 0; i < this->m; i++) {
        _od[this->edges[i].s]++;
        _d[this->edges[i].s]++;
        _d[this->edges[i].t]++;
    }
    this->cd[0] = 0;
    for (i = 1; i < this->n + 1; i++) {
        this->cd[i] = this->cd[i - 1] + _od[i - 1];
        max_deg = (_od[i - 1] > max_deg) ? _od[i - 1] : max_deg;
        _od[i - 1] = 0;
        _sub[i - 1] = i - 1;
        this->lab[i - 1] = height;
    }
    for (i = 0; i < this->m; i++) {
        int s = this->edges[i].s, t = this->edges[i].t;
        this->adj[this->cd[s] + _od[s]++] = t;
    }

    this->ns = new int [height + 1];
    this->ns[height] = this->n;

    this->d = new int* [height + 1];
    this->od = new int* [height + 1];
    this->sub = new int* [height + 1];
    for (i = 0; i < height; i++) {
        this->d[i] = new int [this->n];
        this->od[i] = new int [this->n];
        this->sub[i] = new int [max_deg];
    }
    this->d[height] = _d;
    this->od[height] = _od;
    this->sub[height] = _sub;

    // for early-termination
    this->plex_new2old = new int [this->n];
    this->plex_f = new int [this->n];
    this->plex_p = new int [this->n];
    this->plex_d_num = new int [this->n];
    this->plex_d = new int* [this->n];
    for (i = 0; i < this->n; i++) this->plex_d[i] = new int [P + 1];
    this->plex_lev = new int [this->n];
    this->plex_loc = new int [this->n];
}

void VBBkC_Graph_t::list_clique(int l, unsigned long long* cliques) {
    int u, v, w, i, j, k, end;

    if (this->ns[l] < l)    return;

    if (l == 1) {       // for ep
        (*cliques) += this->ns[l];
        return;
    }

    if (l == 2) {
        for (i = 0; i < this->ns[l]; i++) {
            u = this->sub[l][i];
            end = this->cd[u] + this->od[l][u];
            for (j = this->cd[u]; j < end; j++) {
                (*cliques)++;
            }
        }
        return;
    }

    for (i = 0; i < this->ns[l]; i++) {
        u = this->sub[l][i];

        if (this->color[u] < l) continue;

        this->ns[l - 1] = 0;
        end = this->cd[u] + this->od[l][u];
        for (j = this->cd[u]; j < end; j++) {
            v = this->adj[j];
            if (this->lab[v] == l) {
                this->lab[v] = l - 1;
                this->sub[l - 1][this->ns[l - 1]++] = v;
                this->od[l - 1][v] = 0;
                this->d[l - 1][v] = 0;
            }
        }

        this->m = 0;

        for (j = 0; j < this->ns[l - 1]; j++) {
            v = this->sub[l - 1][j];
            end = this->cd[v] + this->od[l][v];
            for (k = this->cd[v]; k < end; k++) {
                w = this->adj[k];
                if (this->lab[w] == l - 1) {
                    this->od[l - 1][v]++;
                    this->d[l - 1][v]++;
                    this->d[l - 1][w]++;
                    this->m++;
                }
                else {
                    this->adj[k--] = this->adj[--end];
                    this->adj[end] = w;
                }
            }
        }

        list_clique(l - 1, cliques);

        for (j = 0; j < this->ns[l - 1]; j++) {
            v = this->sub[l - 1][j];
            this->lab[v] = l;
        }
    }
}

void VBBkC_Graph_t::list_clique_plus(int l, unsigned long long *cliques) {
    int u, v, w, i, j, k, end;

    if (this->ns[l] < l)    return;

    if (l == 2) {
        for (i = 0; i < this->ns[l]; i++) {
            u = this->sub[l][i];
            end = this->cd[u] + this->od[l][u];
            for (j = this->cd[u]; j < end; j++) {
                (*cliques)++;
            }
        }
        return;
    }

    if (l == 1) {
        (*cliques) += this->ns[l];
        return;
    }

    if (list_clique_in_plex(l, cliques)) return;

    for (i = 0; i < this->ns[l]; i++) {
        u = this->sub[l][i];

        if (this->color[u] < l) continue;

        this->ns[l - 1] = 0;
        end = this->cd[u] + this->od[l][u];
        for (j = this->cd[u]; j < end; j++) {
            v = this->adj[j];
            this->lab[v] = l - 1;
            this->sub[l - 1][this->ns[l - 1]++] = v;
            this->od[l - 1][v] = 0;
            this->d[l - 1][v] = 0;

        }

        this->m = 0;

        for (j = 0; j < this->ns[l - 1]; j++) {
            v = this->sub[l - 1][j];
            end = this->cd[v] + this->od[l][v];
            for (k = this->cd[v]; k < end; k++) {
                w = this->adj[k];
                if (this->lab[w] == l - 1) {
                    this->od[l - 1][v]++;
                    this->d[l - 1][v]++;
                    this->d[l - 1][w]++;
                    this->m++;
                }
                else {
                    this->adj[k--] = this->adj[--end];
                    this->adj[end] = w;
                }
            }
        }

        list_clique_plus(l - 1, cliques);

        for (j = 0; j < this->ns[l - 1]; j++) {
            v = this->sub[l - 1][j];
            this->lab[v] = l;
        }
    }
}

void VBBkC_Comb_list(int *list, int list_size, int start, int l, int k, unsigned long long *cliques) {
    if (l == k) {
        (*cliques)++;
        return;
    }

    for (int i = start; i < list_size; i++) {
        VBBkC_Comb_list(list, list_size, i + 1, l + 1, k, cliques);
    }
}

void VBBkC_Graph_t::list_in_plex(int p, int q, unsigned long long *cliques) {
    if (p == 0) {
        VBBkC_Comb_list(this->plex_f, this->plex_f_size, 0, 0, min(q, this->plex_f_size - q), cliques);
        return;
    }

    int i, j, u, v, prev, vis = 0;

    for (i = this->plex_sub; i < this->plex_size && (this->plex_p_size >= p && this->plex_f_size>= q); i++) {
        u = this->plex_p[i];
        if (this->plex_lev[u]) continue;

        for (j = 0; j < this->plex_d_num[u]; j++) {
            v = this->plex_d[u][j];
            if (this->plex_loc[v] >= i && this->plex_lev[v] == 0) {
                this->plex_lev[v] = p;
                this->plex_p_size--;
                assert(this->plex_p_size >= 0);
            }
        }

        this->plex_sub = i + 1;
        list_in_plex(p - 1, q, cliques);

        for (j = 0; j < this->plex_d_num[u]; j++) {
            v = this->plex_d[u][j];
            if (this->plex_loc[v] >= i && this->plex_lev[v] == p) {
                this->plex_lev[v] = 0;
                this->plex_p_size++;
            }
        }

        this->plex_p_size--;
        vis++;
    }
    this->plex_p_size += vis;
}

bool VBBkC_Graph_t::list_clique_in_plex(int l, unsigned long long* cliques) {
    int i, j, u, v, end, dis_m = 0, dis_num;

    if (this->m < this->ns[l] * (this->ns[l] - P) / 2) return false;

    if (this->m == this->ns[l] * (this->ns[l] - 1) / 2) {
        VBBkC_Comb_list(this->sub[l], this->ns[l], 0, 0, min(l, this->ns[l] - l), cliques);
        return true;
    }

    for (i = 0; i < this->ns[l]; i++) {
        u = this->sub[l][i];
        dis_num = this->ns[l] - this->d[l][u];
        dis_m = dis_m > dis_num ? dis_m : dis_num;
        if (dis_m > P) return false;
    }

    if (dis_m >= 2) {
        int *old2new = new int [this->n];
        int *e = new int [this->ns[l] * this->ns[l]]();

        for (i = 0; i < this->ns[l]; i++) {
            u = this->sub[l][i];
            this->plex_new2old[i] = u;
            this->plex_d_num[i] = 0;
            this->plex_lev[i] = 0;
            old2new[u] = i;
        }

        this->plex_f_size = 0;
        this->plex_p_size = 0;

        for (i = 0; i < this->ns[l]; i++) {
            u = this->sub[l][i];
            end = this->cd[u] + this->od[l][u];
            for (j = this->cd[u]; j < end; j++) {
                v = this->adj[j];
                e[old2new[u] * this->ns[l] + old2new[v]] = 1;
                e[old2new[v] * this->ns[l] + old2new[u]] = 1;
            }
            if (this->d[l][u] == this->ns[l] - 1) this->plex_f[this->plex_f_size++] = i;
            else {
                this->plex_loc[i] = this->plex_p_size;
                this->plex_p[this->plex_p_size++] = i;
            }
        }

        this->plex_size = this->plex_p_size;

        for (i = 0; i < this->ns[l]; i++) {
            u = this->sub[l][i];
            if (this->d[l][u] == this->ns[l] - 1) continue;

            for (j = 0; j < this->ns[l]; j++) {
                if (e[i * this->ns[l] + j]) continue;

                this->plex_d[i][this->plex_d_num[i]++] = j;
                assert(this->plex_d_num[i] <= dis_m);
            }
        }

        delete [] old2new;
        delete [] e;

        for(i = 0; i <= l; i++) {
            this->plex_sub = 0;
            list_in_plex(i, l - i, cliques);
        }

        return true;
    }

    return false;
}

VBBkC_Graph_t::~VBBkC_Graph_t() {
    delete [] this->edges;
    delete [] this->ns;
    delete [] this->cd;
    delete [] this->adj;
    delete [] this->rank;
    delete [] this->new2old;
    delete [] this->lab;
    delete [] this->color;

    for (int i = 0; i < this->depth + 1; ++i) {
        if (this->d) delete [] this->d[i];
        if (this->od) delete [] this->od[i];
        if (this->sub) delete [] sub[i];
    }

    delete [] this->d;
    delete [] this->od;
    delete [] this->sub;

    delete [] this->plex_p;
    delete [] this->plex_f;
    for (int i = 0; i < this->n; i++) if(this->plex_d) delete [] this->plex_d[i];
    delete [] this->plex_d;
    delete [] this->plex_d_num;
    delete [] this->plex_lev;
    delete [] this->plex_new2old;
    delete [] this->plex_loc;
}

void VBBkC_Graph_t::read_edges_from_file(const char *file_name) {      // read from .clean file
    FILE *fp;
    if ((fp = fopen(file_name, "r")) == nullptr) {
        printf("Cannot open file %s.\n", file_name);
        exit(0);
    }

    int u, v, i;
    int *old2new = new int [N_NODES];
    for (i = 0; i < N_NODES; i++) old2new[i] = -1;
    this->n = this->m = 0;
    this->edges = new Edge_t [N_EDGES];
    this->new2old = new int [N_NODES];

    while (fscanf(fp, "%u %u%*[^\n]%*c", &u, &v) == 2) {
        if (u > N_NODES || v > N_NODES) {
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
            exit(0);
        }
        if (old2new[u] == -1) {
            old2new[u] = this->n;
            this->new2old[this->n++] = u;
        }
        if (old2new[v] == -1) {
            old2new[v] = this->n;
            this->new2old[this->n++] = v;
        }

        this->edges[this->m++] = Edge_t(old2new[u], old2new[v], false);

        if (this->m > N_EDGES) {
            printf("Enlarge N_EDGES to at least %d.\n", this->m);
            exit(0);
        }
    }

    printf("|V| = %d, |E| = %d\n", this->n, this->m);

    fclose(fp);

    delete [] old2new;
}

void VBBkC_Graph_t::core_DAG() {
    int i, j, s, t;
    KeyVal_t kv;
    Heap_t heap;

    this->core_num = 0;
    this->rank = new int [this->n];
    auto *_d = new int [this->n]();
    auto *_cd = new int [this->n + 1];
    auto *_adj = new int [2 * this->m];

    for (i = 0; i < this->m; i++) {
        _d[this->edges[i].s]++;
        _d[this->edges[i].t]++;
    }
    _cd[0] = 0;
    for (i = 1; i < this->n + 1; i++) {
        _cd[i] = _cd[i - 1] + _d[i - 1];
        _d[i - 1] = 0;
    }
    for (i = 0; i < this->m; i++) {
        _adj[_cd[this->edges[i].s] + _d[this->edges[i].s]++] = this->edges[i].t;
        _adj[_cd[this->edges[i].t] + _d[this->edges[i].t]++] = this->edges[i].s;
    }

    heap.make_heap(_d, this->n);
    for (i = 0; i < this->n; i++) {
        kv = heap.pop();
        this->core_num = (kv.val > this->core_num) ? kv.val : this->core_num;
        this->rank[kv.key] = n - i - 1;
        for (j = _cd[kv.key]; j < _cd[kv.key + 1]; j++)
            heap.update(_adj[j]);
    }
    heap.release_heap();

    delete [] _d;
    delete [] _cd;
    delete [] _adj;

    for (i = 0; i < this->m; i++) {
        s = this->edges[i].s;
        t = this->edges[i].t;

        this->edges[i].s = (this->rank[s] > this->rank[t]) ? s : t;
        this->edges[i].t = (this->rank[s] > this->rank[t]) ? t : s;
    }
    printf("Core number = %d\n", this->core_num);
}

void VBBkC_Graph_t::vertex_oriented_branch(VBBkC_Graph_t *sg, int node) {
    int i, j, u, v, w, end;
    static int *old2new = nullptr, *_lab = nullptr;

#pragma omp threadprivate(old2new, _lab)
    if (!_lab) {
        old2new = new int [this->n];
        _lab = new int [this->n]();
    }

    sg->n = 0;
    sg->new2old = new int [this->core_num];

    u = this->sub[K][node];
    end = this->cd[u] + this->od[K][u];
    for (i = this->cd[u]; i < end; i++) {
        v = this->adj[i];
        _lab[v] = K - 1;
        old2new[v] = sg->n;
        sg->new2old[sg->n++] = v;
    }

    sg->m = 0;
    sg->edges = new Edge_t [sg->n * sg->n];

    for (i = 0; i < sg->n; i++) {
        v = sg->new2old[i];
        end = this->cd[v] + this->od[K][v];
        for (j = this->cd[v]; j < end; j++) {
            w = this->adj[j];

            if (_lab[w] == K - 1) {
                sg->edges[sg->m].s = old2new[v];
                sg->edges[sg->m++].t = old2new[w];
            }
        }
    }

    end = this->cd[u] + this->od[K][u];
    for (i = this->cd[u]; i < end; i++) {
        v = this->adj[i];
        _lab[v] = K;
    }
}

void VBBkC_Graph_t::vertex_oriented_twice_branch(VBBkC_Graph_t *sg, int edge) {
    int i, j, u, v, w, end;
    static int *old2new = nullptr, *_lab = nullptr;

#pragma omp threadprivate(old2new, _lab)
    if (!_lab) {
        old2new = new int [this->n];
        _lab = new int [this->n]();
    }

    sg->n = 0;
    sg->new2old = new int [this->core_num];

    u = this->edges[edge].s;
    v = this->edges[edge].t;

    end = this->cd[u] + this->od[K][u];
    for (i = this->cd[u]; i < end; i++) {
        w = this->adj[i];
        _lab[w] = K - 1;
    }

    end = this->cd[v] + this->od[K][v];
    for (i = this->cd[v]; i < end; i++) {
        w = this->adj[i];
        if (_lab[w] == K - 1) {
            _lab[w] = K - 2;
            old2new[w] = sg->n;
            sg->new2old[sg->n++] = w;
        }
    }

    sg->m = 0;
    sg->edges = new Edge_t [sg->n * sg->n];

    for (i = 0; i < sg->n; i++) {
        v = sg->new2old[i];
        end = this->cd[v] + this->od[K][v];
        for (j = this->cd[v]; j < end; j++) {
            w = this->adj[j];
            if (_lab[w] == K - 2) {
                sg->edges[sg->m].s = old2new[v];
                sg->edges[sg->m++].t = old2new[w];
            }
        }
    }

    end = this->cd[u] + this->od[K][u];
    for (i = this->cd[u]; i < end; i++) {
        w = this->adj[i];
        _lab[w] = K;
    }
}

VBBkC_t::VBBkC_t() = default;

double VBBkC_t::list_k_clique(const char* file_name, int type) {
    double runtime;
    struct rusage start, end;
    VBBkC_Graph_t G, g;

    G.read_edges_from_file(file_name);

    if (type == 0) {
        GetCurTime(&start);
        G.core_DAG();
        G.build_from_DAG(K);
        VBBkC_Graph_t* gg = new VBBkC_Graph_t [G.n];
        for (int i = 0; i < G.n; i++) {
            G.vertex_oriented_branch(gg+i, i);
            gg[i].color_DAG();
            gg[i].build_from_DAG(K - 1);
            gg[i].list_clique(K - 1, &N);
        }
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 1) {
        GetCurTime(&start);
        G.core_DAG();
        G.build_from_DAG(K);
        for (int i = 0; i < G.n; i++) {
            G.vertex_oriented_branch(&g, i);
            g.color_DAG();
            g.build_from_DAG(K - 1);
            g.list_clique_plus(K - 1, &N);
        }
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 2) {
        GetCurTime(&start);
        G.core_DAG();
        G.build_from_DAG(K);
        VBBkC_Graph_t* gg = new VBBkC_Graph_t [G.m];
        for (int i = 0; i < G.m; i++) {
            G.vertex_oriented_twice_branch(gg+i, i);
            gg[i].color_DAG();
            gg[i].build_from_DAG(K - 2);
            gg[i].list_clique(K - 2, &N);
        }
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 3) {
        GetCurTime(&start);
        G.core_DAG();
        G.build_from_DAG(K);
        for (int i = 0; i < G.m; i++) {
            G.vertex_oriented_twice_branch(&g, i);
            g.color_DAG();
            g.build_from_DAG(K - 2);
            g.list_clique_plus(K - 2, &N);
        }
        GetCurTime(&end);
        runtime = GetTime(&start, &end);
    }

    else if (type == 4) {
#pragma omp parallel default(none) private(g, start, end) shared(G, K) reduction(max:runtime) reduction(+:N)
        {
            GetCurTime(&start);
            G.core_DAG();
            G.build_from_DAG(K);
#pragma omp for schedule(dynamic, 1) nowait
            for (int i = 0; i < G.n; i++) {
                G.vertex_oriented_branch(&g, i);
                g.color_DAG();
                g.build_from_DAG(K - 1);
                g.list_clique(K - 1, &N);
            }
            GetCurTime(&end);
            runtime = GetTime(&start, &end);
        }
    }

    else if (type == 5) {
#pragma omp parallel default(none) private(g, start, end) shared(G, K) reduction(max:runtime) reduction(+:N)
        {
            GetCurTime(&start);
            G.core_DAG();
            G.build_from_DAG(K);
#pragma omp for schedule(dynamic, 1) nowait
            for (int i = 0; i < G.n; i++) {
                G.vertex_oriented_branch(&g, i);
                g.color_DAG();
                g.build_from_DAG(K - 1);
                g.list_clique_plus(K - 1, &N);
            }
            GetCurTime(&end);
            runtime = GetTime(&start, &end);
        }
    }

    else if (type == 6) {
#pragma omp parallel default(none) private(g, start, end) shared(G, K) reduction(max:runtime) reduction(+:N)
        {
            GetCurTime(&start);
            G.core_DAG();
            G.build_from_DAG(K);
#pragma omp for schedule(dynamic, 1) nowait
            for (int i = 0; i < G.m; i++) {
                G.vertex_oriented_twice_branch(&g, i);
                g.color_DAG();
                g.build_from_DAG(K - 2);
                g.list_clique(K - 2, &N);
            }
            GetCurTime(&end);
            runtime = GetTime(&start, &end);
        }
    }

    else {
#pragma omp parallel default(none) private(g, start, end) shared(G, K) reduction(max:runtime) reduction(+:N)
        {
            GetCurTime(&start);
            G.core_DAG();
            G.build_from_DAG(K);
#pragma omp for schedule(dynamic, 1) nowait
            for (int i = 0; i < G.m; i++) {
                G.vertex_oriented_twice_branch(&g, i);
                g.color_DAG();
                g.build_from_DAG(K - 2);
                g.list_clique_plus(K - 2, &N);
            }
            GetCurTime(&end);
            runtime = GetTime(&start, &end);
        }
    }

    return runtime;
}
