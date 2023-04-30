//
// Created by Wang Kaixin on 2022/6/8.
//

#ifndef DESCOL_DEFS_H
#define DESCOL_DEFS_H

#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <functional>

#define N_EDGES 200000000
#define N_NODES 10000000

using namespace std;

class Edge_t {
public:
    int s = 0;
    int t = 0;

    Edge_t();
    Edge_t(int s, int t, bool directed);

    bool operator<(const Edge_t& e) const;
    bool operator==(const Edge_t& e) const;

    struct Hash_Edge_t {
        size_t operator()(const Edge_t& e) const;
    };
};

class HashMap_t {

    struct HashItem_t {
        Edge_t key;
        int val;
    };

private:
    const size_t TAB_SIZE = 0x7fffffff;
    const size_t MAX_COLL = 100;
    int* table_size;
    HashItem_t** table;

public:
    HashMap_t();
    ~HashMap_t();

    int exist(const Edge_t& key);
    void insert(const Edge_t& key, const int val);

};

class KeyVal_t {
public:
    int key = 0;
    int val = 0;

    KeyVal_t();
    KeyVal_t(int key, int val);
    bool operator<(const KeyVal_t& kv) const;
};

class Heap_t {
private:
    unsigned n{};
    KeyVal_t *kv_list{};
    int *pt{};

    void swap(unsigned i, unsigned j);
    void bubble_up(unsigned i);
    void bubble_down(unsigned i);

public:
    Heap_t();
    ~Heap_t();

    bool empty();
    void insert(KeyVal_t kv);
    KeyVal_t pop();
    void update(unsigned key);
    KeyVal_t min_element();

    void make_heap(const int *v_list, unsigned v_size);
    void release_heap();
};

void clean_edges(const char* r_file_name, const char* w_file_name);

void GetCurTime(struct rusage* curTime);

double GetTime(struct rusage* start, struct rusage* end);

//
//template <class T>
//void hash_combine(std::size_t& seed, const T& v);

#endif //DESCOL_DEFS_H


/*
 *
//    void er_model(const char* write, long n, double p);
//    void ff_model(const char* write, long n, double p);
void graph_t::er_model(const char* write, long n, double p) {
    int s, t;
    long long k;
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> dist(0, 1);

    FILE *fp = fopen(write, "w");
    this->n = n;
    this->m = 0;

    k = 1 + log(dist(gen)) / log(1 - p);

    while (log(k) <= 2 * log(n)) {
        s = (k - 1) / n;
        t = (k - 1) % n;

        if (s < t) {
            this->m++;
            fprintf(fp, "%u %u\n", s, t);
        }

        k += 1 + log(dist(gen)) / log(1 - p);
    }

    fclose(fp);
}

void graph_t::ff_model(const char* write, long n, double p) {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> random(0, 1);
    map<long, vector<long>> V;

    queue<long> vertex_q;
    map<long, bool> visited;

//    edge_t *E = (edge_t*) malloc(N_EDGES * sizeof(edge_t));
    this->n = n;
    this->m = 0;
    long u, v, ambassador, index;
    double q;

//    string filename = "FF" + to_string((int)(n / 1e6)) + "M.edges";
    FILE *fp = fopen(write, "w");

    for (int i = 1; i < n; i++) {
        fprintf(stdout, "%.2lf%%\r", (double) i / n * 100);
        fflush(stdout);
        uniform_int_distribution<> random_vertex(0, i-1);
        ambassador = random_vertex(gen);
        vertex_q.push(ambassador);
        visited[i] = true;
        visited[ambassador] = true;

        while (!vertex_q.empty()) {
            u = vertex_q.front();
            vector<long> neighbors (V[u]);

            while (true) {
                q = random(gen);
                if (q > p || neighbors.empty()) break;

                uniform_int_distribution<> random_neighbor(0, neighbors.size() - 1);
                index = random_neighbor(gen);
                v = neighbors[index];
                neighbors[index] = neighbors[neighbors.size() - 1];
                neighbors.pop_back();
                if (!visited[v]) {
                    vertex_q.push(v);
                    visited[v] = true;
                }
            }

            V[i].push_back(u);
            V[u].push_back(i);

            fprintf(fp, "%lu %u\n", u, i);

            this->m++;
            vertex_q.pop();
        }
        visited.clear();
    }

    fclose(fp);

}
*/