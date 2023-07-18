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
        int val = -1;
    };

private:
    const size_t TAB_SIZE = 0x7fffff;
    const size_t MAX_COLL = 100;
    int* table_size;
    HashItem_t** table;

public:
    HashMap_t();
    ~HashMap_t();

    int exist(const Edge_t& key);
    void insert(const Edge_t& key, int val);
    void remove(const Edge_t& key);
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

#endif //DESCOL_DEFS_H
