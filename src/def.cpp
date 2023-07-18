#include <cassert>
#include "def.h"
#include <set>

size_t h(const Edge_t& e) {
    int s_ = e.s < e.t ? e.s : e.t;
    int t_ = e.s < e.t ? e.t : e.s;
    size_t hash = 1;
    std::hash<int> H;

    hash ^= H(s_) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    hash ^= H(t_) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    return hash;
}

Edge_t::Edge_t() = default;

Edge_t::Edge_t(int s, int t, bool directed) {
    if (directed) {
        this->s = s;
        this->t = t;
    }
    else {
        this->s = s < t ? s : t;
        this->t = s < t ? t : s;
    }
}

bool Edge_t::operator<(const Edge_t &e) const {
    return (this->s < e.s) || (this->s == e.s && this->t < e.t);
}

bool Edge_t::operator==(const Edge_t &e) const {
    return this->s == e.s && this->t == e.t;
}

size_t Edge_t::Hash_Edge_t::operator()(const Edge_t &e) const {
    return h(e);
}

HashMap_t::HashMap_t() {
    table = new HashMap_t::HashItem_t* [TAB_SIZE];
    table_size = new int [TAB_SIZE]();
}

HashMap_t::~HashMap_t() {
    delete [] table;
    delete [] table_size;
}

int HashMap_t::exist(const Edge_t &key) {
    size_t pos = h(key) % TAB_SIZE;
    int i, val;
    Edge_t e;

    if (table_size[pos] > 0) {
        for (i = 0; i < table_size[pos]; i++) {
            e = table[pos][i].key;
            val = table[pos][i].val;
            if ((e.s == key.s && e.t == key.t) || (e.s == key.t && e.t == key.s)) return val;
        }
    }

    return -1;
}

void HashMap_t::insert(const Edge_t &key, int val) {
    size_t pos = h(key) % TAB_SIZE;

    if (table_size[pos] == 0)
        table[pos] = new HashMap_t::HashItem_t [MAX_COLL];

    assert(table_size[pos] < MAX_COLL);     // Enlarge max_coll or redefine hash function h(e)

    table[pos][table_size[pos]].key = key;
    table[pos][table_size[pos]++].val = val;
}

void HashMap_t::remove(const Edge_t &key) {
    size_t pos = h(key) % TAB_SIZE;
    int i;
    Edge_t e;

    if (table_size[pos] > 0) {
        for (i = 0; i < table_size[pos]; i++) {
            e = table[pos][i].key;
            if ((e.s == key.s && e.t == key.t) || (e.s == key.t && e.t == key.s)) {
                table[pos][i] = table[pos][--table_size[pos]];
                if (table_size[pos] == 0) delete [] table[pos];
                break;
            }
        }
        return;
    }

    exit(-1);
}

KeyVal_t::KeyVal_t() = default;

KeyVal_t::KeyVal_t(int key, int val) {
    this->key = key;
    this->val = val;
}

bool KeyVal_t::operator<(const KeyVal_t &kv) const {
    return this->val > kv.val || (this->val == kv.val && this->key < kv.key);
}

Heap_t::Heap_t() = default;
Heap_t::~Heap_t() = default;

void Heap_t::swap(unsigned i, unsigned j) {
    KeyVal_t kv_tmp = this->kv_list[i];
    int pt_tmp = this->pt[kv_tmp.key];
    this->pt[this->kv_list[i].key] = this->pt[this->kv_list[j].key];
    this->kv_list[i] = this->kv_list[j];
    this->pt[this->kv_list[j].key] = pt_tmp;
    this->kv_list[j] = kv_tmp;
}

void Heap_t::bubble_up(unsigned int i) {
    unsigned j = (i - 1) >> 1;
    while (i > 0) {
        if (this->kv_list[j].val > this->kv_list[i].val) {
            this->swap(i, j);
            i = j;
            j = (i - 1) >> 1;
        }
        else break;
    }
}

void Heap_t::bubble_down(unsigned int i) {
    unsigned l = (i << 1) + 1, r = l + 1, j;
    while (l < this->n) {
        j = ((r < this->n) && (this->kv_list[r].val < this->kv_list[l].val)) ? r : l;
        if (this->kv_list[i].val > this->kv_list[j].val) {
            this->swap(i, j);
            i = j;
            l = (i << 1) + 1;
            r = l + 1;
        }
        else break;
    }
}

bool Heap_t::empty() {
    return this->n == 0;
}

void Heap_t::insert(KeyVal_t kv) {
    this->pt[kv.key] = this->n;
    this->kv_list[this->n] = kv;
    this->bubble_up(this->n);
    this->n++;
}

KeyVal_t Heap_t::pop() {
    assert(!this->empty());
    KeyVal_t min = this->kv_list[0];
    this->pt[min.key] = -1;
    this->kv_list[0] = this->kv_list[--(this->n)];
    this->pt[this->kv_list[0].key] = 0;
    this->bubble_down(0);
    return min;
}

KeyVal_t Heap_t::min_element() {
    assert(!this->empty());
    return this->kv_list[0];
}

void Heap_t::update(unsigned int key) {
    int i = this->pt[key];
    if (i != -1) {
        ((this->kv_list[i]).val)--;
        this->bubble_up(i);
    }
}

void Heap_t::make_heap(const int *v_list, unsigned int v_size) {
    unsigned i;
    KeyVal_t kv;

    this->n = 0;
    this->pt = new int [v_size];
    for (i = 0; i < v_size; i++) this->pt[i] = -1;
    this->kv_list = new KeyVal_t [v_size];

    for (i = 0; i < v_size; i++) {
        kv.key = i;
        kv.val = v_list[i];
        this->insert(kv);
    }
}

void Heap_t::release_heap() {
    delete [] this->kv_list;
    delete [] this->pt;
}

void clean_edges(const char *r_file_name, const char *w_file_name) {
    set<Edge_t> E;
    int s, t;
    FILE *f_r = fopen(r_file_name, "r");
    FILE *f_w = fopen(w_file_name, "w");

    while (fscanf(f_r, "%u %u%*[^\n]%*c", &s, &t) == 2) {
        if (s == t) {
            printf("%u Self loop.\n", s);
            continue;
        }
        if (E.find(Edge_t(s, t, false)) != E.end()) {
            printf("(%u, %u) duplicates.\n", s, t);
            continue;
        }

        E.insert(Edge_t(s, t, false));
        fprintf(f_w, "%u %u\n", s, t);
    }

    fclose(f_r);
    fclose(f_w);
}

void GetCurTime(struct rusage* curTime) {
    if (getrusage(RUSAGE_THREAD, curTime) != 0) {
        fprintf(stderr, "The running time info couldn't be collected successfully.\n");
        exit(0);
    }
}

double GetTime(struct rusage* start, struct rusage* end) {
    return ((float)(end->ru_utime.tv_sec - start->ru_utime.tv_sec)) * 1e3 +
           ((float)(end->ru_utime.tv_usec - start->ru_utime.tv_usec)) * 1e-3;
}   // unit: ms

