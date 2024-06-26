#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <set>
#include <algorithm>
#include <string>
#include <omp.h>
#include <cassert>
#include "set_operation.h"
#include "def.h"
#include "edge_oriented.h"


using namespace std;
int K, L;
unsigned long long N = 0;

int main(int argc, char** argv) {
    double runtime;
    string act(argv[1]);

    if (act == "p") {                           // pre-process
        printf("Pre-process %s\n", argv[2]);
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".edges" && suffix != ".mtx") exit(0);
        string prefix = src_filename.substr(0, src_filename.find_last_of('.'));
        string clean_filename = prefix + ".clean";
        clean_edges(src_filename.c_str(), clean_filename.c_str());
        string index_filename = prefix + ".index";
        runtime = EBBkC_t::truss_order(clean_filename.c_str(), index_filename.c_str());
        printf("Pre-processed in %.2lf ms\n\n", runtime);
    }

    else if (act == "e") {                      // EBBkC+ET
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);

        K = atoi(argv[3]);
        L = atoi(argv[4]);

        runtime = EBBkC_t::list_k_clique(argv[2]);
        printf("Number of %u-cliques: %llu\n", K, N);
        printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime);
    }

    else if (act == "ep") {                      // EBBkC+ET (parallel)
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);

        K = atoi(argv[3]);
        L = atoi(argv[4]);

        omp_set_num_threads(atoi(argv[5]));

        runtime = EBBkC_t::list_k_clique_parallel(argv[2]);
        printf("Number of %u-cliques: %llu\n", K, N);
        printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime);
    }

    else {
        printf("Wrong usage.\n");
        exit(0);
    }

    return 0;
}
