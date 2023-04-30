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
#include "def.h"
#include "VBBkC.h"
#include "EBBkC.h"


using namespace std;
int K, P;
unsigned long long N = 0;

int main(int argc, char** argv) {
    double runtime;
    string act(argv[1]);

    /*
     * act: ["clean", "index", "ff", "er", "list"]
     * "clean": input is an edge file (.edges or .mtx), output is a cleaned edge file (.clean)
     * "index": input is a cleaned edge file (.clean), output is an index file (.index)
     * "list": input is an index file (.index) and k, output is k-cliques.
     */

    if (act == "clean") {
        printf("Cleaning %s\n", argv[2]);
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".edges" && suffix != ".mtx") exit(0);
        string prefix = src_filename.substr(0, src_filename.find_last_of('.'));
        string tgt_filename = prefix + ".clean";
        clean_edges(src_filename.c_str(), tgt_filename.c_str());
    }

    else if (act == "index") {
        EBBkC_t frame;
        printf("Indexing %s\n", argv[2]);
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".clean") exit(0);
        string prefix = src_filename.substr(0, src_filename.find_last_of('.'));
        string tgt_filename = prefix + ".index";
        runtime = frame.truss_order(src_filename.c_str(), tgt_filename.c_str());
        printf("Indexed in %.2lf ms\n", runtime);
    }
    else if (act == "v") {
        VBBkC_t frame;
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);
        K = atoi(argv[3]);
        P = atoi(argv[4]);

        omp_set_num_threads(24);

        printf("Iterate over all cliques\n");
        runtime = frame.list_k_clique(argv[2], atoi(argv[5]));
        printf("Number of %u-cliques: %llu\n", K, N);
        printf("Runtime %.2lf ms\n\n", runtime);
    }

    else if (act == "e") {            // serial
        EBBkC_t frame;
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);
        K = atoi(argv[3]);
        P = atoi(argv[4]);

        omp_set_num_threads(24);

        printf("Iterate over all cliques\n");
        runtime = frame.list_k_clique(argv[2], atoi(argv[5]));
        printf("Number of %u-cliques: %llu\n", K, N);
        printf("Runtime %.2lf ms\n\n", runtime);
    }

    else {
        printf("Wrong usage.\n");
        exit(0);
    }

    return 0;
}
