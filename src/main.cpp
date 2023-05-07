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
#include "vertex_oriented.h"


using namespace std;
int K, L;
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
        printf("Indexing %s\n", argv[2]);
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".clean") exit(0);
        string prefix = src_filename.substr(0, src_filename.find_last_of('.'));
        string tgt_filename = prefix + ".index";
        runtime = EBBkC_t::truss_order(src_filename.c_str(), tgt_filename.c_str());
        printf("Indexed in %.2lf ms\n", runtime);
    }
    else if (act == "v") {
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);

        K = atoi(argv[3]);
        L = atoi(argv[4]);

        int type = atoi(argv[5]);
        assert(type == 0 || type == 1 || type == 2 || type == 3);

        runtime = VBBkC_t::list_k_clique(argv[2], type);
        printf("Number of %u-cliques: %llu\n", K, N);

        if (type == 0) {
            printf("DDegCol Runtime %.2lf ms\n\n", runtime);
        }
        else if (type == 1) {
            printf("DDegree Runtime %.2lf ms\n\n", runtime);
        }
        else if (type == 2) {
            printf("SDegree Runtime %.2lf ms\n\n", runtime);
        }
        else {
            printf("BitCol Runtime %.2lf ms\n\n", runtime);
        }
    }

    else if (act == "e") {            // serial
        string src_filename(argv[2]);
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".index") exit(0);

        K = atoi(argv[3]);
        L = atoi(argv[4]);

        int type = atoi(argv[5]);

        assert(type == 0 || type == 1 || type == 2);

        runtime = EBBkC_t::list_k_clique(argv[2], type);
        printf("Number of %u-cliques: %llu\n", K, N);

        if (type == 0) {
            printf("EBBkC Runtime %.2lf ms\n\n", runtime);
        }
        else if (type == 1) {
            printf("EBBkC+ Runtime %.2lf ms\n\n", runtime);
        }
        else {
            printf("EBBkC++ (%d) Runtime %.2lf ms\n\n", L, runtime);
        }
    }

//    else if (act == "set") {
//        int i;
//        int a[10] = {10, 9, 4, 2, 6}, b[10] = {8, 3, 5, 2, 4, 9};
//        int size_a = 5, size_b = 6;
//        int c[10];
//        sort(a, a + size_a);
//        for (i = 0; i < size_a; i++) {
//            printf("%d ", a[i]);
//        }
//        printf("\n");
//        sort(b, b + size_b);
//        for (i = 0; i < size_b; i++) {
//            printf("%d ", b[i]);
//        }
//        printf("\n");
//
//        int size_c = intersect_simd4x(a, size_a, b, size_b, c);
//
//        for (i = 0; i < size_a; i++) {
//            printf("%d ", a[i]);
//        }
//        printf("\n");
//
//
//        for (i = 0; i < size_b; i++) {
//            printf("%d ", b[i]);
//        }
//        printf("\n");
//
//        for (i = 0; i < size_c; i++) {
//            printf("%d ", c[i]);
//        }
//        printf("\n");
//    }

    else {
        printf("Wrong usage.\n");
        exit(0);
    }

    return 0;
}
