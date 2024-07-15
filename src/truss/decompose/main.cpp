//#include <cassert>
//#include <chrono>
//#include <omp.h>
//
//#include "../util/graph/graph.h"
//#include "../util/log/log.h"
//#include "../util/timer.h"
//#include "parallel_all_edge_cnc.h"
//#include "../util/reordering/reorder_utils.h"
//#include "iter_helper.h"
//
//int main(int argc, char *argv[]) {
//    setlocale(LC_NUMERIC, "");
//    if (argc < 2) {
//        fprintf(stderr, "%s <Graph file>\n", argv[0]);
//        exit(1);
//    }
//
////    read_env();
//
//    graph_t g;
//
//    //load the graph from file
//    Graph yche_graph(argv[1]);
//    g.adj = yche_graph.edge_dst;
//    g.num_edges = yche_graph.node_off;
//    g.n = yche_graph.nodemax;
//    g.m = yche_graph.edgemax;
//
//    string reorder_method(argv[2]);
//
//    vector<int32_t> new_vid_dict;
//    vector<int32_t> old_vid_dict;
//    ReorderWrapper(g, argv[1], reorder_method, new_vid_dict, old_vid_dict);
//
//    /************   Compute k - truss *****************************************/
//    //edge list array
//    Timer get_eid_timer;
//
//    Edge *edgeIdToEdge = (Edge *) malloc((g.m / 2) * sizeof(Edge));
//    assert(edgeIdToEdge != nullptr);
//    log_info("Malloc Time: %.9lf s", get_eid_timer.elapsed());
//    get_eid_timer.reset();
//
//    //Populate the edge list array
//    getEidAndEdgeList(&g, edgeIdToEdge);
//    log_info("Init Eid Time: %.9lf s", get_eid_timer.elapsed());
//    get_eid_timer.reset();
//
//    int *EdgeSupport = (int *) malloc(g.m / 2 * sizeof(int));
//    assert(EdgeSupport != nullptr);
//
//    auto max_omp_threads = omp_get_max_threads();
//    log_info("Max Threads: %d", max_omp_threads);
//#pragma omp parallel for
//    for (auto i = 0; i < max_omp_threads; i++) {
//        auto avg = g.m / 2 / max_omp_threads;
//        auto iter_beg = avg * i;
//        auto iter_end = (i == max_omp_threads - 1) ? g.m / 2 : avg * (i + 1);
//        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
//    }
//    log_info("Init EdgeSupport Time: %.9lf s", get_eid_timer.elapsed());
//    get_eid_timer.reset();
//
//    Timer global_timer;
//    PKT_intersection(&g, EdgeSupport, edgeIdToEdge);
//
//
//    FILE *fp = fopen("./edge.index", "w");
//    for (eid_t i = 0; i < g.m / 2; i++) {
//        Edge e = edgeIdToEdge[i];
//        assert(g.edge_truss[i] != 0);
//        fprintf(fp, "%d %d %d %d %ld ", e.u, e.v, g.edge_truss[i], g.edge_rank[i], g.v_set[i].size());
//        for (vid_t j = 0; j < g.v_set[i].size(); j++) fprintf(fp, "%d ", g.v_set[i][j]);
//        fprintf(fp, "\n");
//    }
//    fclose(fp);
//
//    //Free memory
//    free_graph(&g);
//    free(edgeIdToEdge);
//    free(EdgeSupport);
//
//    return 0;
//}
