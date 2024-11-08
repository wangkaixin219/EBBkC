#include <functional>
#include <numeric>

#include "../util/graph/graph.h"
#include "../util/log/log.h"
#include "../util/containers/boolarray.h"
#include "../util/intersection/set_utils.h"
#include "pkt_support_update_utils.h"
#include "iter_helper.h"

//Process a sub-level in a level using intersection based approach
void PKT_processSubLevel_intersection_handling_skew(graph_t *g, eid_t *curr,
#ifndef BMP_QUEUE
        bool *InCurr,
#else
                                                    BoolArray<word_type> &InCurr,
#endif
                                                    long currTail, int *EdgeSupport, int level, eid_t *next,
#ifndef BMP_QUEUE
        bool *InNext,
#else
                                                    BoolArray<word_type> &InNext,
#endif
                                                    long *nextTail,
#ifndef BMP_PROCESSED
        bool *processed,
#else
                                                    BoolArray<word_type> &processed,
#endif
                                                    Edge *edgeIdtoEdge, eid_t *off_end,
                                                    bool *is_vertex_updated, IterHelper &iter_helper
) {
    const long BUFFER_SIZE_BYTES = 2048;
    const long BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(vid_t);

    eid_t buff[BUFFER_SIZE];
    LocalWriteBuffer<eid_t, long> local_write_buffer(buff, BUFFER_SIZE, next, nextTail);

    eid_t bkt_buff[BUFFER_SIZE];
    LocalWriteBuffer<eid_t, size_t> local_bucket_buf(bkt_buff, BUFFER_SIZE, iter_helper.bucket_buf_,
                                                     &iter_helper.window_bucket_buf_size_);

    SupportUpdater sup_updater(EdgeSupport, InNext, level, local_write_buffer, local_bucket_buf,
                               iter_helper.bucket_level_end_, iter_helper.in_bucket_window_);


#pragma omp single
    for (long i = 0; i < currTail; i++) {
        eid_t e1 = curr[i];
        g->edge_rank[e1] = (++g->__rank);
        g->edge_truss[e1] = level + 3;
    }

#pragma omp for schedule(dynamic, 4)
    for (long i = 0; i < currTail; i++) {
        //process edge <u,v>
        eid_t e1 = curr[i];
        Edge edge = edgeIdtoEdge[e1];

        vid_t u = edge.u;
        vid_t v = edge.v;

        eid_t uStart = g->num_edges[u], uEnd = off_end[u + 1];
        eid_t vStart = g->num_edges[v], vEnd = off_end[v + 1];

        eid_t off_nei_u = uStart, off_nei_v = vStart;

        static thread_local vector<pair<eid_t, eid_t >> intersection_res(1024 * 1024);
        size_t beg = 0;
        if (uEnd - uStart > 0 && vEnd - vStart > 0) {
            SetInterSectionLookup(g, off_nei_u, uEnd, off_nei_v, vEnd, intersection_res, beg);
#ifdef PAPER_FIGURE
            {
                static thread_local int tid = omp_get_thread_num();
                vector<int> ws;
                for (auto i = 0; i < beg; i++) {
                    ws.emplace_back((iter_helper.g)->adj[intersection_res[i].first]);
                }
                stringstream ss;
                ss << "tid: " << tid << ", " << make_pair(u, v) << ":" << ws;
                log_info("%s", ss.str().c_str());
            }
#endif
            for (auto iter = 0u; iter < beg; iter++) {
                std::tie(off_nei_v, off_nei_u) = intersection_res[iter];
                eid_t e2 = g->eid[off_nei_v];  //<v,w>
                eid_t e3 = g->eid[off_nei_u];  //<u,w>

                assert(g->adj[off_nei_u] == g->adj[off_nei_v]);
                eid_t w = g->adj[off_nei_u];

                sup_updater.PeelTriangleByRank(g, e1, e2, e3, processed, InCurr);
                sup_updater.GenerateVSet(g, e1, e2, e3, w, processed, InCurr);
            }
        }
#ifndef NO_SHRINK_GRAPH
        is_vertex_updated[u] = true;
        is_vertex_updated[v] = true;
#endif
    }

    sup_updater.SubmitLocalBufferNext();
    iter_helper.MarkProcessed();
}

/**   Computes the support of each edge in parallel
 *    Computes k-truss in parallel   ****/
int PKT_intersection(graph_t *g, int *&EdgeSupport, Edge *&edgeIdToEdge) {
    IterHelper iter_helper(g, &EdgeSupport, &edgeIdToEdge);
    auto process_functor = [&iter_helper, g](int level) {
        PKT_processSubLevel_intersection_handling_skew(g, iter_helper.curr_, iter_helper.in_curr_,
                                                       iter_helper.curr_tail_,
                                                       *iter_helper.edge_sup_ptr_, level, iter_helper.next_,
                                                       iter_helper.in_next_, &iter_helper.next_tail_,
                                                       iter_helper.processed_, *iter_helper.edge_lst_ptr_,
                                                       iter_helper.off_end_,
                                                       iter_helper.is_vertex_updated_, iter_helper);
    };
    AbstractPKT(g, EdgeSupport, edgeIdToEdge, iter_helper, process_functor);

    return iter_helper.level_size_ + 2;
}

