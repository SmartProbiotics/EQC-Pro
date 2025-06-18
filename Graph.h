/*
 * Graph.h
 *
 * Author: ljchang@outlook.com
 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph {
public:
	std::string dir; //input graph directory
	ui n; //number of nodes of the graph
	ept m; //number of edges of the graph
	ui K; //the value of k in k-plex

	ept *pstart; //offset of neighbors of nodes
	ui *edges; //adjacent ids of edges

	ept *pend;
	ept *pend_buf;
	ui *tri_cnt;
	ui *edges_pointer;
	ui *Qe;
	char *deleted;

	std::vector<ui> kDefectiveClique;


	Graph(const char *_dir, const int _K) ;
	Graph();
	Graph(const Graph&);
	~Graph() ;
	
	void read_graph_binary() ;
	void read_graph() ;

	void output_one_kDefectiveClique() ;
	void verify_kDefectiveClique() ;

	void kDefectiveClique_degen() ;
	int kDefectiveClique_exact(int lb,int ub) ;

private:
	void reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges) ;
	void write_subgraph(ui n, const std::vector<std::pair<ui,ui> > &edge_list) ;
	void extract_subgraph(ui u, const ui *p_rid, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges) ;
	void extract_subgraph_with_prune(ui u, const ui *p_rid, ui *degree, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges) ;

	ui degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output) ;
	void ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output) ;
	void core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output) ;
	void truss_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *out_mapping, ui *rid, ept *pstart, ui *edges, ui *degree, bool output) ;

	void orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid) ;
	void oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj) ;
	bool remove_and_shrink_oriented_tri(ui &n, ui &m, ui *out_mapping, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *rid, ui *degree) ;
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edges_pointer, ui *buf) ;
	void compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted, ept *pstart, ept *pend, ui *edges) ;
	char find(ui u, ui w, ept b, ept e, char *deleted, ept &idx, ui *edges) ;
	void truss_peeling(ui *Qe, ui *tri_cnt, ui *edges_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges) ;
};
#endif
