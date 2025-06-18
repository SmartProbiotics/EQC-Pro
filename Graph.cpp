/*
 * Graph.cpp
 *
 * Author: ljchang@outlook.com
 */

#include "Graph.h"
#include "kDefectiveClique_BB.h"
#include "kDefectiveClique_BB_matrix.h"
#include "CTPrune.h"

using namespace std;

Graph::Graph(const char *_dir, const int _K)
{
	dir = string(_dir);
	K = _K;
	n = m = 0;

	pstart = nullptr;
	edges = nullptr;

	kDefectiveClique.clear();

	pend = pend_buf = nullptr;
	tri_cnt = edges_pointer = Qe = nullptr;
	deleted = nullptr;
}

Graph::Graph(const Graph &other)
	: dir(other.dir), n(other.n), m(other.m), K(other.K)
{

	pstart = new ept[n + 1];
	edges = new ui[m];
	memcpy(pstart, other.pstart, (n + 1) * sizeof(ept));
	memcpy(edges, other.edges, m * sizeof(ui));

	if (other.pend)
	{
		pend = new ept[n];
		memcpy(pend, other.pend, n * sizeof(ept));
	}
	else
	{
		pend = nullptr;
	}

	if (other.pend_buf)
	{
		pend_buf = new ept[n];
		memcpy(pend_buf, other.pend_buf, n * sizeof(ept));
	}
	else
	{
		pend_buf = nullptr;
	}

	if (other.tri_cnt)
	{
		tri_cnt = new ui[n];
		memcpy(tri_cnt, other.tri_cnt, n * sizeof(ui));
	}
	else
	{
		tri_cnt = nullptr;
	}

	if (other.edges_pointer)
	{
		edges_pointer = new ui[m];
		memcpy(edges_pointer, other.edges_pointer, m * sizeof(ui));
	}
	else
	{
		edges_pointer = nullptr;
	}

	if (other.Qe)
	{
		Qe = new ui[m];
		memcpy(Qe, other.Qe, m * sizeof(ui));
	}
	else
	{
		Qe = nullptr;
	}

	if (other.deleted)
	{
		deleted = new char[n];
		memcpy(deleted, other.deleted, n * sizeof(char));
	}
	else
	{
		deleted = nullptr;
	}
}

Graph::Graph()
{
	K = 0;
	n = m = 0;
	pstart = nullptr;
	edges = nullptr;

	kDefectiveClique.clear();

	pend = pend_buf = nullptr;
	tri_cnt = edges_pointer = Qe = nullptr;
	deleted = nullptr;
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != nullptr)
	{
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if (tri_cnt != nullptr)
	{
		delete[] tri_cnt;
		tri_cnt = nullptr;
	}
	if (edges_pointer != nullptr)
	{
		delete[] edges_pointer;
		edges_pointer = nullptr;
	}
	if (Qe != nullptr)
	{
		delete[] Qe;
		Qe = nullptr;
	}
	if (deleted != nullptr)
	{
		delete[] deleted;
		deleted = nullptr;
	}
}

void Graph::read_graph_binary()
{
	FILE *f = Utility::open_file(dir.c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(int), 1, f);
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	// printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);
	if (pstart == nullptr)
		pstart = new ept[n + 1];
	if (edges == nullptr)
		edges = new ui[m];

	pstart[0] = 0;
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] > 0)
		{
			fread(edges + pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			ui *buff = edges + pstart[i];
			sort(buff, buff + degree[i]);
			ui idx = 0;
			for (ui j = 0; j < degree[i]; j++)
			{
				if (buff[j] >= n)
					printf("vertex id %u wrong\n", buff[j]);
				if (buff[j] == i || (j > 0 && buff[j] == buff[j - 1]))
					continue;
				buff[idx++] = buff[j];
			}
			degree[i] = idx;
		}

		pstart[i + 1] = pstart[i] + degree[i];
	}

	fclose(f);
	delete[] degree;
}

void Graph::output_one_kDefectiveClique()
{
	FILE *fout = Utility::open_file("kDefectiveClique.txt", "w");
	fprintf(fout, "%lu\n", kDefectiveClique.size());
	sort(kDefectiveClique.begin(), kDefectiveClique.end());
	for (ui i = 0; i < kDefectiveClique.size(); i++)
		fprintf(fout, " %u", kDefectiveClique[i]);
	fprintf(fout, "\n");
	fclose(fout);
}

void Graph::verify_kDefectiveClique()
{
	char *vis = new char[n];
	memset(vis, 0, sizeof(char) * n);

	FILE *fin = Utility::open_file("kDefectiveClique.txt", "r");

	ui max_kDefectiveClique_size = n, kDefectiveClique_n, idx = 0;
	char ok = 1;
	while (fscanf(fin, "%u", &kDefectiveClique_n) == 1)
	{
		++idx;
		if (max_kDefectiveClique_size == n)
		{
			max_kDefectiveClique_size = kDefectiveClique_n;
			printf("max kDefectiveClique size: %u\n", max_kDefectiveClique_size);
		}
		if (kDefectiveClique_n != max_kDefectiveClique_size)
			printf("!!! WA kDefectiveClique size: %u!\n", kDefectiveClique_n);
		vector<ui> kDefectiveClique;
		for (ui i = 0; i < kDefectiveClique_n; i++)
		{
			ui tmp;
			fscanf(fin, "%u", &tmp);
			kDefectiveClique.pb(tmp);
		}

		for (ui i = 0; i < kDefectiveClique.size(); i++)
		{
			if (vis[kDefectiveClique[i]])
			{
				printf("WA kDefectiveClique! Duplicate vertex: %u\n", idx);
				ok = 0;
				break;
			}
			vis[kDefectiveClique[i]] = 1;
		}
		ui cnt = 0;
		for (ui i = 0; i < kDefectiveClique.size(); i++)
		{
			ui d = 0;
			for (ui j = pstart[kDefectiveClique[i]]; j < pstart[kDefectiveClique[i] + 1]; j++)
				if (vis[edges[j]])
					++d;
			assert(d <= kDefectiveClique.size() - 1);
			cnt += kDefectiveClique.size() - 1 - d;
		}
		if (cnt > K * 2)
		{
			ok = 0;
			printf("WA kDefectiveClique! Missing total %u edges\n", cnt / 2);
		}
		for (ui i = 0; i < kDefectiveClique.size(); i++)
			vis[kDefectiveClique[i]] = 0;
	}
	if (ok)
		printf("Correct kDefectiveClique!\n");
	fclose(fin);

	delete[] vis;
}

void Graph::kDefectiveClique_degen()
{
	Timer t;

	kDefectiveClique.clear();

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;

	if (kDefectiveClique.size() < UB)
		printf("\tHeuristic kDefectiveClique Size: %lu, UB: %u, Time: %s (microseconds)\n", kDefectiveClique.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
	else
		printf("\tMaximum kDefectiveClique Size: %lu, Total Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

int Graph::kDefectiveClique_exact(int lb, int ub)
{
	Timer t;
	kDefectiveClique.clear();

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);
	// UB=ub<UB?ub:UB;
	cout << "heuristic solution: " << kDefectiveClique.size() << "   " << "lb: " << lb << endl;
	if (kDefectiveClique.size() < lb)
	{
		cout << "virtual lb :" << kDefectiveClique.size() << "->" << lb << endl;
		kDefectiveClique.resize(lb);
	}
	if (kDefectiveClique.size() < UB)
	{
		ui old_size = kDefectiveClique.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];
		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);

		if (kDefectiveClique.size() > K + 1)
		{
			CTPrune::truss_pruning(n, m, kDefectiveClique.size() - K - 1, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
			// truss_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
		}
		ego_degen(n, m, peel_sequence, pstart, edges, degree, rid, vis, heap, true);
		bool two_stage = (kDefectiveClique.size() >= K + 1);
		// bool two_stage=false;

		if (kDefectiveClique.size() > old_size)
		{
			old_size = kDefectiveClique.size();
			for (ui i = 0; i < kDefectiveClique.size(); i++)
			{
				assert(kDefectiveClique[i] < n);
				kDefectiveClique[i] = out_mapping[kDefectiveClique[i]];
			}

			if (kDefectiveClique.size() > K + 1)
			{
				CTPrune::truss_pruning(n, m, kDefectiveClique.size() - K - 1, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
				// truss_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
			}
		}

		Timer tt;

		if (two_stage)
		{
			vector<ui> ids;
			vector<pair<ui, ui>> vp;

			ui *peel_sequence_rid = core;
			for (ui i = 0; i < n; i++)
				peel_sequence_rid[peel_sequence[i]] = i;

			memset(vis, 0, sizeof(char) * n);

			kDefectiveClique_BB_matrix *kDefectiveClique_solver_m = new kDefectiveClique_BB_matrix();
			kDefectiveClique_solver_m->allocate_memory(n, K);

			ui search_cnt = 0;
			double min_density = 1, total_density = 0;

			reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
			for (ui i = n; i > 0 && kDefectiveClique.size() < UB; i--)
			{
				ui u = peel_sequence[i - 1];
				if (pend[u] - pstart[u] + K < kDefectiveClique.size() || n - i < kDefectiveClique.size())
					continue;

				if (kDefectiveClique.size() > K + 1)
					extract_subgraph_with_prune(u, peel_sequence_rid, degree, ids, rid, vp, vis, pstart, pend, edges);
				else
					extract_subgraph(u, peel_sequence_rid, ids, rid, vp, vis, pstart, pend, edges);

				if (ids.empty())
					continue;

				double density = vp.size() * 2 / (double)ids.size() / (ids.size() - 1);
				++search_cnt;
				total_density += density;
				if (density < min_density)
					min_density = density;

				kDefectiveClique_solver_m->load_graph(ids.size(), vp);
				ui t_old_size = kDefectiveClique.size();
				kDefectiveClique_solver_m->kDefectiveClique(K, UB, kDefectiveClique, true);
				if (kDefectiveClique.size() > t_old_size)
				{
					for (ui j = 0; j < kDefectiveClique.size(); j++)
						kDefectiveClique[j] = ids[kDefectiveClique[j]];
				}
			}
			delete kDefectiveClique_solver_m;

			if (search_cnt == 0)
				printf("search_cnt: 0, ave_density: 1, min_density: 1\n");
			else
				printf("search_cnt: %u, ave_density: %.5lf, min_density: %.5lf\n", search_cnt, total_density / search_cnt, min_density);
		}

		if (n > kDefectiveClique.size() && (!two_stage || kDefectiveClique.size() < K + 1))
		{
			kDefectiveClique_BB *kDefectiveClique_solver = new kDefectiveClique_BB();
			kDefectiveClique_solver->load_graph(n, pstart, pstart + 1, edges);
			kDefectiveClique_solver->kDefectiveClique(K, UB, kDefectiveClique);
			delete kDefectiveClique_solver;
		}

		if (kDefectiveClique.size() > old_size)
		{
			old_size = kDefectiveClique.size();
			for (ui i = 0; i < kDefectiveClique.size(); i++)
			{
				assert(kDefectiveClique[i] < n);
				kDefectiveClique[i] = out_mapping[kDefectiveClique[i]];
			}
		}

		delete[] out_mapping;
		delete[] rid;

		// printf("*** Search time: %s\n", Utility::integer_to_string(tt.elapsed()).c_str());
	}

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;
	return kDefectiveClique.size();
	printf("\tMaximum kDefectiveClique Size: %lu, Total Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}
	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		for (ui j = pstart[u]; j < pend[u] && rid[edges[j]] > rid[u]; j++)
		{
			ui v = edges[j];
			edges[pend[v]++] = u;
			assert(pend[v] <= pstart[v + 1]);
		}
	}
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend[i] == pstart[i + 1]);
#endif
	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		while (end < pstart[i + 1] && rid[edges[end]] > rid[i])
			++end;
	}
}

void Graph::write_subgraph(ui n, const vector<pair<ui, ui>> &edge_list)
{
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for (ui i = 0; i < edge_list.size(); i++)
		fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::extract_subgraph(ui u, const ui *p_rid, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	rid[u] = 0;
	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 1;
		rid[edges[i]] = ids.size() - 1;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);
	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			ui w = edges[j];
			if (exists[w])
				continue;
			ids.push_back(w);
			exists[w] = 1;
			rid[w] = ids.size() - 1;
		}
	}
	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}
	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_with_prune(ui u, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(!exists[i]);
#endif

	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);

	// Utility::print_array("ids1", ids.data(), 0, ids.size(), 0);

	ui *Q = rid;
	ui Q_n = 0;
	for (ui i = 1; i < ids.size(); i++)
	{
		ui v = ids[i];
		degree[v] = 0;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] >= p_rid[u]; j++)
		{
			if (exists[edges[j]])
				++degree[v];
		}
		if (degree[v] + K < kDefectiveClique.size())
			Q[Q_n++] = v;
	}
	for (ui i = 0; i < Q_n; i++)
	{
		ui v = Q[i];
		exists[v] = 3;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] >= p_rid[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if (degree[edges[j]] + K == kDefectiveClique.size())
					Q[Q_n++] = edges[j];
				--degree[edges[j]];
			}
	}
	assert(Q_n < ids.size());
	if (ids.size() - Q_n + K <= kDefectiveClique.size())
	{
		for (ui i = 0; i < ids.size(); i++)
			exists[ids[i]] = 0;
		ids.clear();
		return;
	}

	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
		if (exists[ids[i]] == 2)
		{
			ui v = ids[i];
			for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] >= p_rid[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids.push_back(edges[j]);
					exists[edges[j]] = 1;
					degree[edges[j]] = 1;
				}
				else
					++degree[edges[j]];
			}
		}

	ui new_size = 1;
	for (ui i = 1; i < old_size; i++)
	{
		if (exists[ids[i]] == 3)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);
	for (ui i = old_size; i < ids.size(); i++)
	{
		if (degree[ids[i]] + K < kDefectiveClique.size())
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for (ui i = 0; i < ids.size(); i++)
		rid[ids[i]] = i;
	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
			{
				assert(rid[v] < ids.size() && rid[edges[j]] < ids.size());
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
			}
	}
	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output)
{
	Timer t;

	ui threshold = (kDefectiveClique.size() > K ? kDefectiveClique.size() - K : 0); // all vertices with degree < threshold can be pruned

	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			peel_sequence[queue_n++] = i;
	for (ui i = 0; i < queue_n; i++)
	{
		ui u = peel_sequence[i];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
			{
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
			}
	}
	ui UB = n;
	if (queue_n == n)
		UB = kDefectiveClique.size();

	ept total_edges = 0;
	memset(vis, 0, sizeof(char) * n);
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
		{
			peel_sequence[queue_n + (new_size++)] = i;
			total_edges += degree[i];
		}
		else
		{
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	if (new_size != 0)
	{
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		UB = 0;
		for (ui i = 0; i < new_size; i++)
		{
			ui u, key;
			heap->pop_min(u, key);
			if (key > max_core)
				max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + K + 1;
			if (new_size - i < t_UB)
				t_UB = new_size - i;
			if (t_UB > UB)
				UB = t_UB;

			if (idx == n && (total_edges + 2 * K) / (new_size - i) >= new_size - i - 1)
				idx = i;
			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 0)
				{
					total_edges -= 2;
					heap->decrement(edges[j], 1);
				}
		}

		if (output)
			printf("*** Degeneracy kDefectiveClique size: %u, max_core: %u, UB: %u, Time: %s (microseconds)\n", new_size - idx, max_core, UB, Utility::integer_to_string(t.elapsed()).c_str());

		if (new_size - idx > kDefectiveClique.size())
		{
			kDefectiveClique.clear();
			for (ui i = idx; i < new_size; i++)
				kDefectiveClique.pb(peel_sequence[queue_n + i]);
			if (!output)
				printf("Degen finds a kDefectiveClique of size: %u\n", new_size - idx);
		}
	}

	return UB;
}

void Graph::ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output)
{
	Timer t;
	if (pend == nullptr)
		pend = new ept[n + 1];
	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);

	if (pend_buf == nullptr)
		pend_buf = new ept[n + 1];
	if (edges_pointer == nullptr)
		edges_pointer = new ui[m];
#ifdef NHEU
	ui *pstart_s = pend_buf;
	ui *pend_s = rid;
	ui *edges_s = edges_pointer;

	vector<ui> Q;
	vector<ui> vs;
	memset(vis, 0, sizeof(char) * n);
	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		if (pend[u] - pstart[u] < kDefectiveClique.size())
			continue;
		// if(core[u]+K < kDefectiveClique.size()) break;

		vs.clear();
		for (ui j = pstart[u]; j < pend[u]; j++)
		{
			vs.push_back(edges[j]);
			vis[edges[j]] = 1;
			degree[edges[j]] = 0;
		}
		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					++degree[vs[j]];
					++degree[edges[k]];
				}
		pend_s[vs[0]] = pstart_s[vs[0]] = 0;
		for (ui j = 1; j < vs.size(); j++)
			pend_s[vs[j]] = pstart_s[vs[j]] = pstart_s[vs[j - 1]] + degree[vs[j - 1]];
		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					edges_s[pend_s[vs[j]]++] = edges[k];
					edges_s[pend_s[edges[k]]++] = vs[j];
				}

		ui threshold = (kDefectiveClique.size() > K + 1 ? kDefectiveClique.size() - K - 1 : 0); // all vertices with degree < threshold can be pruned
		Q.clear();
		for (ui j = 0; j < vs.size(); j++)
			if (degree[vs[j]] < threshold)
			{
				Q.push_back(vs[j]);
				vis[vs[j]] = 0;
			}
		for (ui j = 0; j < Q.size(); j++)
			for (ui k = pstart_s[Q[j]]; k < pend_s[Q[j]]; k++)
				if (vis[edges_s[k]])
				{
					if ((degree[edges_s[k]]--) == threshold)
					{
						Q.push_back(edges_s[k]);
						vis[edges_s[k]] = 0;
					}
				}
		ui cnt = 0, total_edges = 0;
		for (ui j = 0; j < vs.size(); j++)
			if (vis[vs[j]])
			{
				total_edges += degree[vs[j]];
				vs[cnt++] = vs[j];
			}
		assert(cnt + Q.size() == vs.size());
		vs.resize(cnt);
		if (cnt == 0)
			continue;

		heap->init(vs.size(), vs.size() - 1, vs.data(), degree);
		bool found = false;
		for (ui ii = 0; ii < vs.size(); ii++)
		{
			ui v, key;
			heap->pop_min(v, key);
			if (found)
			{
				kDefectiveClique.push_back(v);
				continue;
			}

			if (vs.size() - ii + 1 <= kDefectiveClique.size())
				break;

			if ((total_edges + 2 * K) / (vs.size() - ii) >= vs.size() - ii - 1)
			{
				kDefectiveClique.clear();
				kDefectiveClique.push_back(u);
				kDefectiveClique.push_back(v);
				found = true;
				continue;
			}

			vis[v] = 0;
			for (ept j = pstart_s[v]; j < pend_s[v]; j++)
				if (vis[edges_s[j]])
				{
					total_edges -= 2;
					heap->decrement(edges_s[j], 1);
				}
		}
		for (ui j = 0; j < vs.size(); j++)
			vis[vs[j]] = 0;
	}
#endif
	for (ui i = 0; i < n; i++)
		pend_buf[i] = pend[i];
	for (ui i = 0; i < n; i++)
		for (ept j = pstart[i]; j < pend[i]; j++)
			edges[pend_buf[edges[j]]++] = i;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend_buf[i] == pstart[i + 1]);
#endif

	if (output)
		printf("*** EGo-Degen kDefectiveClique size: %lu, Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// in_mapping and out_mapping can be the same array
void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output)
{
	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] + K >= kDefectiveClique.size())
		{
			rid[i] = cnt;
			if (in_mapping == nullptr)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] + K >= kDefectiveClique.size())
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] + K >= kDefectiveClique.size())
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		// printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + K < kDefectiveClique.size());
		assert(cnt == 0 || core[peel_sequence[n - cnt]] + K >= kDefectiveClique.size());
		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			core[i] = core[out_mapping[i]];
		}

		if (pos > 0 && pos < m / 2)
		{
			ept *pstart_new = new ept[cnt + 1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
			memcpy(edges_new, edges, sizeof(ui) * pos);
			delete[] pstart;
			pstart = pstart_new;
			delete[] edges;
			edges = edges_new;
		}

		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

void Graph::truss_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *out_mapping, ui *rid, ept *pstart, ui *edges, ui *degree, bool output)
{
	if (pend == nullptr)
		pend = new ept[n + 1];
	if (tri_cnt == nullptr)
		tri_cnt = new ui[m];
	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);
	oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

	while (n && remove_and_shrink_oriented_tri(n, m, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree))
	{
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
	}

	if (pend_buf == nullptr)
		pend_buf = new ept[n + 1];
	if (edges_pointer == nullptr)
		edges_pointer = new ui[m];
	reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid);

	if (Qe == nullptr)
		Qe = new ui[m];
	for (ui i = 0; i < n; i++)
	{
		pend[i] = pstart[i + 1];
		degree[i] = pstart[i + 1] - pstart[i];
	}
	if (deleted == nullptr)
		deleted = new char[m];
	memset(deleted, 0, sizeof(char) * m);
	truss_peeling(Qe, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] > 0)
		{
			assert(degree[i] + K + 1 > kDefectiveClique.size());
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
	ui t_cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[peel_sequence[i]] > 0)
			peel_sequence[t_cnt++] = rid[peel_sequence[i]];
	assert(t_cnt == cnt);
	ui pos = 0;
	cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] > 0)
		{
			ui start = pstart[i];
			pstart[cnt] = pos;
			for (ui j = start; j < pend[i]; j++)
				if (!deleted[j])
				{
					assert(degree[edges[j]] > 0);
					edges[pos++] = rid[edges[j]];
				}
			++cnt;
		}
	pstart[cnt] = m = pos;
	n = cnt;

	if (output)
		printf("*** After truss shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

// orient graph
void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for (int i = 0; i < n; i++)
		sum += pend[i] - pstart[i];
	assert(sum * 2 == m);
#endif
}

// oriented triangle counting
void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
{
	memset(adj, 0, sizeof(ui) * n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui) * m);
	for (ui u = 0; u < n; u++)
	{
		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = j + 1;

		for (ept j = pstart[u]; j < pend[u]; j++)
		{
			ui v = edges[j];
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (adj[edges[k]])
				{
					++tri_cnt[j];
					++tri_cnt[k];
					++tri_cnt[adj[edges[k]] - 1];
					++cnt;
				}
		}

		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = 0;
	}

#ifndef NDEBUG
	// printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

bool Graph::remove_and_shrink_oriented_tri(ui &n, ui &m, ui *out_mapping, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *rid, ui *degree)
{
	// printf("begin\n");
	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];
	ept removed_edges = 0;
	for (ui i = 0; i < n; i++)
		for (ui j = pstart[i]; j < pend[i]; j++)
			if (tri_cnt[j] + K + 2 <= kDefectiveClique.size())
			{
				--degree[i];
				--degree[edges[j]];
				++removed_edges;
			}

	// printf("removed_edges: %u, m: %u\n", removed_edges, m/2);

	if (removed_edges <= m / 4)
	{
		// printf("finish\n");
		return false;
	}

	// printf("here1\n");

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] > 0)
		{
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
	ui t_cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[peel_sequence[i]] > 0)
			peel_sequence[t_cnt++] = rid[peel_sequence[i]];
	assert(t_cnt == cnt);

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
	{
		ui cnt = 0;
		for (ui j = pstart[i]; j < pend[i]; j++)
		{
			if (tri_cnt[j] + K + 2 > kDefectiveClique.size())
				++cnt;
			assert(edges[j] < n);
		}
		assert(cnt <= degree[i]);
	}
#endif

	// printf("here2\n");

	ui pos = 0;
	cnt = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] > 0)
		{
			ui start = pstart[i];
			pstart[cnt] = pos;
			for (ui j = start; j < pend[i]; j++)
				if (tri_cnt[j] + K + 2 > kDefectiveClique.size())
					edges[pos++] = rid[edges[j]];
			pend[cnt] = pos;
			// if(pos-pstart[cnt] > degree[i]) printf("pos-pstart[cnt]: %u, degree[i]: %u\n", pos-pstart[cnt], degree[i]);
			assert(pos - pstart[cnt] <= degree[i]);
			pos += degree[i] - (pos - pstart[cnt]);
			// printf("degree[i]: %u, pos: %u, m: %u\n", degree[i], pos, m);
			++cnt;
		}
	pstart[cnt] = m = pos;
	n = cnt;

	// printf("finish\n");
	return true;
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edges_pointer, ui *buf)
{
	for (ui i = 0; i < n; i++)
		pend2[i] = pend[i];
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pstart[i]; j < pend[i]; j++)
		{
			ept &k = pend2[edges[j]];
			edges[k] = i;
			tri_cnt[k] = tri_cnt[j];
			++k;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend2[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pend2[i]; j < pstart[i + 1]; j++)
		{
			ept &k = pend[edges[j]];
			edges[k] = i;
			tri_cnt[k] = tri_cnt[j];
			edges_pointer[k] = j;
			edges_pointer[j] = k;
			++k;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend[i] == pend2[i]);
#endif

	ept *ids = pend2;
	for (ui i = 0; i < n; i++)
	{
		if (pend[i] == pstart[i] || pend[i] == pstart[i + 1])
			continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while (j < pend[i] || k < pstart[i + 1])
		{
			if (k >= pstart[i + 1] || (j < pend[i] && edges[j] < edges[k]))
			{
				ids[pos] = edges[j];
				buf[pos++] = edges_pointer[j++];
			}
			else
			{
				ids[pos] = edges[k];
				buf[pos++] = edges_pointer[k++];
			}
		}
		assert(pos + pstart[i] == pstart[i + 1]);
		for (ept j = 0; j < pos; j++)
		{
			ui idx = pstart[i] + j, k = buf[j];
			edges[idx] = ids[j];
			tri_cnt[idx] = tri_cnt[k];
			edges_pointer[idx] = k;
			edges_pointer[k] = idx;
		}
	}
}

void Graph::compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted, ept *pstart, ept *pend, ui *edges)
{
	ui end = pstart[u];
	for (ui i = pstart[u]; i < pend[u]; i++)
		if (!deleted[i])
		{
			edges[end] = edges[i];
			tri_cnt[end] = tri_cnt[i];
			edges_pointer[end] = edges_pointer[i];
			edges_pointer[edges_pointer[end]] = end;
			deleted[end] = 0;
			++end;
		}
	pend[u] = end;
}

char Graph::find(ui u, ui w, ept b, ept e, char *deleted, ept &idx, ui *edges)
{
	if (b >= e)
		return 0;

	while (b + 1 < e)
	{
		idx = b + (e - b) / 2;
		if (edges[idx] > w)
			e = idx;
		else
			b = idx;
	}

	idx = b;
	if (edges[idx] == w && !deleted[idx])
		return 1;
	return 0;
}

// return the number of peeled edges
void Graph::truss_peeling(ui *Qe, ui *tri_cnt, ui *edges_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	char *exist = deleted;
	for (ui i = 0; i < n; i++)
	{
		assert(pend[i] == pstart[i + 1]);
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			exist[edges[j]] = 1;
		for (ui j = pstart[i] + 1; j < pstart[i + 1]; j++)
			assert(edges[j] > edges[j - 1]);
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
		{
			assert(edges_pointer[edges_pointer[j]] == j);
			assert(tri_cnt[j] == tri_cnt[edges_pointer[j]]);
			assert(edges[j] != i);
			ui cnt = 0, v = edges[j];
			for (ui k = pstart[v]; k < pstart[v + 1]; k++)
				if (exist[edges[k]])
					++cnt;
			assert(cnt == tri_cnt[j]);
		}
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			exist[edges[j]] = 0;
	}
	memset(deleted, 0, sizeof(char) * n);
#endif
	assert(kDefectiveClique.size() >= K + 2);
	ui t_threshold = kDefectiveClique.size() - K - 1;
	ept Qe_n = 0;
	for (ui i = 0; i < n; i++)
		for (ui j = pstart[i]; j < pend[i]; j++)
			if (tri_cnt[j] < t_threshold && edges[j] > i)
			{
				Qe[Qe_n++] = i;
				Qe[Qe_n++] = edges[j];
			}
	for (ept j = 0; j < Qe_n; j += 2)
	{
		ui u = Qe[j], v = Qe[j + 1], idx;
		find(u, v, pstart[u], pend[u], deleted, idx, edges);
		assert(edges[idx] == v);

		ui tri_n = tri_cnt[idx];
		deleted[idx] = deleted[edges_pointer[idx]] = 1;
		--degree[u];
		--degree[v];
		if (pend[u] - pstart[u] > degree[u] * 2)
		{
			compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			// printf("degree[u]: %u, %u\n", degree[u], pend[u]-pstart[u]);
			assert(degree[u] == pend[u] - pstart[u]);
		}
		if (pend[v] - pstart[v] > degree[v] * 2)
		{
			compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			assert(degree[v] == pend[v] - pstart[v]);
		}
		if (pend[u] - pstart[u] < pend[v] - pstart[v])
			swap(u, v);

		if (pend[u] - pstart[u] > (pend[v] - pstart[v]) * 2)
		{	// binary search
			// if(false) {
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (!deleted[k])
				{
					if (tri_n && find(u, edges[k], pstart[u], pend[u], deleted, idx, edges))
					{
						--tri_n;
						--tri_cnt[edges_pointer[idx]];
						if ((tri_cnt[idx]--) == t_threshold)
						{
							Qe[Qe_n++] = u;
							Qe[Qe_n++] = edges[idx];
						}
						--tri_cnt[edges_pointer[k]];
						if ((tri_cnt[k]--) == t_threshold)
						{
							Qe[Qe_n++] = v;
							Qe[Qe_n++] = edges[k];
						}
					}
				}
		}
		else
		{ // sorted_merge
			ept ii = pstart[u], jj = pstart[v];
			while (ii < pend[u] && jj < pend[v])
			{
				if (edges[ii] == edges[jj])
				{
					if (!deleted[ii] && !deleted[jj])
					{
						--tri_n;
						--tri_cnt[edges_pointer[ii]];
						if ((tri_cnt[ii]--) == t_threshold)
						{
							Qe[Qe_n++] = u;
							Qe[Qe_n++] = edges[ii];
						}
						--tri_cnt[edges_pointer[jj]];
						if ((tri_cnt[jj]--) == t_threshold)
						{
							Qe[Qe_n++] = v;
							Qe[Qe_n++] = edges[jj];
						}
					}

					++ii;
					++jj;
				}
				else if (edges[ii] < edges[jj])
					++ii;
				else
					++jj;
			}
		}
		// if(tri_n != 0) printf("tri_n: %u\n", tri_n);
		assert(tri_n == 0);
	}
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		for (ui j = pstart[i]; j < pend[i]; j++)
			assert(deleted[j] || tri_cnt[j] >= t_threshold);
#endif
}
