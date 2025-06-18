/*
 * Author: ljchang@outlook.com
 */

#ifndef _KDEFECTIVE_CLIQUE_BB_MATRIX_
#define _KDEFECTIVE_CLIQUE_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"

class kDefectiveClique_BB_matrix {
private:
	ui n,m,edges_in_s;
	ui max_n;
	char *matrix;
	long long matrix_size;
    ui *removed_level;

    ui *degree;
    ui *degree_in_S;

    ui K, UB;
    ui *best_solution;
    ui best_solution_size;

    ui *SR; // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;

    ui *neighbors;
    ui *nonneighbors;
    ui *buf;
    ui *buf1;
    ui *buf2;
    ui *buf3;
    ui *buf4;
    char *vis;

public:
    kDefectiveClique_BB_matrix() {
    	n = max_n = 0;
    	matrix = NULL;
    	matrix_size = 0;
        degree = degree_in_S = NULL;
        
        best_solution = NULL;
        K = best_solution_size = UB = 0;

        SR = SR_rid = NULL;
        removed_level = NULL;

        neighbors = nonneighbors = NULL;
        buf = buf1 = buf2 = buf3 = buf4 = NULL;
        vis = NULL;
    }

    ~kDefectiveClique_BB_matrix() {
    	if(matrix != NULL) {
    		delete[] matrix;
    		matrix = NULL;
    	}
        if(degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
        if(degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if(best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if(SR != NULL) {
            delete[] SR;
            SR = NULL;
        }
        if(SR_rid != NULL) {
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if(removed_level != NULL) {
        	delete[] removed_level;
        	removed_level = NULL;
        }
        if(neighbors != NULL) {
        	delete[] neighbors;
        	neighbors = NULL;
        }
        if(nonneighbors != NULL) {
        	delete[] nonneighbors;
        	nonneighbors = NULL;
        }
        if(buf != NULL) {
        	delete[] buf;
        	buf = NULL;
        }
        if(buf1 != NULL) {
        	delete[] buf1;
        	buf1 = NULL;
        }
        if(buf2 != NULL) {
        	delete[] buf2;
        	buf2 = NULL;
        }
        if(buf3 != NULL) {
        	delete[] buf3;
        	buf3 = NULL;
        }
        if(buf4 != NULL) {
        	delete[] buf4;
        	buf4 = NULL;
        }
        if(vis != NULL) {
        	delete[] vis;
        	vis = NULL;
        }
    }

    void allocate_memory(ui _max_n, ui _K) {
    	max_n = _max_n; K = _K;

        assert(degree == NULL);
		degree = new ui[max_n+1];
		degree_in_S = new ui[max_n+1];
		best_solution = new ui[max_n+1];
		SR = new ui[max_n+1];
		SR_rid = new ui[max_n+1];
		removed_level = new ui[max_n+1];
		if(max_n > K) {
			neighbors = new ui[max_n+1];
			nonneighbors = new ui[max_n+1];
			buf = new ui[max_n+1];
			buf1 = new ui[max_n+1];
			buf2 = new ui[max_n+1];
			buf3 = new ui[max_n+1];
			buf4 = new ui[max_n+1];
		}
		else {
			neighbors = new ui[K+2];
			nonneighbors = new ui[K+2];
			buf = new ui[K+2];
			buf1 = new ui[K+2];
			buf2 = new ui[K+2];
			buf3 = new ui[K+2];
			buf4 = new ui[K+2];
		}
		vis = new char[max_n+1];
    }

    void load_graph(ui _n, const std::vector<std::pair<ui,ui> > &vp) {
		n = _n; m=vp.size(), edges_in_s=0;
		assert(n > 0&&n <= max_n);

		if(matrix == NULL) {
			matrix_size = n*(long long)n;
			matrix = new char[matrix_size];
		}

		if(((long long)n)*n > matrix_size) {
			do {
				matrix_size *= 2;
			} while(((long long)n)*n > matrix_size);
			delete[] matrix; matrix = new char[matrix_size];
		}

		memset(matrix, 0, sizeof(char)*((long long)n)*n);
		for(ui i = 0;i < n;i ++) degree[i] = 0;
		for(ui i = 0;i < vp.size();i ++) {
			assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			assert(!matrix[a*n+b]);
			degree[a] ++;
			degree[b] ++;
			matrix[a*n + b] = matrix[b*n + a] = 1;
		}
	}

    void kDefectiveClique(ui _K, ui _UB, std::vector<ui> &kDC, bool must_include_0) {
        K = _K;
        UB = _UB;
        best_solution_size = kDC.size();

		memset(vis, 0, sizeof(char)*(n+1));
		memset(degree_in_S, 0, sizeof(ui)*(n+1));
		for(ui i = 0;i < n;i ++) removed_level[i] = n;
		for(ui i = 0;i < n;i ++) SR[i] = SR_rid[i] = i;
		while(!Qv.empty()) Qv.pop();

		for(ui i = 0;i < n;i ++) {
			if(degree[i] + K < best_solution_size) {
				removed_level[i] = 0;
				Qv.push(i);
			}
		}

		ui R_end = n;
		remove_vertices_and_prune(0, R_end, 0);

		// if(R_end != n) printf("Initially pruned %u vertices!\n", n - R_end);

		if(!must_include_0||SR_rid[0] < R_end) BB_search(0, R_end, 1, must_include_0);

        if(best_solution_size > kDC.size()) {
            kDC.clear();
            for(int i = 0;i < best_solution_size;i ++) kDC.push_back(best_solution[i]);
        }
    }

private:
    void compute_a_heuristic_solution_and_prune(ui R_end, ui level) {
    	// the following computes the degeneracy ordering and a heuristic solution
#ifndef NDEBUG
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	for(ui i = 0;i < R_end;i ++) assert(!vis[SR[i]]);
#endif
		ui *core = neighbors;
		ui *t_degree = nonneighbors;
		ui total_edges = 0;
		for(ui i = 0;i < R_end;i ++) {
			t_degree[SR[i]] = degree[SR[i]];
			assert(t_degree[SR[i]] < R_end);
			total_edges += degree[SR[i]];
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end;i ++) {
			ui u, min_degree = n;
			for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]&&t_degree[SR[j]] < min_degree) {
				u = SR[j];
				min_degree = t_degree[u];
			}
			assert(min_degree < n);
			if(t_degree[u] > max_core) max_core = t_degree[u];
			core[u] = max_core;

			long long t_n = R_end - i;

			if(t_n*(t_n-1)/2 <= total_edges/2 + K&&R_end - i > best_solution_size) {
				best_solution_size = R_end - i;
				ui cnt = 0;
				for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]) best_solution[cnt ++] = SR[j];
				assert(cnt == best_solution_size);
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			vis[u] = 1;
			if(t_degree[u] == 0) continue;
			char *t_matrix = matrix + u*n;
			for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]&&t_matrix[SR[j]]) {
				-- t_degree[SR[j]];
				total_edges -= 2;
			}
		}
		for(ui i = 0;i < R_end;i ++) vis[SR[i]] = 0;

		assert(Qv.empty());
		for(ui i = 0;i < R_end;i ++) if(core[SR[i]] + K < best_solution_size) {
			assert(removed_level[SR[i]] > level);
			removed_level[SR[i]] = level;
			Qv.push(SR[i]);
		}
    }

    ui degeneracy_ordering_and_coloring(ui S_end, ui R_end, ui *color) {
		ui *rid = buf;
		ui *peel_sequence = buf1;
		ui *t_degree = color;
		for(ui i = S_end;i < R_end;i ++) {
			t_degree[SR[i]] = degree[SR[i]] - degree_in_S[SR[i]];
#ifndef NDEBUG
			ui d = 0;
			char *t_matrix = matrix + SR[i]*n;
			assert(!vis[SR[i]]);
			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ d;
			assert(d == t_degree[SR[i]]);
#endif
		}

		for(ui i = S_end;i < R_end;i ++) {
			ui u, min_degree = n;
			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_degree[SR[j]] < min_degree) {
				u = SR[j];
				min_degree = t_degree[u];
			}
			assert(min_degree < n);
			peel_sequence[i-S_end] = u;
			rid[u] = i-S_end;

			vis[u] = 1;
			if(t_degree[u] == 0) continue;
			char *t_matrix = matrix + u*n;
			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_matrix[SR[j]]) -- t_degree[SR[j]];
		}
		for(ui i = S_end;i < R_end;i ++) vis[SR[i]] = 0;

    	ui max_color = 0;
    	for(ui i = R_end-S_end;i > 0;i --) {
    		ui u = peel_sequence[i-1];
			char *t_matrix = matrix + u*n;
			for(ui j = S_end;j < R_end;j ++) if(rid[SR[j]] >= i&&t_matrix[SR[j]]) vis[color[SR[j]]] = 1;
    		for(ui j = 0;;j ++) if(!vis[j]) {
    			color[u] = j;
    			if(j > max_color) max_color = j;
    			break;
    		}
			for(ui j = S_end;j < R_end;j ++) if(rid[SR[j]] >= i&&t_matrix[SR[j]]) vis[color[SR[j]]] = 0;
    	}

    	return max_color + 1;
    }

    bool is_kDefectiveClique(ui R_end) {
    	ui total_edges = 0;
    	for(ui i = 0;i < R_end;i ++) total_edges += degree[SR[i]];
    	long long all_edges = R_end;
    	return all_edges*(R_end-1)/2 <= total_edges/2 + K;
    }

    void store_a_kDefectiveClique(ui S_end) {
    	assert(S_end > best_solution_size);
		best_solution_size = S_end;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
		printf("Find a kDefectiveClique of size: %u\n", best_solution_size);
#endif
    }

    void BB_search(ui S_end, ui R_end, ui level, bool must_include_0) {
#ifndef NDEBUG
    	assert(compute_missing_edges_in_S(S_end) <= K);
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	for(ui i = 0;i < R_end;i ++) assert(!vis[SR[i]]);
#endif

        if(best_solution_size >= UB||R_end <= best_solution_size) return ;

        fflush(stdout);

    	if(S_end > best_solution_size) store_a_kDefectiveClique(S_end);
        if(R_end > best_solution_size&&is_kDefectiveClique(R_end)) store_a_kDefectiveClique(R_end);
        if(R_end <= best_solution_size+1) return ;

        // printf("level %u\n", level);

#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			char *t_matrix = matrix + u*n;
			for(ui j = 0;j < S_end;j ++) if(t_matrix[SR[j]]) ++ cnt;
			assert(degree_in_S[u] == cnt);
			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ cnt;
			assert(degree[u] == cnt);
		}
#endif

        ui old_R_end = R_end;

        if(S_end == 0) {
        	assert(Qv.empty());
        	compute_a_heuristic_solution_and_prune(R_end, level);
        	if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
        		return ;
        	}
        }

        assert(Qv.empty());
        bool terminate = collect_removable_vertices_based_on_degree_in_S(S_end, R_end, level);
        if(terminate||remove_vertices_and_prune(S_end, R_end, level)) {
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }

#ifdef _RR_OPT_
        assert(Qv.empty());
        terminate = collect_removable_vertices_based_on_degree_in_S_opt(S_end, R_end, level);
        if(terminate||remove_vertices_and_prune(S_end, R_end, level)) {
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }
#else
        if(S_end > 0) {
        	assert(Qv.empty());
        	collect_removable_vertices_based_on_vertex_pair(S_end, R_end, level); // remove vertices from R
			if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
				return ;
			}
        }
#endif

		if(R_end > best_solution_size&&is_kDefectiveClique(R_end)) store_a_kDefectiveClique(R_end);
		if(R_end <= best_solution_size+1) {
        	restore_R(S_end, R_end, old_R_end, level);
			return ;
		}

#ifndef _RR_OPT_
		if(degree_based_prune(S_end, R_end, level)||coloring_based_prune(S_end, R_end, level)) {
			restore_R(S_end, R_end, old_R_end, level);
			return ;
		}
#else
		// if(coloring_based_prune(S_end, R_end, level)) {
		// 	restore_R(S_end, R_end, old_R_end, level);
		// 	return ;
		// }
#endif

		//printf("here2\n");

        ui u = n; // u is the branching vertex
        bool must_include = false;
        if(must_include_0) {
        	assert(S_end == 0);
        	if(SR_rid[0] >= R_end) {
        		restore_R(S_end, R_end, old_R_end, level);
        		return ;
        	}
        	u = 0;
        	must_include = true;
        }
        if(u == n) {
			for(ui i = S_end;i < R_end;i ++) if(degree[SR[i]] + 2 >= R_end) {
				must_include = true;
				u = SR[i];
				break;
			}
        }
		if(u == n) {
			u = choose_branch_vertex_based_on_degree(S_end, R_end, level);
			// u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end, level);
		}
		assert(u != n&&SR_rid[u] >= S_end);
		assert(SR_rid[u] < R_end);
		assert(SR[SR_rid[u]] == u);
        assert(degree[u] + K >= best_solution_size);

#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			char *t_matrix = matrix + u*n;
			for(ui j = 0;j < S_end;j ++) if(t_matrix[SR[j]]) ++ cnt;
			assert(degree_in_S[u] == cnt);
			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ cnt;
			assert(degree[u] == cnt);
		}
#endif

        //printf("here3\n");

        // the first branch includes u into S
        ui pre_best_solution_size = best_solution_size, t_old_R_end = R_end;
        if(!move_u_from_R_to_S(u, S_end, R_end, level)) {
        	if(!remove_vertices_and_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1, false);
        }
        restore_R(S_end, R_end, t_old_R_end, level);

        if(must_include) {
        	move_u_from_S_to_R(S_end, R_end, level);
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }
        //printf("back to level: %u\n", level);

        // the second branch excludes u from S
        assert(Qv.empty());
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
        if(!pruned&&best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_and_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1, false);
        }
        restore_R(S_end, R_end, old_R_end, level);
    }

    bool collect_removable_vertices_based_on_degree_in_S(ui S_end, ui R_end, ui level) {
    	if(S_end >= R_end) return false;
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++) ++ cnt[S_end - degree_in_S[SR[i]]];
    	assert(R_end > best_solution_size);

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui remaining_vertices_n = best_solution_size - S_end, idx = S_end;
    	assert(missing_edges <= K&&S_end <= best_solution_size);
    	assert(Qv.empty());
    	for(ui i = 0;i <= S_end;i ++) {
    		if(cnt[i] < remaining_vertices_n) {
    			remaining_vertices_n -= cnt[i];
    			missing_edges += i*cnt[i];
    		}
    		else {
    			missing_edges += i*remaining_vertices_n;
    			idx = i;

    			ui next_value = i;
    			if(cnt[i] == remaining_vertices_n) {
    				for(i ++;i <= S_end;i ++) if(cnt[i]) {
    					next_value = i;
    					break;
    				}
    			}
    			if(missing_edges + next_value > K) return true;

    			break;
    		}
    	}
    	for(ui i = S_end;i < R_end;i ++) if(S_end-degree_in_S[SR[i]] > idx&&missing_edges+S_end-degree_in_S[SR[i]] > K) {
    		removed_level[SR[i]] = level;
    		Qv.push(SR[i]);
    	}

    	return false;
    }

    bool collect_removable_vertices_based_on_degree_in_S_opt(ui S_end, ui R_end, ui level) {
    	if(S_end >= R_end) return false;
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++) ++ cnt[S_end - degree_in_S[SR[i]]];
    	assert(R_end > best_solution_size);

    	ui *candidates = buf1;
    	for(ui i = 1;i <= S_end;i ++) cnt[i] += cnt[i-1];
    	for(ui i = S_end;i < R_end;i ++) candidates[-- cnt[S_end-degree_in_S[SR[i]]]] = SR[i];
#ifndef NDEBUG
    	for(ui i = 1;i < R_end-S_end;i ++) assert(degree_in_S[candidates[i]] <= degree_in_S[candidates[i-1]]);
#endif

    	ui *ids = buf;
    	ui ids_n = 0;

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui remaining_vertices_n = best_solution_size - S_end;
    	assert(missing_edges <= K&&S_end <= best_solution_size);
    	assert(Qv.empty());
    	if(remaining_vertices_n >= R_end-S_end) return true;
    	ui t_missing_edges = missing_edges;
    	for(ui i = 0;i <= remaining_vertices_n;i ++) t_missing_edges += S_end - degree_in_S[candidates[i]];
    	if(t_missing_edges > K) return true;

    	ui *suffix_sum = buf2;
    	suffix_sum[R_end-S_end] = 0;
    	for(ui i = 0;i < R_end-S_end;i ++) {
    		ui j = R_end-S_end-1-i;
    		suffix_sum[j] = suffix_sum[j+1] + S_end - degree_in_S[candidates[j]];
    	}

    	for(ui i = 0;i < R_end-S_end&&remaining_vertices_n;i ++) {
    		ui u = candidates[i];
    		if(missing_edges + S_end - degree_in_S[u] > K||ids_n + R_end-S_end-i <= remaining_vertices_n) {
    			removed_level[u] = level;
    			Qv.push(u);
    			continue;
    		}

    		t_missing_edges = missing_edges;
    		ui end_idx = i + 1;
    		if(ids_n < remaining_vertices_n) end_idx = i + 1 + remaining_vertices_n - ids_n;
    		t_missing_edges += suffix_sum[i] - suffix_sum[end_idx];
    		if(t_missing_edges > K) {
    			removed_level[u] = level;
    			Qv.push(u);
    			continue;
    		}

    		if(t_missing_edges+R_end-S_end-1+degree_in_S[u]-degree[u] > K) {
    			char *t_matrix = matrix + u*n;
    			if(ids_n >= remaining_vertices_n) {
    				ui cut_off = degree_in_S[ids[remaining_vertices_n-1]];
    				ui same_degree_nonneighbors = 0;
    				for(ui j = 0;j < remaining_vertices_n&&t_missing_edges+remaining_vertices_n-j > K;j ++) if(!t_matrix[ids[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[ids[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = remaining_vertices_n;j < ids_n&&t_missing_edges > K;j ++) {
    					if(degree_in_S[ids[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[ids[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    				for(ui j = i+1;j < R_end-S_end&&t_missing_edges > K;j ++) {
    					if(degree_in_S[candidates[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[candidates[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    			}
    			else {
    				ui cut_off = degree_in_S[candidates[i+remaining_vertices_n-ids_n]];
    				ui same_degree_nonneighbors = 0;
    				for(ui j = 0;j < ids_n&&t_missing_edges+remaining_vertices_n-j > K;j ++) if(!t_matrix[ids[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[ids[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = i+1;j <= i+remaining_vertices_n-ids_n&&t_missing_edges+remaining_vertices_n-ids_n-j+i+1>K;j ++) if(!t_matrix[candidates[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[candidates[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = i+1+remaining_vertices_n-ids_n;j < R_end-S_end&&t_missing_edges > K;j ++) {
    					if(degree_in_S[candidates[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[candidates[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    			}

    			if(t_missing_edges > K) {
    				removed_level[u] = level;
    				Qv.push(u);
    				continue;
    			}
    		}

    		if(ids_n < remaining_vertices_n) missing_edges += S_end - degree_in_S[u];
    		ids[ids_n++] = u;
    	}

    	return false;
    }

    void collect_removable_vertices_based_on_vertex_pair(ui S_end, ui R_end, ui level) {
       	assert(S_end >= 1&&Qv.empty());
       	ui u = SR[S_end-1], neighbors_n = 0, non_neighbors_n = 0;
       	char *t_matrix = matrix + u*n;
       	for(ui i = S_end;i < R_end;i ++) {
       		if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
       		else nonneighbors[non_neighbors_n ++] = SR[i];
       	}

       	ui missing_edges_in_S = compute_missing_edges_in_S(S_end);
       	assert(missing_edges_in_S <= K);

       	for(ui i = 0;i < neighbors_n;i ++) {
       		ui v = neighbors[i];
       		ui common_neighbors = 0, common_non_neighbors = 0;
       		t_matrix = matrix + v*n;
       		for(ui j = 0;j < neighbors_n;j ++) if(t_matrix[neighbors[j]]) ++ common_neighbors;
       		for(ui j = 0;j < non_neighbors_n;j ++) if(!t_matrix[nonneighbors[j]]) ++ common_non_neighbors;
       		ui exclusive_non_neighbors = neighbors_n + non_neighbors_n - 1 - common_neighbors - common_non_neighbors;

       		ui local_UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			local_UB += mmin(tmp, common_non_neighbors);
       		}
       		if(local_UB <= best_solution_size) {
       			removed_level[v] = level;
       			Qv.push(v);
       			neighbors[i] = neighbors[-- neighbors_n];
       			-- i;
       		}
       	}

       	++ missing_edges_in_S;
       	assert(missing_edges_in_S <= K||non_neighbors_n == 0);
       	for(ui i = 0;i < non_neighbors_n;i ++) {
   			ui v = nonneighbors[i];
   			ui common_neighbors = 0, common_non_neighbors = 0;
   			t_matrix = matrix + v*n;
       		for(ui j = 0;j < neighbors_n;j ++) if(t_matrix[neighbors[j]]) ++ common_neighbors;
       		for(ui j = 0;j < non_neighbors_n;j ++) if(j != i&&!t_matrix[nonneighbors[j]]) ++ common_non_neighbors;
       		ui exclusive_non_neighbors = neighbors_n + non_neighbors_n - 1 - common_neighbors - common_non_neighbors;

       		ui local_UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			local_UB += mmin(tmp, common_non_neighbors);
       		}
   			if(local_UB <= best_solution_size) {
   				removed_level[v] = level;
   				Qv.push(v);
   				nonneighbors[i] = nonneighbors[-- non_neighbors_n];
   				-- i;
   			}
   		}
    }

    bool collect_removable_vertices(ui S_end, ui R_end, ui level) {
		for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K < best_solution_size) return true;

		for(ui i = S_end;i < R_end;i ++) if(removed_level[SR[i]] > level){
			ui v = SR[i];
			if(degree[v] + K < best_solution_size) {
				removed_level[v] = level;
				Qv.push(v);
			}
		}

		return false;
	}

    bool degree_based_prune(ui S_end, ui R_end, ui level) {
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++) ++ cnt[S_end - degree_in_S[SR[i]]];

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui t_UB = S_end;
    	for(ui i = 0;i <= S_end;i ++) for(ui j = 0;j < cnt[i];j ++) {
    		if(missing_edges + i > K) break;

    		missing_edges += i;
    		++ t_UB;
    	}

    	return t_UB <= best_solution_size;
    }

    bool coloring_based_prune(ui S_end, ui R_end, ui level) {
#ifndef NDEBUG
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif
    	ui *color = neighbors;
    	ui color_n = degeneracy_ordering_and_coloring(S_end, R_end, color);
    	// printf("color_n: %u\n", color_n);

    	ui missing_edges_n = compute_missing_edges_in_S(S_end);
    	assert(missing_edges_n <= K);
    	if(S_end + color_n + K - missing_edges_n <= best_solution_size) return true;

#ifndef NDEBUG
    	// for(ui i = 0;i < n;i ++) assert(!vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] = 1;
    	for(ui i = 0;i < color_n;i ++) assert(vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] = 0;
    	//for(ui i = S_end;i < R_end;i ++) if(color[SR[i]] == 16) printf("color of %u is 16\n", SR[i]);
#endif

    	ui *head = nonneighbors;
    	ui *next = buf;
    	for(ui i = 0;i <= S_end;i ++) head[i] = n;
    	for(ui i = S_end;i < R_end;i ++) {
    		assert(degree_in_S[SR[i]] <= S_end);
    		ui u = SR[i], non_neighbors_n = S_end - degree_in_S[SR[i]];
    		//if(color[u] == 16) printf("color[%u] is 16\n", u);
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
    	}

    	ui *color_head = buf1;
    	ui *color_next = buf2;
    	for(ui i = 0;i < color_n;i ++) color_head[i] = n;
    	for(ui i = 0;i <= S_end;i ++) for(ui u = head[S_end-i];u != n;u = next[u]) {
    		ui c = color[u];
    		assert(c < color_n);
    		color_next[u] = color_head[c];
    		color_head[c] = u;
    		//printf("color_head[%u] = %u\n", c, u);
    	}

    	for(ui i = 0;i <= K;i ++) head[i] = n;
    	for(ui i = 0;i < color_n;i ++) {
    		ui u = color_head[i];
    		//printf("color_head[%u]: %u\n", i, u);
    		assert(u != n&&color[u] == i);
    		color_head[i] = color_next[u];

    		ui non_neighbors_n = S_end - degree_in_S[u];
    		assert(non_neighbors_n <= K);
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
    	}

    	ui min_key = 0, t_UB = S_end;
    	while(true) {
    		while(min_key <= K&&head[min_key] == n) ++ min_key;
    		if(min_key > K||missing_edges_n + min_key > K) break;

    		missing_edges_n += min_key;
    		++ t_UB;
    		ui u = head[min_key];
    		assert(u != n);
    		head[min_key] = next[u];
    		ui new_u = color_head[color[u]];
    		if(new_u == n) continue;

    		color_head[color[u]] = color_next[new_u];
    		ui non_neighbors_n = S_end - degree_in_S[u];
    		ui new_key = S_end - degree_in_S[new_u] + 1 + min_key - non_neighbors_n;
    		if(new_key > K) continue;
    		next[new_u] = head[new_key];
    		head[new_key] = new_u;
    	}

    	return t_UB <= best_solution_size;
    }

    bool move_u_from_R_to_S(ui u, ui &S_end, ui R_end, ui level) {
    	assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
    	swap_pos(S_end, SR_rid[u]);
    	++ S_end;

    	char *t_matrix = matrix + u*n;
    	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
			++ degree_in_S[SR[i]];
			++ edges_in_s;
		}

        ui missing_edges = compute_missing_edges_in_S(S_end);
        assert(missing_edges <= K);

        assert(Qv.empty());
        for(ui i = S_end;i < R_end;i ++) if(S_end - degree_in_S[SR[i]] + missing_edges > K) {
        	removed_level[SR[i]] = level;
        	Qv.push(SR[i]);
        }

        return false;
    }

    void move_u_from_S_to_R(ui &S_end, ui R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[-- S_end];
		char *t_matrix = matrix + u*n;
		for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
			-- degree_in_S[SR[i]];
			-- edges_in_s;
		}
	}

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		removed_level[u] = level;

		bool ret = false;
		char *t_matrix = matrix + u*n;

		for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
			ui v = SR[i];
			-- degree_in_S[v];
			--edges_in_s;
			-- degree[v];
			m--;
			if(degree[v] + K < best_solution_size) {
				if(SR_rid[v] < S_end) ret = true;
				else {
					assert(removed_level[v] > level);
					removed_level[v] = level;
					Qv.push(v);
				}
			}
		}
		return ret;
	}

    bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			// printf("remove %u with degree %u\n", u, degree[u]);

			bool terminate = false;
			char *t_matrix = matrix + u*n;
			for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
				ui w = SR[i];
				assert(degree[w] > 0);
				-- degree[w];
				m--;
				if(degree[w] + K < best_solution_size) {
					if(SR_rid[w] < S_end) terminate = true; // UB
					else if(removed_level[w] > level) { // RR
						removed_level[w] = level;
						Qv.push(w);
					}
				}
			}
			if(terminate) return true;
		}

        return false;
    }

    void restore_R(ui S_end, ui &R_end, ui old_R_end, ui level) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(removed_level[u] == level&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
            removed_level[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(removed_level[u] == level&&SR_rid[u] == R_end);
            removed_level[u] = n;

            char *t_matrix = matrix + u*n;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]){
				++ degree[SR[i]];
				m++;
			}

            ++ R_end;
            // printf("restore %u with degree %u\n", u, degree[u]);
        }
    }
    
    ui compute_missing_edges_in_S(ui S_end) {
    	ui res = 0;
    	for(ui i = 0;i < S_end;i ++) {
    		assert(degree_in_S[SR[i]] < S_end);
    		res += S_end - 1 - degree_in_S[SR[i]];
    	}
    	return res/2;
    }

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex_based_on_degree(ui S_end, ui R_end, ui level) {
    	//printf("start choosing branching vertices\n");
    	ui u = SR[S_end];
		for(ui i = S_end+1;i < R_end;i ++) {
			ui v = SR[i];
#ifndef NEW_BRANCH
			if(degree_in_S[u] == S_end) {
				if(degree_in_S[v] != S_end||degree[v] > degree[u]) u = v;
			}
			else if(degree_in_S[v] != S_end&&(degree[v] > degree[u]||(degree[v] == degree[u]&&degree_in_S[v] < degree_in_S[u]))) u = v;	
#else
			//good at sparse graph
			ui vn=R_end-S_end, edge_num=m-edges_in_s;
	#ifdef SPARSE
				if(degree_in_S[v]>degree_in_S[u]) u=v;
				else if(degree_in_S[v]==degree_in_S[u]&&degree[v]>degree[u]) u=v; //good at ia-email-EU 0.9 and inf-road-use
	#endif
	#ifdef MID
				if(degree_in_S[v]>degree_in_S[u]) u=v; 
				else if(degree_in_S[v]==degree_in_S[u]&&degree[v]<degree[u]) u=v;
	#endif
	#ifdef DENSE
				if(degree_in_S[v]<degree_in_S[u]) u=v; //mid
				else if(degree_in_S[v]==degree_in_S[u]&&degree[v]<degree[u]) u=v;
	#endif
	#ifdef MIX
				if(degree[v]-degree_in_S[v]<degree[u]-degree_in_S[u]) u=v;
	#endif
#endif
		}

#ifndef NDEBUG
		ui missing_edges_n = compute_missing_edges_in_S(S_end);
		assert(u != n&&missing_edges_n + S_end - degree_in_S[u] <= K);
		assert(degree[u] + 2 < R_end);
#endif

		return u;
    }

    ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end, ui level) {
    	//printf("start choosing branching vertices\n");
    	ui min_degree_in_S = S_end+1, u = n;
		for(ui i = S_end;i < R_end;i ++) {
			ui v = SR[i];
			if(degree_in_S[v] < min_degree_in_S) {
				min_degree_in_S = degree_in_S[v];
				u = v;
			}
			else if(degree_in_S[v] == min_degree_in_S&&degree[v] > degree[u]) u = v;
		}

#ifndef NDEBUG
		ui missing_edges_n = compute_missing_edges_in_S(S_end);
		assert(u != n&&missing_edges_n + S_end - degree_in_S[u] <= K);
		assert(degree[u] + 2 < R_end);
#endif

		return u;
    }
};

#endif
