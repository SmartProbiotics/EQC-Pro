#include "Graph.h"
#include "SuperHeu.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"
using namespace std;
using namespace popl;
int my_ceil(double a)
{
    int l = floor(a), r = ceil(a);
    if (abs(a - l) < 1e-6)
        return l;
    return r;
}
int my_floor(double a)
{
    int l = floor(a), r = ceil(a);
    if (abs(a - r) < 1e-6)
        return r;
    return l;
}

unordered_map<int, int> st;
Graph static_g;
double gama;

string fileName;
int searchTime = 0;

int f(int k, int lb, int ub, vector<ui> &result)
{
    cout << "<<==============================Exact Search==========================>>" << endl;
    cout << "[lb]: " << lb << " [k]: " << k << endl;
    if (st.count(k))
        return st[k];
    Graph graph(static_g);
    graph.K = k;
    int s = clock();
    int n = graph.kDefectiveClique_exact(lb, graph.n);
    searchTime += clock() - s;
    if (n > lb)
    {
        st[k] = n;
        result = graph.kDefectiveClique;
    }
    return n;
}

bool check_sol(Graph *graph, double gamma, vector<ui> &result)
{
    ui n = graph->n;
    bool *st = new bool[n];
    memset(st, 0, sizeof(bool) * n);
    for (auto &t : result)
        st[t] = true;
    ui m = 0;
    for (auto &u : result)
        for (ept i = graph->pstart[u]; i != graph->pstart[u + 1]; i++)
            if (st[graph->edges[i]])
                m++;
    n = result.size();
    delete[] st;
    cout << "Edge Cnt :" << m << " ";
    cout << "Result Density: " << double(m) / n / (n - 1) << endl;
    return (double(m) / n >= (n - 1) * gamma);
}

void degen(Graph *graph, long double gamma, int &LB, int UB, vector<ui> &result, ui *&peel_sequence)
{
    graph->kDefectiveClique.clear();
    int n = graph->n;
    ept *pstart = graph->pstart;
    ui *edges = graph->edges;
    peel_sequence = new ui[n];
    ui *core = new ui[n];
    ui *degree = new ui[n];
    char *vis = new char[n];
    ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
    for (ui i = 0; i < n; i++)
        degree[i] = pstart[i + 1] - pstart[i];
    ept total_edges = 0;
    memset(vis, 0, sizeof(char) * n);
    for (ui i = 0; i < n; i++)
    {
        peel_sequence[i] = i;
        total_edges += degree[i];
    }
    printf("[Total Edges:%d]\n", total_edges);
    heap->init(n, n - 1, peel_sequence, degree);
    ui max_core = 0;
    ui idx = n;
    for (ui i = 0; i < n; i++)
    {
        ui u, key;
        heap->pop_min(u, key);
        if (key > max_core)
            max_core = key;
        core[u] = max_core;
        peel_sequence[i] = u;
        long double a = 1 - gamma;
        if (idx == n && (double)total_edges / (n - i) >= (n - i - 1) * gamma)
            idx = i;
        vis[u] = 1;
        for (ept j = pstart[u]; j < pstart[u + 1]; j++)
            if (vis[edges[j]] == 0)
            {
                total_edges -= 2;
                heap->decrement(edges[j], 1);
            }
    }
    LB = n - idx;
    if (result.size() < LB)
    {
        result.clear();
        for (int i = idx; i < n; i++)
            result.push_back(peel_sequence[i]);
    }
    delete heap;
    delete[] core;
    delete[] vis;
    delete[] degree;
}

#include <vector>
#include <cstring>
#include <algorithm>

void greedyColoring(Graph *graph, const ui *peel_sequence, int *&color_cnt, int &max_color)
{
    ept *pstart = graph->pstart;
    ui *edges = graph->edges;
    int n = graph->n;
    int *color = new int[n];
    std::memset(color, 0, n * sizeof(int));
    max_color = 1;

    std::vector<int> forbiddenTime(n + 2, 0);
    int curTime = 1;

    for (int idx = 0; idx < n; ++idx)
    {
        int u = peel_sequence[n - idx - 1];
        if (color[u] != 0)
            continue;
        for (ept ei = pstart[u]; ei < pstart[u + 1]; ++ei)
        {
            int v = edges[ei];
            int c = color[v];
            if (c)
                forbiddenTime[c] = curTime;
        }
        int c;
        for (c = 1; c <= max_color; ++c)
        {
            if (forbiddenTime[c] != curTime)
                break;
        }
        color[u] = c;
        if (c == max_color)
            ++max_color;
        ++curTime;
    }

    // 分配 color_cnt 并统计每个颜色数量
    color_cnt = new int[max_color];
    std::memset(color_cnt, 0, max_color * sizeof(int));
    for (int i = 0; i < n; ++i)
    {
        ++color_cnt[color[i]];
    }
    delete[] color;
}

int computeUB(int k, int *color_cnt, int max_color)
{
    int max_cnt = 0;
    for (int c = 1; c < max_color; ++c)
    {
        max_cnt = std::max(max_cnt, color_cnt[c]);
    }

    int l = 0, r = max_cnt;
    while (l < r)
    {
        int mid = (l + r + 1) >> 1;
        long long sum_defective = 0;
        for (int c = 1; c < max_color; ++c)
        {
            int cnt = std::min(mid, color_cnt[c]);
            sum_defective += 1LL * cnt * (cnt - 1) / 2;
        }
        if (sum_defective <= k)
        {
            l = mid;
        }
        else
        {
            r = mid - 1;
        }
    }

    int UB = 0;
    for (int c = 1; c < max_color; ++c)
    {
        int cnt = std::min(color_cnt[c], l);
        k -= cnt * (cnt - 1) / 2;
        UB += cnt;
    }
    if (l > 0)
        UB += k / (l + 1);
    return UB;
}

void ego_degen(Graph *graph, double gamma, vector<ui> &result, ui *peel_sequence)
{
    ui n = graph->n, m = graph->m;
    ept *pstart = graph->pstart;
    ui *edges = graph->edges;
    ui *peel_id = new ui[n], *deg = new ui[n];
    ui *subID = new ui[n];
    ui *vset = new ui[n], vnum = 0;
    ui *ids = new ui[n];
    memset(deg, 0, sizeof(ui) * n);
    memset(subID, -1, sizeof(ui) * n);

    for (int i = 0; i < n; i++)
        peel_id[peel_sequence[i]] = i;
    for (int i = 0; i < n - result.size(); i++)
    {
        vnum = 0;
        ui u = peel_sequence[i];
        for (ept j = pstart[u]; j != pstart[u + 1]; j++)
        {
            ui v = edges[j];
#ifndef NHDC
            if (peel_id[v] > peel_id[u])
#endif
            {
                subID[v] = vnum;
                vset[vnum++] = v;
            }
        }
        if (vnum <= result.size())
        {
            for (int j = 0; j != vnum; j++)
                subID[vset[j]] = -1;
            continue;
        }
        vector<vector<ui>> neig;
        neig.resize(vnum);
        ui total_edges = 0;
        for (int j = 0; j < vnum; j++)
        {
            ui v = vset[j];
            for (ept k = pstart[v]; k != pstart[v + 1]; k++)
                if (subID[edges[k]] != -1)
                {
                    ui w = edges[k];
                    neig[subID[v]].push_back(subID[w]);
                }
            deg[j] = neig[j].size();
            total_edges += deg[j];
        }
        // cout<<"Total Edges:"<<total_edges<<" Vnum:"<<vnum<<endl;
        // do degen here

        for (int j = 0; j < vnum; j++)
            ids[j] = j;
        ListLinearHeap *heap = new ListLinearHeap(vnum, vnum - 1);
        heap->init(vnum, vnum - 1, ids, deg);
        ui res_idx = vnum;
        for (ui j = 0; j < vnum; j++)
        {
            if (double(total_edges + 2 * (vnum - j)) / (vnum + 1 - j) >= (vnum + 1 - j - 1) * gamma)
            {
                res_idx = j;
                break;
            }
            ui mid, key;
            heap->pop_min(mid, key);
            deg[mid] = 0;
            for (auto v : neig[mid])
            {
                if (deg[v] == 0)
                    continue;
                total_edges -= 2;
                heap->decrement(v, 1);
            }
        }
        if (result.size() < vnum - res_idx + 1)
        {
            result.clear();
            int num = vnum - res_idx;
            while (num--)
            {
                ui mid, k;
                heap->pop_min(mid, k);
                assert(subID[vset[mid]] != -1);
                result.push_back(vset[mid]);
            }
            result.push_back(u);
        }
        delete heap;
        for (int j = 0; j != vnum; j++)
        {
            deg[j] = 0;
            subID[vset[j]] = -1;
        }
    }
    delete[] ids;
    delete[] vset;
    delete[] subID;
    delete[] peel_id;
    delete[] deg;
}

void print_result(vector<ui> result)
{
    cout << "Result Set:" << endl;
    sort(result.begin(), result.end());
    copy(result.begin(), result.end(), ostream_iterator<int>(std::cout, " "));
    cout << endl;
}

bool heuristic_find(int k, vector<ui> &res, ui *to_static, ui *&out_mapping)
{
#ifdef NO_EXTEND
    return false;
#endif
    Graph tmp;
    tmp.K = k;
    int n = static_g.n;
    vector<ui> tmp_res = res;
    bool *st = new bool[n];
    ui *id = new ui[n], *neig = new ui[n], ne_idx = 0;
    memset(st, 0, sizeof(bool) * n);
    memset(id, 0, sizeof(ui) * n);
    for (auto &t : tmp_res)
        t = to_static[t];
    auto pstart = static_g.pstart, edges = static_g.edges;
    int new_n = 0, new_m = 0;

    for (auto &t : tmp_res)
        st[t] = true;
    for (auto &t : tmp_res)
        neig[ne_idx++] = t;
    for (auto &t : tmp_res)
    {
        for (ept i = pstart[t]; i != pstart[t + 1]; i++)
        {
            if (!st[edges[i]])
            {
                neig[ne_idx++] = edges[i];
                st[edges[i]] = true;
            }
        }
    }

    for (int i = 0; i < ne_idx; i++)
    {
        ui t = neig[i];
        for (ept j = pstart[t]; j != pstart[t + 1]; j++)
            if (st[edges[j]])
                new_m++;
    }

    for (int i = 0; i < ne_idx; i++)
        id[neig[i]] = new_n++;

    ui *npstart = new ui[new_n + 1];
    ept *nedges = new ept[new_m];

    int idx = 0;
    for (int i = 0; i < new_n; i++)
    {
        ui v = neig[i];
        npstart[i] = idx;
        for (ept j = pstart[v]; j != pstart[v + 1]; j++)
            if (st[edges[j]])
                nedges[idx++] = id[edges[j]];
        if (i == new_n - 1)
            npstart[new_n] = idx;
    }
    cout << "Sub graph V:" << new_n << " " << " M:" << new_m << endl;
    assert(idx == new_m);
    tmp.pstart = npstart;
    tmp.edges = nedges;
    tmp.n = new_n;
    tmp.m = new_m;

    // tmp.kDefectiveClique_exact(check-1,new_n);
    // tmp_res=tmp.kDefectiveClique;
    for (auto &t : tmp_res)
        t = id[t];
    print_result(tmp_res);
    HeuExtractor he(tmp, tmp_res, gama);
    he.extract();
    tmp_res = he.best_res;
    cout << "Solution Density: " << he.sol_density << endl;
    if (he.sol_density >= gama)
    {
        cout << "<<Successful Heuristic!>>" << endl;
        // print_result(tmp_res);
        cout << "Find Result:" << tmp_res.size() << endl;
        print_result(tmp_res);
        for (auto &t : tmp_res)
            t = out_mapping[neig[t]];
        res = tmp_res;
    }
    else
    {
        cout << "<<failed Heuristic!>>" << endl;
        cout << "Find Result:" << tmp_res.size() << endl;
    }
    delete[] st;
    delete[] id;
    delete[] neig;
    return he.sol_density >= gama;
}

void remove_vertice(Graph *graph, vector<ui> &vset, ui *&out_mapping)
{
    bool *st = new bool[graph->n];
    memset(st, 0, sizeof(bool) * graph->n);
    ui *deg = new ui[graph->n];
    memset(deg, 0, sizeof(ui) * graph->n);
    ept *pstart = new ept[graph->n - vset.size() + 1];
    int m = 0;
    for (auto &t : vset)
        st[t] = true;
    for (ui i = 0; i < graph->n; i++)
    {
        if (!st[i])
            for (ept j = graph->pstart[i]; j != graph->pstart[i + 1]; j++)
                if (!st[graph->edges[j]])
                    deg[i]++;
        m += deg[i];
    }

    int n = graph->n - vset.size();
    ui *edges = new ui[m + 1];
    pstart[0] = 0;
    int idx = 0;
    ui *id = new ui[graph->n];
    ui *rid = new ui[n];
    for (int i = 0; i < graph->n; i++)
        if (!st[i])
        {
            id[i] = idx;
            rid[idx++] = i;
        }
    out_mapping = rid;
    idx = 0;
    for (int i = 1; i <= n; i++)
        pstart[i] = pstart[i - 1] + deg[rid[i - 1]];
    for (ui i = 0; i < n; i++)
    {
        for (ept j = graph->pstart[rid[i]]; j != graph->pstart[rid[i] + 1]; j++)
            if (!st[graph->edges[j]])
            {
                edges[idx++] = id[graph->edges[j]];
            }
    }
    delete[] graph->pstart;
    delete[] graph->edges;
    graph->pstart = pstart;
    graph->edges = edges;
    graph->n = n;
    graph->m = m;
    delete[] deg;
    delete[] st;
    delete[] id;
}

void extend(double alpha, int k_cur, vector<ui> &result, ui *&to_static_mapping, ui *&out_mapping, int &lb, bool do_first_heuris_n = false)
{
    if (alpha < 1e-7)
        return;
    cout << "Heu Division Val: " << my_floor(1 / alpha + 0.5) << endl;
    if (do_first_heuris_n)
    {
        puts("===================Heuris N==============================");
        cout << "[N]" << lb << " [K]" << k_cur << endl;
        if (heuristic_find(k_cur, result, to_static_mapping, out_mapping))
            lb = result.size();
    }

    while (1)
    {
        bool n_increase_slow = (lb <= my_floor(1 / alpha + 0.5));
#ifdef ONLY_N
        n_increase_slow = false;
#endif
        if (!n_increase_slow)
        {
            int k_next = my_floor((lb + 1) * lb / 2 * alpha);
            if (k_next == k_cur)
                break;
            lb++;
            k_cur = k_next;
            puts("===================Heuris N==============================");
            cout << "[N]" << lb << " [K]" << k_cur << endl;
            heuristic_find(k_cur, result, to_static_mapping, out_mapping);
            lb = result.size();
        }
        else
        {
            k_cur++;
            k_cur = max(k_cur, my_floor((lb + 1) * lb / 2 * alpha)); // sufficiently use lb. use lb+1 to create a k value.
            int n_limit = my_ceil((1 + sqrt(1 + 8 * k_cur / alpha)) / 2);
            puts("===================Heuris K==============================");
            cout << "[N]" << n_limit << " [K]" << k_cur << endl;
            if (heuristic_find(k_cur, result, to_static_mapping, out_mapping))
                lb = result.size();
            else
                break;
        }
    }
}

int main(int argc, char *argv[])
{
    bool output = false;
    bool binary_input = true;
    long long s1 = clock();
    vector<ui> result;
    gama = stod(string(argv[2]));
    fileName = string(argv[1]);

    long double alpha = 1 - gama;
    Graph *graph = new Graph(fileName.c_str(), 0);
    cout << fileName.substr(fileName.rfind('/')) << endl;
    graph->read_graph_binary();
    int m = graph->m;
    cout << "Graph Infomation: " << "n: " << graph->n << " m: " << m << endl;
    int n = (1 + sqrt(1 + 8 * m)) / 2;
    cout << "gama: " << gama << endl;
    int RES_LB = 0, ub = n;
    ui *peel_sequence;
    degen(graph, gama, RES_LB, ub, result, peel_sequence);
#ifndef NHEU
    ego_degen(graph, gama, result, peel_sequence);
#endif
    n = ub;
    long long k = -1;
    n = result.size();
    printf("[Get LB:%d]\n", n);
    print_result(result);
    function<int(int, int, int, int, vector<ui> &)> get_defect = [&](int n, int k, int lb, int ub, vector<ui> &result)
    {
        int res = f(k, lb, ub, result);
        cout << " find res:" << res << endl;
        return res;
    };

    static_g = Graph(fileName.c_str(), 0);
    static_g.read_graph_binary();
    ui *out_mapping = new ui[static_g.n];
    ui *to_static_mapping = new ui[static_g.n];
    for (int i = 0; i < static_g.n; i++)
        out_mapping[i] = i, to_static_mapping[i] = i;
    int lb = n;
    int k_cur = my_floor((n - 1) * n / 2 * alpha);
    bool n_increase_slow = (n <= my_floor(1 / alpha + 0.5));
    bool do_first_heuris_n = true; // at least do once because the current result dose not guarrantee a maximal kDC in one hop neighbor subgraph.
#ifndef NHEU
    double before_local_search = clock();
    extend(alpha, k_cur, result, to_static_mapping, out_mapping, n, do_first_heuris_n);
    cout << "Local Search Time: " << (double)(clock() - before_local_search) / CLOCKS_PER_SEC << "s" << endl;
#endif
    cout << "Start Binary With LB:" << result.size() << endl;
    print_result(result);
    assert(check_sol(&static_g, gama, result));
    // return 0;
    lb = n;
    ub = (1 + sqrt(1 + 8 * m / gama)) / 2;
    ;
    k = alpha * (ub - 1) * ub / 2;
    int max_color = 0, *color_cnt;
    greedyColoring(&static_g, peel_sequence, color_cnt, max_color);
    ub = computeUB(k, color_cnt, max_color);
    cout << "Upper Bound: " << ub << endl;
    while (true)
    {
        long long new_k = alpha * (ub - 1) * ub / 2;
        cout << "Last k: " << k << " new k: " << new_k << " ub:" << ub << endl;
        if (new_k >= k)
            break;
        k = new_k;
        ub = computeUB(k, color_cnt, max_color);
    }
    cout << "Color Count: " << max_color << endl;
    for (int i = 0; i < max_color; i++)
        cout << color_cnt[i] << " ";
    delete[] color_cnt;
    vector<ui> deg;
    for (int i = 0; i < static_g.n; i++)
        deg.push_back(static_g.pstart[i + 1] - static_g.pstart[i]);
    sort(deg.begin(), deg.end());
    reverse(deg.begin(), deg.end());
    ui deg_lb = 1, deg_ub = ub;
    while (deg_ub != deg_lb)
    {
        ui mid = (deg_ub + deg_lb + 1) / 2;
        ui sum = 0;
        for (int i = 0; i < mid; i++)
            sum += min(deg[i], mid - 1);
        // cout<<"Mid: "<<mid<<" Sum: "<<sum<<" Need: "<<my_floor(gama*mid*(mid-1))<<endl;
        if (sum >= my_floor(gama * mid * (mid - 1)))
            deg_lb = mid;
        else
            deg_ub = mid - 1;
    }

    cout << "Deg UB: " << deg_ub << endl;
    if (deg_ub < ub)
        ub = deg_ub;

    while (ub > lb)
    {
        // cout<<"UB: "<<ub<<" K: "<<k<<endl;
        ub = get_defect(ub, k, lb, ub, result);
        long long new_k = alpha * (ub - 1) * ub / 2;
        if (new_k == k)
            break;
        k = new_k;
    }
    print_result(result);

    delete[] out_mapping;
    delete[] peel_sequence;
    delete[] to_static_mapping;

    long long s2 = clock();
    cout << "Running Time: " << (double)(s2 - s1) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Search Time: " << (double)searchTime / CLOCKS_PER_SEC << "s" << endl;
    cout << "res: " << ub << endl;
    return 0;
}
