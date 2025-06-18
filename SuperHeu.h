#include "Graph.h"
#include <bits/stdc++.h>
using namespace std;
// #define Log
class LinearSet
{
private:
    ui *data;
    ui id;
    ui *pos;

public:
    LinearSet(ui len, ui max_value)
    {
        data = new ui[len];
        pos = new ui[max_value + 1];
        memset(pos, -1, sizeof(ui) * (max_value + 1));
        id = 0;
    }
    bool count(ui v) { return pos[v] != -1; }
    ~LinearSet()
    {
        delete[] data;
        delete[] pos;
    }
    void clear()
    {
        for (int i = 0; i < id; i++)
            pos[data[i]] = -1;
        id = 0;
    }

    void insert(ui v)
    {
        data[id] = v;
        pos[v] = id++;
    }

    void erase(ui v)
    {
        if (pos[v] == -1)
            return;
        ui u = data[--id];
        pos[u] = pos[v];
        data[pos[u]] = u;
        pos[v] = -1;
    }

    bool empty() const
    {
        return id == 0;
    }

    ui *begin()
    {
        return data;
    }

    ui *end()
    {
        return data + id;
    }
};
class HeuExtractor
{
public:
    ept *pstart; // 指向每个顶点邻接表的起始位置
    ui *edges;   // 存储所有边的数组
    ui m, n;     // 图中边数、点数
    int k;       // 允许缺失边的阈值
    double gamma;
    std::vector<ui> best_res; // 当前得到的最优解集合
    double sol_density = 0;

    HeuExtractor(const Graph &other, const std::vector<ui> &result, double gamma)
        : best_res(result), n(other.n), m(other.m), pstart(other.pstart), edges(other.edges), k(other.K), gamma(gamma)
    {
    }
    void extract()
    {
        vector<ui> init_best_res = best_res;
        for (int i = 0; i < init_best_res.size(); i++)
        {
            vector<ui> tmp = {init_best_res[i]};
            init(tmp);
            extract_once_v2();
        }
    }

    void extract_once_v2()
    {
        while (miss_edge <= k)
        {
            if (idS > (int)best_res.size())
            {
                best_res.clear();
                for (int j = 0; j < idS; j++)
                    best_res.push_back(S[j]);
                double a = (idS + 1) * idS / 2 * (1 - gamma);
                int l = floor(a), r = ceil(a);
                if (abs(a - r) < 1e-6)
                    k = r;
                else
                    k = l;
                sol_density = (double)((idS - 1) * idS / 2 - miss_edge) / ((idS - 1) * idS / 2);
                cout << "Better Solution:" << best_res.size() << " Missing Edges:" << miss_edge << " Next Allow K:" << k << endl;
            }
            getBestV();
            if (bestV == (ui)-1)
                break;
            moveToSet(bestV, S, idS);
            // updateScore(); //v3 without this line
            getBestV();
            if (bestV == (ui)-1)
                break;
            moveToSet(bestV, S, idS);
            getWorstV();
            moveToSet(worstV, C, idC);
            updateScore();
        }
    }

    ~HeuExtractor()
    {
        delete[] score;
        delete[] degInS;
        delete[] degree;
        delete[] S;
        delete[] C;
        delete[] mC;
        delete[] posS;
        delete[] posC;
    }

private:
    ui *S = nullptr, *C = nullptr, *mC = nullptr;
    ui *degree = nullptr, *degInS = nullptr;
    int *posS = nullptr, *posC = nullptr;
    long long *score = nullptr;
    LinearSet *neig = nullptr;
    int idS = 0, idC = 0, idmC = 0, miss_edge = 0;
    ui timeStamp = 0, worstV = 0, bestV = 0, max_deg = 0;

    void updateScore()
    {
        for (int i = 0; i < idC; i++)
        {
            ui u = C[i];
            score[u] += degInS[u];
        }
        for (int i = 0; i < idS; i++)
        {
            ui u = S[i];
            score[u] -= max_deg - (degree[u] - degInS[u]);
        }
    }

    void getWorstV()
    {
        ui idx = 0, mi_deg = n;
        for (int i = 0; i < idS; i++)
        { // rule 1:min_deg in S
            ui u = S[i];
            if (mi_deg > degInS[u])
            {
                idx = 0;
                mC[idx++] = u;
                mi_deg = degInS[u];
            }
            else if (mi_deg == degInS[u])
                mC[idx++] = u;
        }
        ui len = 0;
        ui mi_score = 0 + 1;
        for (int i = 0; i < idx; i++)
        { // rule 1:min_score in S
            ui u = mC[i];
            if (mi_score > score[u])
            {
                len = 0;
                mC[len++] = u;
                mi_score = score[u];
            }
            else if (mi_score == score[u])
                mC[len++] = u;
            ;
        }
        assert(len != 0);
        worstV = mC[rand() % len];
    }

    ui addRule(ui *st, ui len)
    {
        ui idx = 0, mx_deg = 0;
        for (int i = 0; i < len; i++)
        { // rule 1:max_deg in S
            ui u = st[i];
            if (mx_deg < degInS[u])
            {
                idx = 0;
                st[idx++] = u;
                mx_deg = degInS[u];
            }
            else if (mx_deg == degInS[u])
                st[idx++] = u;
        }
        len = 0;
        ui mx_score = 0;
        for (int i = 0; i < idx; i++)
        { // rule 2:max_score in S
            ui u = st[i];
            assert(posC[u] != -1);
            if (mx_score < score[u])
            {
                len = 0;
                st[len++] = u;
                mx_score = score[u];
            }
            else if (mx_score == score[u])
                st[len++] = u;
        }
        assert(len != 0);
        return st[rand() % len];
    }

    void getBestV()
    {
        idmC = 0;
        for (auto u : *neig)
            mC[idmC++] = u;
        if (!idmC)
        {
            bestV = (ui)-1;
            return;
        } // all v in S with degInC[u]==0
        bestV = addRule(mC, idmC);
    }

    inline bool tryAdd(ui &u)
    {
        return (miss_edge + idS - degInS[u] <= k);
    }

    void removeFromSet(ui u, ui *st, int &idx)
    {
        if (idx == 0)
            return;
        int *pos = (st == S) ? posS : posC;
        int curPos = pos[u];
        if (curPos == -1)
            return;
        bool fromS = st == S;
        st[curPos] = st[idx - 1];
        pos[st[idx - 1]] = curPos;
        idx--;

        if (fromS)
        {
            for (ept i = pstart[u]; i != pstart[u + 1]; i++)
            {
                ui v = edges[i];
                degInS[v]--;
                if (!degInS[v] && neig->count(v))
                    neig->erase(v);
            }
            miss_edge -= idS - degInS[u];
        }
        pos[u] = -1;
    }

    void moveToSet(ui u, ui *st, int &idx)
    {
#ifdef Log
        if (st == S)
            cout << "Add Vertice: " << u << " Score:" << score[u] << endl;
#endif
        bool toS = st == S;
        st[idx] = u;
        int *pos = (toS) ? posS : posC;
        pos[u] = idx;
        idx++;
        ui *other = (toS) ? C : S;
        int &other_idx = (toS) ? idC : idS;
        removeFromSet(u, other, other_idx);
        for (ept i = pstart[u]; i != pstart[u + 1]; i++)
        {
            ui v = edges[i];
            if (toS)
                degInS[v]++;
        }
        if (toS)
        {
            miss_edge += idS - degInS[u] - 1;
            if (neig->count(u))
                neig->erase(u);
        }
        else
        { // To C
            if (degInS[u])
                neig->insert(u);
        }
        score[u] = 0;
    }

    void init(vector<ui> &initS)
    {
        miss_edge = 0;
        if (!score)
            score = new long long[n]();
        else
            memset(score, 0, sizeof(long long) * n);
        if (!degree)
            degree = new ui[n];
        if (!degInS)
            degInS = new ui[n]();
        else
            memset(degInS, 0, sizeof(ui) * n);
        if (!S)
            S = new ui[n]();
        else
            memset(S, 0, sizeof(ui) * n);
        if (!C)
            C = new ui[n]();
        else
            memset(C, 0, sizeof(ui) * n);
        if (!posS)
            posS = new int[n];
        else
            memset(posS, -1, sizeof(int) * n);
        if (!posC)
            posC = new int[n];
        else
            memset(posC, -1, sizeof(int) * n);
        if (!mC)
            mC = new ui[n]();
        else
            memset(mC, 0, sizeof(ui) * n);
        if (!neig)
            neig = new LinearSet(n, n);
        else
            neig->clear();
        for (int i = 0; i < n; i++)
            degree[i] = pstart[i + 1] - pstart[i];
        idS = idC = 0;
        max_deg = 0;
        for (int i = 0; i < n; i++)
            max_deg = max(max_deg, pstart[i + 1] - pstart[i]);
        for (auto u : initS)
            moveToSet(u, S, idS);
        for (int i = 0; i < n; i++)
            if (posS[i] == -1)
                moveToSet(i, C, idC);
    }
};