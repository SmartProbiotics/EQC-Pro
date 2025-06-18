#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]) {
    string fileName = string(argv[1]);
    freopen(fileName.c_str(), "r", stdin);

    int idx = 0;
    for (int i = fileName.size() - 1; i >= 0; i--) {
        if (fileName[i] == '/') {
            idx = i + 1;
            break;
        }
    }

    fileName = string(argv[2]) + "/" + fileName.substr(idx);
    cout << fileName << endl;
    freopen(fileName.c_str(), "w", stdout);

    string p; cin >> p;
    string edges; cin >> edges;
    
    int n, m; cin >> n >> m;
    vector<vector<int>> e(n); // È¥µô+10£¬Ò²È¥µôset

    cout << n << " " << m << endl;

    int dec = 0;
    for (int i = 0; i < m; i++) {
        char c; cin >> c;
        int u, v; cin >> u >> v;
        if (i == 0) dec = (u != 0 && v != 0);
        u -= dec;
        v -= dec;
        if(v>=n||u>=n) continue;
        e[u].push_back(v);
        e[v].push_back(u);
    }

    for (int i = 0; i < n; i++) {
        cout << i;
        if (!e[i].empty()) cout << " ";
        for (int j = 0; j < e[i].size(); j++) {
            cout << e[i][j] << (j == e[i].size() - 1 ? "" : " ");
        }
        cout << endl;
    }

    return 0;
}
