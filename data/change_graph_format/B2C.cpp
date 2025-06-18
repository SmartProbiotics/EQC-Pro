//处理点的坐标从1开始（都怪real world graphs）处理后下标从0开始 ./C-\>A grpah save_directory

#include<bits/stdc++.h>
using namespace std;
#define pii pair<int,int>
int main(int argc,char* argv[]){
    ofstream of;
    string fileName = string(argv[1]);
    freopen(fileName.c_str(),"r",stdin);
    int idx=0;
    for(int i=fileName.size()-1;i>=0;i--){
        if(fileName[i]=='/'){
            idx=i+1;
            break;
        }
    }
    
    fileName=string(argv[2])+"/"+fileName.substr(idx);
    cout<<fileName<<endl;
    freopen(fileName.c_str(),"w",stdout);

    int n,m;cin>>n>>m;
    cout<<"p edge "<<n<<" "<<m<<endl;
    vector<vector<int>> e(n+10);
    vector<set<int>> st(n+10);
    int dec=0;
    for(int i=0;i<m;i++){
        int v,u;cin>>u>>v;
        cout<<"e "<<u<<" "<<v<<endl;
    }
    of.close();
}