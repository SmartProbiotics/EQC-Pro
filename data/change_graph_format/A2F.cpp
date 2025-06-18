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
    fileName=fileName.substr(idx);
    int tp_idx=fileName.size();
    for(int i=fileName.size()-1;i>=0;i--){
        if(fileName[i]=='.'){
            tp_idx=i;
            break;
        }
    }
    fileName=fileName.substr(0,tp_idx)+".fpcetp";
    fileName=string(argv[2])+"/"+fileName;
    
    cout<<fileName<<endl;
    of.open(fileName.c_str(),ios::out);
    int m,n;cin>>n>>m;
    getchar();
    string s;
    while(getline(cin,s)){
        stringstream ss(s);
        int v;
        ss>>v;
        while(ss>>v) of<<v<<" ";
        of<<endl;
    }

    of.close();
}