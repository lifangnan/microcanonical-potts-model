#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <random>
#include <set>
#include <time.h>
#include <algorithm>
#include <fstream>

using namespace std;

ofstream opt;
ofstream opt2;

void write_uint(int x){
    int i;
    static char buf[12];
    if(x==0){
        putchar('0');
        return;
    }

    for(i=10; x; i--){
        buf[i] = '0' + x % 10;
        x /= 10;
    }
    fputs(buf+1+i, stdout);
}

void write_uints(int x, int y){
        write_uint(x);
        putchar(' ');
        write_uint(y);
        putchar(',');
}


set<pair<int,int> > RRG(int N, int K){
    int seed;
    set<pair<int,int> > edges;
    vector<int> a;
    vector<int> r;
    size_t os, ns;

    //scanf("%d %d",&N,&K); //n is number of vertexs, d is number of vertexs that each vertex links to.
    if((N*K)&0x1){
        cout<<"Cannot generate RRG, please change N and K"<<endl;
        return edges;
    }

    seed = (int)time(0);
    mt19937 gen(seed);

    do {
        edges.clear();
        a.clear();
        for(int i = 0; i < N*K; i++){
        a.push_back(i%N);
        }
        os = ns = N*K;
        do {
            vector<int>::iterator it;

        //    cout << a.size() << endl;

        //    random_shuffle(a.begin(), a.end());
            shuffle(a.begin(), a.end(), gen);
            for(it = a.begin(); it != a.end(); it+=2){
                pair<int, int> edge(minmax(*it, *(it+1)));
                if(*it == *(it+1) || edges.count(edge)) { r.push_back(*it); r.push_back(*(it+1)); }
                else { edges.insert(edge); }
            }
            a.clear();
            r.swap(a);
            os = ns;
            ns = a.size();
        } while(ns != 0 && ns != os);
    } while(ns != 0);

    // for(set<pair<int, int> >::iterator it = edges.begin(); it != edges.end(); ++it){
    // //    cout << (*it).first << " " << (*it).second << endl;
    //     write_uints((*it).first, (*it).second);
    // }

    // for(int i=0;i<a.size();++i){
    //     cout<<i<<":"<<a[i]<<endl;
    // }
    // for(set<pair<int,int> >::iterator it=edges.begin();it!=edges.end();++it){
    //     cout<<(*it).first<<(*it).second<<endl;
    // }

    return edges;
}

class potts_demon
{
private:
    int N,K,Q;
    set<pair<int,int> > edges;
    vector<int> config;
    vector<set<int> > link;
    
public:
    potts_demon(int NN,int KK,int QQ);
    float energy_density_u();
    void initial_with_u(float std_u);
    void demon_update(int MCstep);
    vector<float> rho_1();
};

potts_demon::potts_demon(int NN,int KK,int QQ)
{
    N = NN;
    K = KK;
    Q = QQ;
    edges = RRG(NN,KK);

    for(int i=0;i<NN;++i){
        set<int> temp;
        temp.insert(-1);
        link.push_back(temp);
    }
    for(set<pair<int,int> >::iterator it = edges.begin();it!=edges.end();it++){
        link[(*it).first].insert((*it).second);
        link[(*it).second].insert((*it).first);
    }
    for(int i=0;i<NN;++i)link[i].erase(-1);
}

float potts_demon::energy_density_u(){
    int total_E = 0;
    float u;
    for(set<pair<int,int> >::iterator it = edges.begin();it!=edges.end();it++){
        if(config[(*it).first]==config[(*it).second])total_E--;
    }
    return (float)total_E/N;
}

void potts_demon::initial_with_u(float std_u){
    config.assign(N,1);
    //for(int i=0;i<N;++i) cout<<config[i]<<" ";
    std::default_random_engine generator1;
    std::uniform_int_distribution<int> posi(0,N-1);
    std::default_random_engine generator2;
    std::uniform_int_distribution<int> color(1,Q);
    while(energy_density_u()<std_u){
        for(int i=0;i<N/1500;++i)config[posi(generator1)] = color(generator2);
       // opt2<<energy_density_u()<<endl;
    }
    cout<<"Initialized with u = "<<energy_density_u()<<endl;
    char filename[20];
    sprintf(filename,"u=%.5f.csv",energy_density_u());
    opt.open(filename,ios::out|ios::trunc);
}

void potts_demon::demon_update(int MCstep){
    float Ed = 0;
    std::default_random_engine gen1;
    std::uniform_int_distribution<int> posi(0,N-1);
    std::uniform_int_distribution<int> color(1,Q);
    for(int j=0;j<MCstep/10;++j){
        for(int i=0;i<N*10;++i){
            int site = posi(gen1);
            int Q2 = color(gen1);
            float localE0 = 0;
            float localE1 = 0;
            for(set<int>::iterator it = link[site].begin();it!=link[site].end();it++){
                if(config[site]==config[*it])localE0--;
                if(Q2==config[*it])localE1--;
            }
            float Ed2 = Ed+localE0-localE1;
            if(Ed2>0){
                Ed = Ed2;
                config[site] = Q2;
            }
        }
	vector<float> rho = rho_1();
        opt<<j*10<<","<<energy_density_u()<<","<<rho[1]<<","<<rho[2]<<","<<rho[3]<<","<<rho[4]<<","<<rho[5]<<","<<rho[6]<<","<<-111<<","<<Ed<<endl;
        //opt2<<j<<" ";
    }
}

vector<float> potts_demon::rho_1(){
    vector<float> rho(Q+1,0);
    for( vector<int>::iterator it=config.begin();it!=config.end();it++){
	rho[(*it)]++;
    }
    return rho;
}


int main(){
    int N = 65536;
    int K = 4;
    int Q = 6;
    
    opt2.open("log.txt");
   
    for(int i=275;i<288;++i){
        float u = -2+i*0.004;

        potts_demon model1(N,K,Q);
        model1.initial_with_u(u);
        model1.demon_update(16000);

        opt.close();
    }
    opt2.close();
    return 0;
}
