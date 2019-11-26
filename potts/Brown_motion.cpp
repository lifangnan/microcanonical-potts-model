#include<iostream>
#include<fstream>
#include<stdarg.h>
#include<vector>
#include<math.h>
#include<string>
#include<random>
#include<time.h>

using namespace std;

class vec{
    static int _nDim;


public:
    vector<double> _data;
    vec();
    vec(int count, ...);
    void nDim(int n);
    int nDim() const;
    double norm();
    void showvec();

    double& operator[](int);
    double operator[](int) const;

    vec operator+(const vec& v) const;
    vec operator-(const vec& v) const;
    vec& operator+=(const vec& v);
    vec& operator-=(const vec& v);
    vec operator*(double c) const;

    friend vec operator*(const double c,const vec& v);
    //vec operator*(const vec& v) const;
};

int vec::_nDim=3;
vec::vec(){
    _data.assign(_nDim,0.0);
}

vec::vec(int count, ...){
    va_list pvar;
    va_start(pvar,count);
    for(int i=0;i<count;i++){
        _data.insert(_data.begin()+i,va_arg(pvar,double));
    }
    va_end(pvar);
    _nDim = count;
}

void vec::nDim(int n){
    if(n >=_nDim){
        for(int i=0;i<(n-_nDim);i++){
            _data.push_back(0.0);
        }
    }else{
        for(int i=0;i<(_nDim-n);i++){
            _data.pop_back();
        }
    }
    _nDim = n;
}

int vec::nDim() const{
    return _nDim;
}

double vec::norm(){
    double sum = 0.0;
    vector<double>::iterator it = _data.begin();
    for(it;it!=_data.end();it++){
        sum += (*it)*(*it);
    }
    return sqrt(sum);
}

void vec::showvec(){
    cout<<"The vector is:(";
    vector<double>::iterator it = _data.begin();
    for(it;it!=_data.end();it++){
        cout<<*it<<" ";
    }
    cout<<")"<<endl;
}

vec vec::operator+(const vec& v)const{
    vec vec1;
    vec1.nDim(v.nDim());
    for(int i=0;i<this->_nDim;++i){
        vec1._data[i] = this->_data[i]+v._data[i];
    }
    return vec1;
}

vec vec::operator-(const vec& v)const{
    vec vec1;
    vec1.nDim(v.nDim());
    for(int i=0;i<this->_nDim;++i){
        vec1._data[i] = this->_data[i]-v._data[i];
    }
    return vec1;
}

vec& vec::operator+=(const vec& v){
    for(int i=0;i<this->_nDim;++i){
        this->_data[i] = this->_data[i]+v._data[i];
    }
}

vec& vec::operator-=(const vec& v){
    for(int i=0;i<this->_nDim;++i){
        this->_data[i] = this->_data[i]-v._data[i];
    }
}

vec vec::operator*(double c) const{
    vec vec1;
    vec1.nDim(this->nDim());
    for(int i=0;i<this->_nDim;++i){
        vec1._data[i] = this->_data[i]*c;
    }
    return vec1;
}

vec operator*(const double c,const vec& v){
    vec vec1;
    vec1.nDim(v.nDim());
    for(int i=0;i<v.nDim();++i){
        vec1._data[i] = v._data[i]*c;
    }
    return vec1;
}

default_random_engine e((int)time(0));
normal_distribution<float> generate_randnum(0,1);

class particle
{
private:
    /* data */
    float time = 0;
    
public:
    float m,gamma,T;
    float k = 1.38*pow(10,-23);
    float delta_t = 0.0005;
    vec position;
    vec velocity;
    particle(float mm);
    void initial_position(float xx,float yy);
    void initial_velocity(float vx,float vy);
    void set_environment(float ggamma,float TT);
    float dist_to_origin();
    void move_for_delta_t();
    float number_velocity();
};

particle::particle(float mm)
{
    m = mm;
}

void particle::initial_position(float xx,float yy){
    position = vec(2,xx,yy);
}
void particle::initial_velocity(float vx,float vy){
    velocity = vec(2,vx,vy);
}


void particle::set_environment(float ggamma,float TT){
    gamma = ggamma;
    T = TT;
}

float particle::dist_to_origin(){
    return position.norm();
}

float particle::number_velocity(){
    return velocity.norm();
}

void particle::move_for_delta_t(){
    // /*use 4-order Runge_Kutta method to calculate*/
    // vec Sv1 = acceleration(position,velocity);
    // vec Sp1 = velocity;
    // vec mid_velo = velocity + delta_t/2*Sv1;
    // vec mid_posi = position + delta_t/2*Sp1;
    // vec Sv2 = acceleration(mid_posi,mid_velo);
    // vec Sp2 = mid_velo;
    // mid_velo = velocity + delta_t/2*Sv2;
    // mid_posi = position + delta_t/2*Sp2;
    // vec Sv3 = acceleration(mid_posi,mid_velo);
    // vec Sp3 = mid_velo;
    // vec end_velo = velocity + delta_t*Sv3;
    // vec end_posi = position + delta_t*Sp3;
    // vec Sv4 = acceleration(end_posi,end_velo);
    // vec Sp4 = end_velo;
    // end_velo = velocity + delta_t/6*(Sv1+2*Sv2+2*Sv3+Sv4);
    // end_posi = position + delta_t/6*(Sp1+2*Sp2+2*Sp3+Sp4);
    // position = end_posi;
    // velocity = end_velo;
    float rand1 = generate_randnum(e);

    vec random_move(2,generate_randnum(e),generate_randnum(e));
    vec acceleration = (1/m/(1+0.5*gamma*delta_t))*(random_move*sqrt(2*m*gamma*k*T/delta_t)-m*gamma*velocity);

    vec new_velo = velocity + delta_t*acceleration;
    vec new_posi = position + delta_t*0.5*(velocity+new_velo);
    velocity = new_velo;
    position = new_posi;
}

/*separate line------------------------*/

int main(){
    ofstream opt;
    opt.open("data.csv");

    float mass = 1;
    float gamma = 1;
    float T = 300;

    particle brown(mass);
    brown.set_environment(gamma,T);
    brown.initial_position(0,0);
    brown.initial_velocity(0,0);

    opt<<"x,y,v_x,v_y,distance_to_origin,velocity"<<endl;
    for(int i=0;i<10000;++i){
        //cout<<i;
        opt<<brown.position._data[0]<<","<<brown.position._data[1]<<\
        ","<<brown.velocity._data[0]<<","<<brown.velocity._data[1]<<\
        ","<<brown.dist_to_origin()<<","<<brown.number_velocity()<<endl;
        //cout<<brown.dist_to_origin()<<endl;
        brown.move_for_delta_t();
    }

    return 0;
}