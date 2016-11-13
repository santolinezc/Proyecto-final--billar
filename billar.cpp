#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;
using Eigen::MatrixXd;
///constantes
const int N=1;//numero de bolas
const int dim=2;//dimension
const double DT=1e-2;

//estructura de datos para el cuerpo
struct body {
  Vector2d r,v,F,rold,r2old,vold;
  double m,E;
  void b_init();
  void timestep(double dt);
  void rold_inicial(double dt);
  void stepback();
};

//funciones de la estructura
void body::b_init()
{
  for(int i=0;i<r.size();++i){
    r(i)=0,rold(i)=0,r2old(i)=0;
    v(i)=0,vold(i)=0,F(i)=0;
  }
  m=1;
  E=0;
}
void body::stepback()
{
  for (int ii = 0;ii<r.size();++ii){
    r(ii) = rold(ii);
    rold(ii) = r2old(ii);
    vold(ii) = v(ii);
  }
}
void body::rold_inicial(double dt)
{
  int ii;
  for(ii=0;ii<r.size();++ii){
    rold(ii) = r(ii) - dt*v(ii) + F(ii)*dt*dt/(2*m);
  }
}
void body::timestep(double dt)
{
  int ii;
  double w;
  for(ii=0;ii<r.size();++ii){
    w = r(ii);
    r2old(ii) = rold(ii);
    vold(ii) = v(ii);
    r(ii) = 2*r(ii)-rold(ii)+F(ii)*dt*dt/m;
    v(ii) = (r(ii)-rold(ii))/(2*dt);
    rold(ii) = w;
  }
}

//Declaracion de funciones
void set_table(body  billar);
/////////////////////////
int main()
{
  body cuerpo;
  cuerpo.b_init();
  cuerpo.rold_inicial(DT);
  cuerpo.r<<1,2;
  cuerpo.v<<4,5;
  for(int ii=0;ii<3;++ii){
    cuerpo.timestep(DT);
    std::cout<<"r "<<cuerpo.r<<std::endl<<std::endl;
    std::cout<<"rold "<<cuerpo.rold<<std::endl<<std::endl;
    std::cout<<"v "<<cuerpo.v<<std::endl<<std::endl;
    std::cout<<"vold"<< cuerpo.vold<<std::endl<<std::endl;
  }
  return 0;
}
////////////////////////
//sets  rectangular table
void set_table(body  billar)
{
  //set rectangular table
  double lx = 12.23,ly = 14.24;
  double delta;
  double dt1 = DT/100;
  delta = billar.r(1)-lx;
  if(delta>0){
    billar.stepback();
   
    do{
      billar.timestep(dt1);
      delta = billar.r(1)-lx;
    }while(delta>0);
    
  }
}
