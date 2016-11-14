#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>

using namespace Eigen;
using Eigen::MatrixXd;

///constantes
const int N=1;//numero de bolas
const int dim=2;//dimension
const double DT=1e-3;
const double lx=2.00;
const double ly=1.00;
const int pasos=500;
const double rad=1e-3;
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
    //r(ii) = 2*r(ii)-rold(ii)+F(ii)*dt*dt/m;
    /*v(ii) = (r(ii)-rold(ii))/(2*dt);*/
      r(ii) = dt*v(ii)+rold(ii);
    rold(ii) = w;
  }
}

//Declaracion de funciones

void set_table(body &  billar);
void set_conditions(body & billar);

void init_gnuplot(void);
void print_table(void);
void print_gnuplot(body billar);

/////////////////////////
int main()
{
  std::ofstream fout("datos.dat");
  srand(0);
  body cuerpo;
  cuerpo.b_init();
  set_conditions(cuerpo);
  cuerpo.rold_inicial(DT);
  init_gnuplot();
  
  /*std::cout << cuerpo.r << std::endl;
  std::cout << cuerpo.v << std::endl;
  set_conditions(cuerpo);
  std::cout << cuerpo.r << std::endl;
  std::cout << cuerpo.v << std::endl;
  cuerpo.rold_inicial(DT);
  cuerpo.timestep(DT);
  std::cout << cuerpo.r << std::endl;
  std::cout << cuerpo.v << std::endl;
  cuerpo.timestep(DT);
  std::cout << cuerpo.r << std::endl;
  std::cout << cuerpo.v << std::endl;*/

  for(int ii=0;ii<pasos;++ii){
    cuerpo.timestep(DT);
    fout<< ii*DT << " " << cuerpo.r(0) << " " << cuerpo.r(1)
	<< " " << cuerpo.v(0) << " " << cuerpo.v(1) << std::endl; 
    set_table(cuerpo);
    print_gnuplot(cuerpo);
    }
  
  fout.close();
  return 0;
}
////////////////////////
//sets  rectangular table
void set_table(body  &billar)
{
  double delta,delta1;
  double dt1= DT/100;
  // top side
  delta = billar.r(1)-ly;
  delta1 = -billar.r(1);
  if(delta > 0){
    billar.stepback();
    do{
      billar.timestep(dt1);
      delta = billar.r(1) - ly;
    }while(delta < 0);
    billar.stepback();
    billar.v(1) = -billar.v(1);
    billar.timestep(dt1);
  }
  //bottom side
  if(delta1 > 0){
    billar.stepback();
    do{
      billar.timestep(dt1);
      delta1 = -billar.r(1);
    }while(delta1 < 0);
    billar.stepback();
    billar.v(1) = -billar.v(1);
    billar.timestep(dt1);
  }
  //right side
   delta = billar.r(0)-lx;
  delta1 = -billar.r(0);

  if(delta > 0){
    billar.stepback();
    do{
      billar.timestep(dt1);
      delta = billar.r(0) - lx;
    }while(delta < 0);
    billar.stepback();
    billar.v(0) = -billar.v(0);
    billar.timestep(dt1);
  }
  //left side
  if(delta1 > 0){
    billar.stepback();
    do{
      billar.timestep(dt1);
      delta1 = -billar.r(0);
    }while(delta1 < 0);
    billar.stepback();
    billar.v(0) = -billar.v(0);
    billar.timestep(dt1);
  } 
    
}
//condiciones iniciales
void set_conditions(body & billar)
{
  Vector2d L;
  L<< lx,ly;
  for (int ii = 0; ii < billar.r.size();++ii){
    billar.r(ii) =L(ii)*double(rand())/RAND_MAX;
    billar.v(ii) = 10*double(rand())/RAND_MAX;
    billar.m = 1+double(rand())/RAND_MAX;
  }
}


//sets up gnuplot printing
void init_gnuplot(void)
{
  std::cout << "unset key" <<std::endl;
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange [0:1]" << std::endl;
  std::cout << "set xrange [-0.5:" << lx+0.5 << "]" << std::endl;
  std::cout << "set yrange [-0.5:" << ly+0.5 << "]" << std::endl;
  print_table();
}
void print_table(void)
{
  std::cout << "plot " << lx << ","<< ly <<"* t , " ;
  std::cout << lx << "*t ," << ly << "," ;
  std::cout << 0 << "," <<  ly << "*t , " ;
  std::cout << lx << "*t ," << 0 << std::endl;
}
void print_gnuplot(body billar)
{
  std::cout <<"replot "; 
  std::cout << billar.r(0) << " + " << rad << "*cos(2*pi*t),"
            << billar.r(1) << " + " << rad << "*sin(2*pi*t) " << std::endl;
    
}
