#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include "Auxiliar.h"

using namespace Eigen;


// Funciones de la estructura
void Body::b_init()
{
  for (int ii = 0; ii < r.size(); ++ii){
    r(ii) = 0, rold(ii) = 0, r2old(ii) = 0;
    v(ii) = 0, vold(ii) = 0; F(ii) = 0;
  }
  m = 0, E = 0, Rad = 0;
}
void Body::stepback()
{
  for (int ii = 0; ii < r.size(); ++ii){
    r(ii) = rold(ii);
    rold(ii) = r2old(ii);
    v(ii) = vold(ii);
  }
}
void Body::rold_inicial(double dt)
{
  int ii;
  for(ii=0; ii < r.size(); ++ii){
    rold(ii) = r(ii) - dt*v(ii) + F(ii)*dt*dt/(2*m);
  }
}
void Body::timestep(double dt)
{
  int ii;
  double w;
  for(ii=0; ii < r.size(); ++ii){
    w = r(ii);
    r2old(ii) = rold(ii);
    vold(ii) = v(ii);
    r(ii) = dt*v(ii)+rold(ii);
    rold(ii) = w;
  }
}

//sets up gnuplot printing
Plot::Plot()
{
}


void Plot::init_gnuplot(void)
{
  std::cout << "unset key" <<std::endl;
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange [0:1]" << std::endl;
  std::cout << "set term gif animate delay -0.5" << std::endl;
  std::cout << "set output 'Billar.gif' " << std::endl;
    
}
/*
void Plot::print_table1(void)
{
  std::cout << "set xrange [-0.5-" << rad << ":" << lx+rad+0.5 << "]" << std::endl;//for square table
  std::cout << "set yrange [-0.5-" << rad << ":" << ly+rad+0.5 << "]" << std::endl;
  std::cout << "plot " << lx+rad << ","<< -rad<< "+" <<ly +2*rad<<"* t , " ;
  std::cout << -rad << "+" << lx+2*rad << "*t ," << ly+rad << "," ;
  std::cout << 0-rad << "," << -rad << "+" << ly+2*rad << "*t , " ;
  std::cout << -rad << "+" << lx+2*rad << "*t ," << 0-rad << std::endl;
  }*/

void Plot::print_table2(double R, double alpha)
{
  std::cout << "set xrange [" <<-0.5 - R << ":" << R <<  "]" << std::endl;// for stadium
  std::cout << "set yrange [" <<-0.5-alpha-R <<":"<< alpha+R << "]" << std::endl;
  std::cout << "plot " << R << "*cos(pi*t),"<< alpha << "+" << R << "*sin(pi*t)," ;
  std::cout << R << "*cos(pi*t),-("<< alpha << "+" << R << "*sin(pi*t))," ;
  std::cout << R << "," << -alpha << "+ "<<  2*alpha << "*t , " ;
  std::cout << -R << "," << - alpha << "+" << 2*alpha << "*t " <<std::endl ;
}

void Plot::print_gnuplot(Body billar[])
{
  std::cout <<"plot "; 
  for (int ii = 0; ii < billar[ii].r.size(); ++ii){
    std::cout << billar[ii].r(0) << " + " << billar[ii].Rad << "*cos(2*pi*t),"
	      << billar[ii].r(1) << " + " << billar[ii].Rad << "*sin(2*pi*t)  lt " << ii+1 << ",";
  }
  std::cout << " 0, 0";
  std::cout << std::endl;
}
