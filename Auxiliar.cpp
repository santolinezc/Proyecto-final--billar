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
  std::cout << "set term gif animate delay 2" << std::endl;
  std::cout << "unset ytics" << std::endl;
  std::cout << "unset xtics" << std::endl;
  std::cout << "unset border" << std::endl;
  std::cout << "set output 'Billar.gif' " << std::endl;
}

void Plot::print_table1(double rad, double lx, double ly)
{
  std::cout << "set xrange [-0.5-" << rad << ":" << lx+rad+0.5 << "]" << std::endl;//for square table
  std::cout << "set yrange [-0.5-" << rad << ":" << ly+rad+0.5 << "]" << std::endl;
  std::cout << "plot " << lx+rad << ","<< -rad<< "+" <<ly +2*rad<<"* t , " ;
  std::cout << -rad << "+" << lx+2*rad << "*t ," << ly+rad << "," ;
  std::cout << 0-rad << "," << -rad << "+" << ly+2*rad << "*t , " ;
  std::cout << -rad << "+" << lx+2*rad << "*t ," << 0-rad << std::endl;
  }

void Plot::print_table2(double R, double alpha)
{
  std::cout << R << "*cos(pi*t),"<< alpha << "+" << R << "*sin(pi*t) lt 2 lw 2, " ;
  std::cout << R << "*cos(pi*t),-("<< alpha << "+" << R << "*sin(pi*t)) lt 2 lw 2," ;
  std::cout << R << "," << -alpha << "+ "<<  2*alpha << "*t lt 2 lw 2, " ;
  std::cout << -R << "," << - alpha << "+" << 2*alpha << "*t lt 2 lw 2; ";
}

void Plot::print_gnuplot(Body billar[], int N, double R, double alpha)
{
  std::cout <<"plot "; 
  for (int ii = 0; ii < N; ++ii){
    std::cout << billar[ii].r(0) << " + " << billar[ii].Rad << "*cos(2*pi*t),"
	      << billar[ii].r(1) << " + " << billar[ii].Rad << "*sin(2*pi*t)  lt " << ii+1 << ",";
  }
  print_table2(R, alpha);  
  std::cout << std::endl;
}

void Plot::plot_trajectories(int M, int steps, double R, double alpha)
{

  std::ofstream fout("Data_trajectories.gp");

  fout << "reset" << std::endl;
  fout << "unset key" << std::endl;
  fout << "set size ratio -1" << std::endl;
  fout << "set parametric" << std::endl;
  fout << "set trange [0:1]" << std::endl;
  fout << "set term gif animate delay 3" << std::endl;
  fout << "unset ytics" << std::endl;
  fout << "unset xtics" << std::endl;
  fout << "unset border" << std::endl;
  fout << "set output 'Billar_trajectories.gif'" << std::endl;

  fout << "counter = 0" << std::endl;
  fout << "M = " << M-1 << std::endl;
  fout << "alpha = " << alpha << std::endl;
  fout << "R = " << R << std::endl;

  fout << "do for[ii = 1:" << steps-3 << ":3 ]{if( ii < " << steps-3 << "/100 ){plot for [n = 0:M] 'Data.txt' using 2+2*n:3+2*n every ::1::ii w l ls n lt n dashtype 2, R*cos(pi*t), alpha + R*sin(pi*t) lt 2 lw 2, R*cos(pi*t), -(alpha + R*sin(pi*t)) lt 2 lw 2, R, -alpha + 2*alpha*t lt 2 lw 2, -R, -alpha + 2*alpha*t lt 2 lw 2, for [n = 0:M] 'Data.txt'  using 2+2*n:3+2*n every ::ii::ii w p ps 2 pt n} else {plot for [n = 0:M] 'Data.txt' using 2+2*n:3+2*n every ::(1+counter)::ii w l ls n lt n, for [n = 0:M] 'Data.txt'  using 2+2*n:3+2*n every::ii::ii w p ps 2 pt n, R*cos(pi*t), alpha + R*sin(pi*t) lt 2 lw 2, R*cos(pi*t), -(alpha + R*sin(pi*t)) lt 2 lw 2, R, -alpha + 2*alpha*t lt 2 lw 2, -R, -alpha + 2*alpha*t lt 2 lw 2; counter = counter + 1}}" << std::endl;
  fout.close();

  std::cout << "load 'Data_trajectories.gp' " << std::endl;

}
