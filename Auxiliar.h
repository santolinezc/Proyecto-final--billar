#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>

using namespace Eigen;

class Body
{
 public:
  Vector2d r, v, F, rold, r2old, vold;
  double Rad, m, E;
  void b_init();
  void timestep(double dt);
  void rold_inicial(double dt);
  void stepback(); 
};

class Plot
{
 public:

  Plot();
  
  void init_gnuplot(void);
  void print_table1(void);
  void print_table2(double R, double alpha);
  void print_gnuplot(Body billar[]);
};
