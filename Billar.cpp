#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>
#include "Auxiliar.h"

using namespace Eigen;

// Constantes
const int N = 2;         // Numero de bolas
const int dim = 2;        // Dimension
const double lx = 10.00;  // Longitud mesa rectangular en x
const double ly = 10.00;  // Longitud mesa rectangular en y
const int steps = 100000;   // Numero de iteraciones
const double rad = 0;  // Radio pelotas
const double alpha = 0.1; // Parametro de deforamcion mesa estadio
const double R = 0.5;     // Radio de la mesa estadio
const double DT = 1.0/steps;  // Dt

const int choose = 1;     // 1 : Stadium; 0 : Rectangular table

//Declaracion de funciones
void set_table1(Body billar[]);
void set_conditions_table1(Body billar[]);
void set_table2(Body billar[]);
void set_conditions_table2(Body billar[]);
void timestep_all(Body billar[], double dt);
void iniciar_cuerpos(Body billar[]);
void tables(Body billar[]);
void set_conditions(Body billar[]);

// Opciones
const int OP = 1; // OP = 0 Condiciones iniciales aleatorias, OP = 1 Condiciones iniciales para el calculo del coeficiente de Lyapunov


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
  std::ofstream f1out("Data.txt");
  std::ofstream f2out("Lyapunov.txt");
  srand(0);

  
  Body cuerpo[N];
  iniciar_cuerpos(cuerpo);
  set_conditions(cuerpo);

  Plot Plot;
  Plot.init_gnuplot();
  
  for(int ii=0; ii < steps; ++ii){
    timestep_all(cuerpo, DT);

    f1out << ii*DT;
    for(int k = 0; k < N; ++k){
      f1out << " " <<cuerpo[k].r(0) << " " << cuerpo[k].r(1); //Datos de trajectoria de cad particula
    }
    f1out << std::endl;
    if (OP == 1){
      f2out << ii*DT << " " << (cuerpo[0].r - cuerpo[1].r).norm() << std::endl; //Separacion entre la primera y segunda bola
    }
    tables(cuerpo);
    //Plot.print_gnuplot(cuerpo, N, R, alpha, rad, lx, ly, choose);
    }
  
    f1out.close();
    f2out.close();

    Plot.plot_trajectories(N, steps, R, alpha, lx, ly, choose);  //Tomar condiciones tal que DT sea pequeño, sino, el algoritmo fallara
    
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sets  rectangular table
void set_table1(Body  billar[])
{
  double delta,delta1;
  double dt1= DT/100;
  for(int k = 0; k < N;++k){
    // Top and bottom sides
    delta = billar[k].r(1)-ly;
    delta1 = -billar[k].r(1);
    if(delta > 0 or delta1 > 0){
      billar[k].stepback();
      do{
	billar[k].timestep(dt1);
	delta = billar[k].r(1) - ly;
	delta1 = -billar[k].r(1);
      }while(delta < 0 and delta1 < 0);
      billar[k].stepback();
      billar[k].v(1) = -billar[k].v(1);
      billar[k].timestep(dt1);
    }
    // Right and left sides
    delta = billar[k].r(0)-lx;
    delta1 = -billar[k].r(0);

    if(delta > 0 or delta1 > 0){
      billar[k].stepback();
      do{
	billar[k].timestep(dt1);
	delta = billar[k].r(0) - lx;
	delta1 = -billar[k].r(0);
      }while(delta < 0 and delta1 < 0);
      billar[k].stepback();
      billar[k].v(0) = -billar[k].v(0);
      billar[k].timestep(dt1);
    }
  }    
}
// Condiciones iniciales para mesa rectangular
void set_conditions_table1(Body  billar[])
{
  double L[2];
  L[0] = lx;
  L[1] = ly;
  double Random;
  for (int i = 0; i < N; ++i){
    if (OP == 0){
      for(int ii = 0; ii < 2; ++ii){
	Random = double(rand())/RAND_MAX;
	if (Random == 0.0){Random += 0.1;}
	billar[i].r(ii) = L[ii]*Random;
	billar[i].v(ii) = 300*Random;
	billar[i].m = 1 + double(rand())/RAND_MAX;
	billar[i].r2old(ii) = 0;
	billar[i].F(ii) = 0;
	billar[i].E = 0;
	billar[i].Rad = rad; 
      }
      billar[i].rold_inicial(DT);
    }
    else if (OP == 1){
      double Delta = 1e-2;
      for(int ii = 0; ii < 2; ++ii){
	Random = double(rand())/RAND_MAX;
	if (Random == 0.0){Random += 0.1;}
	if (i == 0){	
	  billar[i].r(ii) = L[ii]*Random;
	  billar[i].v(ii) = 1200*Random;
	  billar[i].m = 1 + double(rand())/RAND_MAX;
	  billar[i].r2old(ii) = 0;
	  billar[i].F(ii) = 0;
	  billar[i].E = 0;
	  billar[i].Rad = rad; 
	}
	else {
	  billar[i].r(ii) = billar[0].r(ii) - Delta*i;
	  billar[i].v(ii) = billar[0].v(ii) - Delta*i;
	  billar[i].m = billar[0].m;
	  billar[i].r2old(ii) = 0;
	  billar[i].F(ii) = 0;
	  billar[i].E = 0;
	  billar[i].Rad = rad; 
	}
      }
    }
    billar[i].rold_inicial(DT);
  }
}

// Sets stadium table
void set_table2(Body billar[])
{
  double delta,LX,LY;
  double dt1 = DT/100;
  for( int ii = 0; ii <N; ++ii){
    LY = billar[ii].r(1);
    if(fabs(LY) < alpha ){
      LX = billar[ii].r(0);
      delta = fabs(LX) - (R - billar[ii].Rad);
      if(delta >= 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  delta = fabs(billar[ii].r(0)) - (R - billar[ii].Rad);
	  
	}while(delta <= 0);
	billar[ii].stepback();
	billar[ii].v(0) = -billar[ii].v(0);
	billar[ii].timestep(dt1);
      }
    }else{
      Vector2d Y, rn,runit,vll,vp,rvec; // In order rn:rnew, r unitario, v paralelo, v perpendicular, r vector
      Y << 0,alpha;
      if(billar[ii].r(1) < 0){
	rn = billar[ii].r + Y;
      }else{
	rn = billar[ii].r - Y;
      }
      delta =rn.norm() - (R - billar[ii].Rad);
      if(delta >= 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  if(billar[ii].r(1) < 0){
	    rn = billar[ii].r + Y;
	  }else{
	    rn = billar[ii].r-Y;
	  } 
    	  delta =rn.norm() - (R - billar[ii].Rad);
	}while(delta < 0);
	billar[ii].stepback();
	rvec = billar[ii].r;
	runit = rn.normalized();
	vp = (billar[ii].v.dot(runit))*runit;
	vll = billar[ii].v-vp;
	vp = -vp;
	billar[ii].v = vp+vll;
      }
    }
  }
}
//sets initial contitions for stadium table
void set_conditions_table2(Body billar[])
{
  double Ri;
  double Random;
  for (int i = 0; i < N; ++i){
    if (OP == 0){
      for(int ii = 0; ii < 2; ++ii){
	Random = R*((1.95*double(rand())/RAND_MAX)-1.0);
	if (Random == 0.0){Random += 0.1;}
	if(ii == 0){
	  Ri = abs(Random);
	}
	else{Ri = -Ri + alpha;}
	billar[i].r(ii) = Ri*((1.95*double(rand())/RAND_MAX)-1.0);
	billar[i].v(ii) = 1200*Random;
	billar[i].m = 1 + double(rand())/RAND_MAX;
	billar[i].r2old(ii) = 0;
	billar[i].F(ii) = 0;
	billar[i].E = 0;
	billar[i].Rad = rad; 
      }
      billar[i].rold_inicial(DT);
    }
    else if (OP == 1){
      double Delta = 1e-2;
      for(int ii = 0; ii < 2; ++ii){
	Random = R*((1.95*double(rand())/RAND_MAX)-1.0);
	if (Random == 0.0){Random += 0.1;}
	if(ii == 0){
	  Ri = abs(Random);
	}
	else{Ri = -Ri + alpha;}
	if (i == 0){	
	  billar[i].r(ii) = Ri*((1.95*double(rand())/RAND_MAX)-1.0);
	  billar[i].v(ii) = 1200*Random;
	  billar[i].m = 1 + double(rand())/RAND_MAX;
	  billar[i].r2old(ii) = 0;
	  billar[i].F(ii) = 0;
	  billar[i].E = 0;
	  billar[i].Rad = rad; 
	}
	else {
	  billar[i].r(ii) = billar[0].r(ii) - Delta*i;
	  billar[i].v(ii) = billar[0].v(ii) - Delta*i;
	  billar[i].m = billar[0].m;
	  billar[i].r2old(ii) = 0;
	  billar[i].F(ii) = 0;
	  billar[i].E = 0;
	  billar[i].Rad = rad; 
	}
      }
    }
    billar[i].rold_inicial(DT);
  }
}
  
// Implements timestep for all bodies
void timestep_all(Body billar[],double dt)
{
  for (int ii = 0; ii < N; ++ii){
    billar[ii].timestep(dt);
  }
}

// Iniciates all bodies
void iniciar_cuerpos(Body billar[])
{
  for (int i = 0; i < N; ++i){
    billar[i].b_init();
  }
}

//
void tables(Body billar[])
{
  if (choose == 1){set_table2(billar);}
  else {set_table1(billar);}
}

void set_conditions(Body billar[])
{
 if (choose == 1){set_conditions_table2(billar);}
 else {set_conditions_table1(billar);}

}
