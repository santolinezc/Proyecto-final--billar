#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>
#include "Auxiliar.h"

using namespace Eigen;

// Constantes
const int N=1;         // Numero de bolas
const int dim=2;       // Dimension
const double DT=1e-3;  // Dt
const double lx=10.00; // Longitud mesa rectangular en x
const double ly=10.00; // Longitud mesa rectangular en y
const int pasos=1000;   // Numero de iteraciones
const double rad=1e-1; // Radio pelotas
const double alpha=0;  // Parametro de deforamcion mesa estadio
const double R=3.0;    // Radio de la mesa estadio


//Declaracion de funciones
void set_table1(Body billar[]);
void set_conditions_table1(Body billar[]);
void set_table2(Body billar[]);
void set_conditions_table2(Body billar[]);
void timestep_all(Body billar[], double dt);
void iniciar_cuerpos(Body billar[]);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
  std::ofstream fout("datos.dat");
  srand(0);
  Body cuerpo[N];
  iniciar_cuerpos(cuerpo);
  set_conditions_table2(cuerpo);

  Plot Plot;
  Plot.init_gnuplot();
  
  for(int ii=0; ii < pasos; ++ii){
    timestep_all(cuerpo,DT);
    fout<< ii*DT << " " << cuerpo[0].r(0) << " " << cuerpo[0].r(1)
    	<< " " << cuerpo[0].v(0) << " " << cuerpo[0].v(1) << std::endl; 
    set_table2(cuerpo);
    Plot.print_gnuplot(cuerpo);
    }
  
    fout.close();
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
  for(int i= 0;i<N; ++i){
    Vector2d L;
    L<< lx,ly;
    Vector2d Y = billar[i].r;
    for (int ii = 0; ii < Y.size();++ii){
      billar[i].r(ii) =L(ii)*double(rand())/RAND_MAX;
      billar[i].v(ii) = 500*(2*double(rand())/RAND_MAX-1);
      billar[i].m = 1+double(rand())/RAND_MAX;
      billar[i].r2old(ii) = 0;
      billar[i].F(ii) = 0;
      billar[i].E = 0;
      billar[i].Rad = rad; 
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
      delta = fabs(LX) - (R-billar[ii].Rad);
      if(delta > 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  delta = fabs(billar[ii].r(0)) - (R-billar[ii].Rad);
	  
	}while(delta < 0);
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
      delta =rn.norm() - (R-billar[ii].Rad);
      if(delta > 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  if(billar[ii].r(1) < 0){
	    rn = billar[ii].r + Y;
	  }else{
	    rn = billar[ii].r-Y;
	  } 
    	  delta =rn.norm() - (R-billar[ii].Rad);
	}while(delta< 0);
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
  for (int i = 0; i < N;++i){
    Vector2d L,Y,rnew,Diff;
    L << R,R+alpha;
    rnew = billar[i].r;
    for(int ii = 0; ii < rnew.size(); ++ii){
      billar[i].r(ii) =L(ii)*(2*double(rand())/RAND_MAX-1);
      billar[i].v(ii) = 500*(2*double(rand())/RAND_MAX-1);
      billar[i].m = 1+double(rand())/RAND_MAX;
      billar[i].r2old(ii) = 0;
      billar[i].F(ii) = 0;
      billar[i].E = 0;
      billar[i].Rad = rad; 
    }
    Y << 0, alpha;
    Diff = rnew -Y;
    if(Diff.norm() > R){
      for(int ij = 0;ij < Y.size(); ++ij){
	billar[i].r(ij) = rnew(ij)-Y(ij)*(1+double(rand())/RAND_MAX); 
      }
    }
    billar[i].rold_inicial(DT);
  }
}
//implements timestep for all bodies
void timestep_all(Body billar[],double dt)
{
  for (int ii = 0; ii < N; ++ii){
    billar[ii].timestep(dt);
  }
}

//iniciates all bodies
void iniciar_cuerpos(Body billar[])
{
  for (int i = 0; i < N; ++i){
    billar[i].b_init();
  }
}
