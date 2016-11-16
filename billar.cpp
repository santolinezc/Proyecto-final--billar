#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>

using namespace Eigen;
//using Eigen::MatrixXd;

///constantes
const int N=1;// numero de bolas
const int dim=2;// dimension
const double DT=1e-3;// dt
const double lx=1.00;// longitud mesa rectangular en x
const double ly=1.00;// longitud mesa rectangular en y
const int pasos=200; // numero de iteraciones
const double rad=1e-2;// radio pelotas
const double alpha=1e-3; //paramttro de deforamcion mesa estadio
const double R=1.0;//radio de la mesa estadio

//estructura de datos para el cuerpo
struct body {
  Vector2d r,v,F,rold,r2old,vold;
  double m,E;
  void timestep(double dt);
  void rold_inicial(double dt);
  void stepback();
};


//funciones de la estructura

void body::stepback()
{
  for (int ii = 0;ii<r.size();++ii){
    r(ii) = rold(ii);
    rold(ii) = r2old(ii);
    v(ii) = vold(ii);
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
    r(ii) = dt*v(ii)+rold(ii);
    rold(ii) = w;
  }
}

//Declaracion de funciones

void set_table1(body   billar[]);
void set_conditions_table1(body  billar[]);
void set_table2(body  billar[]);
void set_conditions_table2(body billar[]);
void timestep_all(body billar[], double dt);

void init_gnuplot(void);
void print_table1(void);
void print_table2(void);
void print_gnuplot(body billar[]);

/////////////////////////
int main()
{
  std::ofstream fout("datos.dat");
  srand(0);
  body cuerpo[N];
  set_conditions_table2(cuerpo); 
  init_gnuplot();
  
  for(int ii=0;ii<pasos;++ii){
    timestep_all(cuerpo,DT);
    fout<< ii*DT << " " << cuerpo[0].r(0) << " " << cuerpo[0].r(1)
    	<< " " << cuerpo[0].v(0) << " " << cuerpo[0].v(1) << std::endl; 
    set_table2(cuerpo);
    print_gnuplot(cuerpo);
    }
  
    fout.close();
  return 0;
}
////////////////////////
//sets  rectangular table
void set_table1(body  billar[])
{
  double delta,delta1;
  double dt1= DT/100;
  for(int k = 0; k < N;++k){
    // top and bottom sides
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
    //right and left sides
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
//condiciones iniciales para mesa rectangular
void set_conditions_table1(body  billar[])
{
  for(int i= 0;i<N; ++i){
    Vector2d L;
    L<< lx,ly;
    Vector2d Y = billar[i].r;
    for (int ii = 0; ii < Y.size();++ii){
      billar[i].r(ii) =L(ii)*double(rand())/RAND_MAX;
      billar[i].v(ii) = 50*double(rand())/RAND_MAX;
      billar[i].m = 1+double(rand())/RAND_MAX;
      billar[i].r2old(ii) = 0;
      billar[i].F(ii) = 0;
      billar[i].E = 0;
    }
    billar[i].rold_inicial(DT);
  }
}
//sets stadium table

void set_table2(body billar[])
{
  double delta,LX,LY;
  double dt1 = DT/100;
  for( int ii = 0; ii <N; ++ii){
    LY = billar[ii].r(1);
    if(fabs(LY) < alpha ){
      LX = billar[ii].r(0);
      delta = fabs(LX)-R;
      if(delta > 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  delta = fabs(billar[ii].r(0)) - R;
	  
	}while(delta < 0);
	billar[ii].stepback();
	billar[ii].v(0) = -billar[ii].v(0);
	billar[ii].timestep(dt1);
      }
    }else{
      Vector2d Y, rn,runit,vll,vp,rvec; // in order rn:rnew, r unitario, v paralelo, v perpendicular, r vector
      Y << 0,alpha;
      rn = billar[ii].r-Y;
      delta =rn.norm()-R;
      if(delta > 0){
	billar[ii].stepback();
	do{
	  billar[ii].timestep(dt1);
	  rn = billar[ii].r-Y;
	  delta =rn.norm() - R;
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
void set_conditions_table2(body billar[])
{
  for (int i = 0; i < N;++i){
    Vector2d L,Y,rnew,Diff;
    L << R,R;
    rnew = billar[i].r;
    for(int ii = 0; ii < rnew.size(); ++ii){
      billar[i].r(ii) =L(ii)*double(rand())/RAND_MAX;
      billar[i].v(ii) = 50*double(rand())/RAND_MAX;
      billar[i].m = 1+double(rand())/RAND_MAX;
      billar[i].r2old(ii) = 0;
      billar[i].F(ii) = 0;
      billar[i].E = 0;
    }
    Y = R*rnew.normalized();
    Diff = rnew -Y;
    if(Diff.norm() > 0){
      for(int ij = 0;ij < Y.size(); ++ij){
	billar[i].r(ij) = rnew(ij)-Y(ij)*(1+double(rand())/RAND_MAX); 
      }
    billar[i].rold_inicial(DT);
    }
  }
}
//implements timestep for all bodies
void timestep_all(body billar[],double dt)
{
  for (int ii = 0; ii < N; ++ii){
    billar[ii].timestep(dt);
  }
}
//sets up gnuplot printing
void init_gnuplot(void)
{
  std::cout << "unset key" <<std::endl;
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange [0:1]" << std::endl;
  //std::cout << "set xrange [-0.5:" << lx+0.5 << "]" << std::endl;//for square table
  //std::cout << "set yrange [-0.5:" << ly+0.5 << "]" << std::endl;
  std::cout << "set xrange [" <<-0.5 - R << ":" << R +0.5 <<  "]" << std::endl;// for stadium
  std::cout << "set yrange [" <<-0.5-alpha-R <<":"<< alpha+R+0.5 << "]" << std::endl;
  print_table2();
  
}
void print_table1(void)
{
  std::cout << "plot " << lx << ","<< ly <<"* t , " ;
  std::cout << lx << "*t ," << ly << "," ;
  std::cout << 0 << "," <<  ly << "*t , " ;
  std::cout << lx << "*t ," << 0 << std::endl;
}

void print_table2(void)
{
  std::cout << "plot " << R << "*cos(pi*t),"<< alpha << "+" << R << "*sin(pi*t)," ;
  std::cout << R << "*cos(pi*t),-("<< alpha << "+" << R << "*sin(pi*t))," ;
  std::cout << R << "," <<  alpha << "*t , " ;
  std::cout << R << "," << - alpha << "*t, " ;
  std::cout << -R << "," <<  alpha << "*t , " ;
  std::cout << -R << "," << - alpha << "*t " <<std::endl ;
}
void print_gnuplot(body billar[])
{
  std::cout <<"replot "; 
  for (int ii = 0; ii < N; ++ii){
    std::cout << billar[ii].r(0) << " + " << rad << "*cos(2*pi*t),"
	      << billar[ii].r(1) << " + " << rad << "*sin(2*pi*t)  lt " << ii+1 << ",";
  }
  std::cout << " 0, 0";
  std::cout << std::endl;
}
