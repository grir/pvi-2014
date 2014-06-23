#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <cmath>

#include <vector>
#include <iostream>

using namespace std;

/////////////////////////////////////
// Global:
gsl_rng * r;
//
/////////////////////////////////////////////
void setupRND(int seed){
   const gsl_rng_type * T;
   
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);
   gsl_rng_set (r, seed); 

}
/////////////////////////////////////////////
double rnd(){

   return gsl_rng_uniform (r); 

}
/////////////////////////////////////////////
double rnd(double a, double b){

   return gsl_rng_uniform (r) * (b - a) + a; 

}

//////////////////////////////////////////////
void destroyRND(){

   gsl_rng_free (r);

}
//////////////////////////////////////////////
vector<double> getRNDC01(unsigned int size){ // size >= 1, i.e. number of inner pts

   double h = 1.0/(size+1);
   vector<double> values;
   values.push_back( rnd(-1,1) );
   for(int i=0; i < size; i++)
         values.push_back( rnd(-1,1) );
   values.push_back( rnd(-1,1) );
   return values;

}

//////////////////////////////////////////////
double dist(vector<double>& a, vector<double>& b){ // equal sizes

     int sz = a.size();
     double mx = 0;
     for (int i=0; i<sz;i++) {
        double diff = abs(a[i] - b[i]); 
        if (mx < diff)  mx = diff;
     }  

     return mx;
}


int main(void)
{

  int i, n = 10;

  setupRND(124);
  for(i=0;i<10000;i++){
     vector<double> a = getRNDC01(1);
     vector<double> b = getRNDC01(1);
     double ss = dist(a, b);
     printf("%.5f\n", ss);
  }

  destroyRND(); 
 
  return 0;
}
