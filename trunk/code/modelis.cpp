#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <cstdlib>


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
         values.push_back( sqrt(h) * rnd(-1,1) );
   values.push_back( rnd(-1,1) );
   return values;

}

//////////////////////////////////////////////
double distSup(vector<double>& a, vector<double>& b){ // equal sizes

     int sz = a.size();
     double mx = 0;
     for (int i=0; i<sz;i++) {
        double diff = abs(a[i] - b[i]); 
        if (mx < diff)  mx = diff;
     }  

     return mx;
}
//////////////////////////////////////////////
double distL2(vector<double>& a, vector<double>& b){ // equal sizes

     int sz = a.size();
     double h = 1.0/sz;
     sz--;  
     double mx = 0;
     double sum = 0.0;
     for (int i=0; i<sz;i++) {
        double d1 = b[i] - a[i];
        double d2 = b[i+1] - a[i+1];
        double a1 = abs(d1);
        double b1 = abs(d2);

        if (d1 * d2 >= 0) sum += (a1+b1) * h/2;
        else  {
           double k = b1 / a1;
           double x = k * h/(1 + k);
           sum += x * b1/2 + (h - x) * a1/2 ;
        } 
     }  

     return sqrt(sum);
}


int main(void)
{

  
  setupRND( rand() );
  printf("sk\n");
  for(int i=0;i<64;i++){
    vector<double> a = getRNDC01(10000);
    for(int j=0;j<64;j++){
    
     vector<double> b = getRNDC01(10000);
     double ss = distL2(a, b);
     printf("%.5f\n", ss);
   }
  }
  destroyRND(); 
 
  return 0;
}
