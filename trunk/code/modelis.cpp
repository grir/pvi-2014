#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <cstdlib>


#include <vector>
#include <bitset>
#include <iostream>

using namespace std;

#define MAX_SIZE 2<<14

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
   for(int i=0; i < (size+1); i++){
         //double dy = sqrt(h) * rnd(-1,1);
         double dy = gsl_ran_gaussian (r, sqrt(h));
         if (dy>1) dy = 1 - dy;
         else if (dy<-1) dy = -1 + dy;
         values.push_back( dy );
   }
   //values.push_back( rnd(-1,1) );
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

/////////////////////////////////////////////////////////////////
void getLine(double x1, double y1, double x2, double y2, double& a, double& b){

    a = (y2 - y1)/(x2 - x1);
    b = y1 - a * x1;

}

//////////////////////////////////////////////////////////////////
double integral(double x1, double x2, double a, double b){

    double i1 = a * a * ( x2 * x2 * x2 -  x1 * x1 * x1) / 3;
    double i2 = a * b * x2 * x2 - a * b * x1 * x1;
    double i3 = b * b * (x2 - x1);
    return i1 + i2 + i3;

}

///////////////////////////////////////////////////////////////////
double distL2(vector<double>& a, vector<double>& b){ // equal sizes

     int sz = a.size();
     double h = 1.0/sz;
     sz--;
     double mx = 0;
     double sum = 0.0;
     for (int i=0; i<sz;i++) {
        double a1, a2, b1, b2, a12, b12;
        getLine(i * h, a[i], (i+1) * h, a[i+1], a1, b1);
        getLine(i * h, b[i], (i+1) * h, b[i+1], a2, b2);
        a12 = a2 - a1;
        b12 = b2 - b1;
           sum += integral(i*h, (i+1)*h, a12, b12);
     }

     return sqrt(sum);
}

///////////////////////////////////////////////////////////////////
double funcMin(vector<double>& a, vector<double>& b){ // equal sizes

     int sz = a.size();
     double mn = abs(a[0] - b[0]);
     for (int i = 0; i < sz; i++) {
        double mn1 = abs(a[i] - b[i]);
        if (mn1 < mn)  mn = mn1;
     }

     return mn;
}


///////////////////////////////////////////////////////////////////
double adjustDistance(vector<double>& rec, vector< vector<double> >& pts,
                      vector<int>& lig, int ligNr, double ds, double kLigLig){ // equal sizes

  //   cout << "Adjusting " << ligNr << endl;
     double deltaEnergy = 0;
     int sz = lig.size();
     double  dnow[sz], dafter[sz];
    // cout << "122" << endl;
     for(int i = 0; i < sz;i++)
        dnow[ i ] = distL2(pts[ lig[i] ], pts[ ligNr ]);

     double alpha = ds / distL2(rec, pts[ ligNr ]);
     double beta = 1 - alpha;
  //   cout << "128" << endl;

     for (int i=0; i<sz;i++) {

        pts[  ligNr  ][ i ] = beta * pts[  ligNr ][ i ] + alpha * rec[ i ];

     }

 //    cout << "136" << endl;
     for(int i=0; i<sz;i++)
        dafter[ i ] = distL2(pts[ lig[i] ], pts[ ligNr ]);

     for(int i=0; i<sz;i++){
        double delta = (dafter[i]-dnow[i]);
        deltaEnergy += kLigLig * delta * delta / 2;
     }
  //  cout << deltaEnergy << endl;
     return deltaEnergy;
}

///////////////////////////////////////////////////////////////////


/*
class Lig {
   public:
      vector< vector<double> > pts;
      int valency;
      int dim;
      double center;
      double radius;
      vector<double> d0;
      double Energy();




}

*/

////////////////////////////////////////////////////////////////////
class InteractionModel {
   public:
      vector< vector<double> > pts;
      int size;
      vector<int> lig;
      int ligSize;
      double kLigLig;         // "elastic constant"
      double eLigRec;         // univalent energy
      double R;               // interaction radius
      double ds;              // bond distance
      double (*distanceFunc)(vector<double>& a, vector<double>& b);

      void recSelection();    // randomly select ligSize ligands
      int energyByNearestReceptor();        //
   //   double energyByOptimalConf();            //
      int valency;
      InteractionModel(int size, int ligSize,
                      double kLigLig, double eLigRec,
                      double R, double ds,
                      double (*distanceFunc)(vector<double>& a, vector<double>& b)
                      );

      void printLigs();

};

///////
void InteractionModel::printLigs(){
     for(int i = 0; i < ligSize; i++)
        cout << lig[ i ] << endl;

}

///////
InteractionModel::InteractionModel(int size, int ligSize,
                      double kLigLig, double eLigRec,
                      double R, double ds,
                      double (*distanceFunc)(vector<double>& a, vector<double>& b)
                      ){

     this->size = size;
     this->ligSize = ligSize;
     this->kLigLig = kLigLig;
     this->eLigRec = eLigRec;
     this->R = R;
     this->ds = ds;
     this->distanceFunc = distanceFunc;

     this->valency = 0;
     // create points
     for(int i = 0; i < size; i++){
        vector<double> a = getRNDC01(100);
        this->pts.push_back(a);
     }
}

///////
void InteractionModel::recSelection(){
     int nums[size];
     int ligs[ligSize];

     //remove all old elements
     lig.clear();

     for(int i = 0; i < size; i++)
        nums[ i ] = i;

     gsl_ran_choose(r, ligs, ligSize, nums, size, sizeof (int));

     for(int i = 0;i < ligSize; i++)      {
        lig.push_back(ligs[i]);
     }

}

///////
int InteractionModel::energyByNearestReceptor(){

     double dL[ligSize][ligSize];
     bitset<MAX_SIZE>  unusedPts;
     unusedPts.set();                     //at start all pts are unused

     double energy = 0.0;

     for(int i = 0; i < ligSize; i++){
        unusedPts.set(lig[ i ],0);
        for(int j = 0; j < ligSize; j++) {
            dL[i][j] = distL2(pts[lig[i]],pts[lig[j]]);
            dL[j][i] = dL[i][j];
        }
     }

     // for each ligand/receptor, we find nearest
     int valency = 0;
     bitset<MAX_SIZE>  unusedLigs;
     unusedLigs.set();                     //at start all Ligs are unused

     for(int i = 0; i < ligSize; i++){
        if (energy > 1E-6) break;         // there is no
        double mindist = 2 * size;
        int ligNr = -1;
        int recNr = -1;
        for(int j = 0; j < ligSize; j++){
           if (unusedLigs[lig[ j ]] == 1) {
              for(int k = 0; k < size; k++) {
                  double d = distanceFunc(pts[lig[j]],pts[k]);
                  if ( (unusedPts[k] == 1) && ( d < R )){
                     if ( d <= mindist ) {
                        mindist = d;
                        ligNr =  j;
                        recNr = k;
                     //         cout << lig[ligNr] << " " << recNr << endl;

                     }
                  }
              }
           }
           else continue;
       }
       if (ligNr == -1) break;



       energy = energy - eLigRec + adjustDistance(pts[recNr],
                         pts, lig, lig[ligNr], ds, kLigLig);
     //  cout << "After adjusting" << endl;
       unusedLigs[lig[ligNr]] = 0;

       if (energy <= 1E-6) valency++;

    }
    return  valency;
}



int main(void)
{
  srand (time(NULL));
  setupRND( rand() );
//  setupRND( 11 );
  for(int i=1; i<=100; i++){
     double R =  i * 0.01;
     cout << R;
     for (int sz=32; sz <= 2<<12; sz = sz * 3 / 2){
      //cout << "size=" << sz << endl;
      //     for(int i=0; i<20; i++){
         double alpha = 0.5;
         int size = sz;
         int ligSize = 16;
         double eLigRec = 1;
         //double R = 0.5 + i * 0.01;

         double kLigLig = 1E2;

         double ds = R * alpha;

         double (*dF)(vector<double>& a, vector<double>& b);
         //dF = distL2;
         dF = funcMin;

         InteractionModel  model( size, ligSize, kLigLig, eLigRec, R, ds, dF);
      //  cout << "Created " << endl;
         model.recSelection();
      //model.printLigs();

      //cout << "Ligands selected " << endl;

         int va = model.energyByNearestReceptor();


         cout << " " << va;
      }
      cout << endl;
  }

  destroyRND();

  return 0;
}
