#Fermi Ulma Pasta problem
#Amit Kumar, 2010
#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;

int main(){

ofstream ofl("alpha.dat",ios::out);
ofstream ofl1("apost.dat",ios::out);
ofstream ofl2("aend.dat",ios::out);

int N = 32;
long double M = 32;
long double beta = 0.03;
long double tmax = 50;
long double h = 0.001;
long double x[N+2],xold[N+2];
long double omega[N];
long double Pi = 3.14159;
long double v[N],force1[N],force2[N],sum1,sum2;
long double xnew[N],vnew[N];
long double omega2[9],energy[9],omegak2[9],Adot[9],A[9];
long double var1,var2,var3;
long double amp;
long double eng;

// boundary conditions //
x[0] = 0;
x[N+1] = 0;

// initialization //

for(int i=1;i<=N;i++){
amp = 10;
x[i] = sqrt(2/(M))*amp*sin(Pi*i/(M));
xold[i] = x[i];
v[i] = 0;
cout <<  x[i] << endl;
}
cout << sqrt(2/(M))  << "     " << pow(2*sin(Pi/(2*(M))),2)  << "          "  <<
 4*sin(Pi/(2*(M)))*sin(Pi/(2*(M))) << endl;

for(int k =1;k<=4;k++){
         sum1 = 0;
         sum2 = 0;
     for(int m=1;m<=N;m++){
         sum1 = sum1 + x[m]*sin(m*k*Pi/(M));
         sum2 = sum2 + v[m]*sin(m*k*Pi/(M));
     }
     A[k] = sum1/sqrt((M)/2);
     Adot[k] = sum2/sqrt((M)/2);
     omegak2[k]= 2*sin(Pi*k/(2*(M)));
     var1 = A[k];
     var2 = Adot[k];
     var3 = omega2[k];
     energy[k] = pow(Adot[k],2) + pow(A[k]*omegak2[k],2);
     cout << var1 << "     " << var2 << "    " << omegak2[k] << "     " <<  var3 <<
"      "  << energy[k] << "    " << eng <<  endl;
     }
    cout << energy[1]  << "      "  <<  energy[2] << "       "  << energy[3] << "   " << energy[4] << endl;



// force integration beta model//

for(long double t=0;t<tmax;t=t+h){                

// finding energy //

for(int k =1;k<=8;k++){
         sum1 = 0;
         sum2 = 0;
     for(int m=1;m<=N;m++){
         sum1 = sum1 + x[m]*sin(m*k*Pi/(M));
         sum2 = sum2 + v[m]*sin(m*k*Pi/(M));
     }
     A[k] = sum1/sqrt((M)/2);
     Adot[k] = sum2/sqrt((M)/2);
     omegak2[k]= 2*sin(Pi*k/(2*(M)));
     var1 = A[k];
     var2 = Adot[k];
     var3 = omega2[k];
     energy[k] = pow(Adot[k],2) + pow(A[k]*omegak2[k],2);   
    }
     ofl << t <<  "       "  << energy[1]  << "      "  <<  energy[2] << "       " << energy[3] << "        "   << energy[4] << "      " << energy[5] << "     "<< energy[6] << "     " << energy[7] << "     " << energy[8] << endl;
  
// updating position and velocty //

     for(int i=1;i<=N;i++){
     force1[i] = (x[i+1] - 2*x[i] + x[i-1]) + beta*(pow((x[i+1] - x[i]),2) - pow((x[i] - x[i-1]),2));
     ofl2 << force1[i] << endl;
     xnew[i] = 2*x[i] - xold[i]  + h*h*force1[i];
     ofl1 << t << "     " << x[i] << "    " << v[i] << endl;
     v[i] = (xnew[i] - xold[i])/(2*h);
     xold[i] = x[i];
     x[i] = xnew[i];
     }
}
}   
   
