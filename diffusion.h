#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;

double k=1.; //conductivity
double L=1.; //Length
double Sp= 0.;
double Sc= 0.; //linearised source term S = Sc+Sp*Tp
double q = 0.; //neumann right (heat flux)
double h = 0.; //mixed right  (convective boundary condition)
double T_inf = 0.; //Temperature of convective medium

class one_D{
	public:
    double *T;
    int n;
	double delx;
    double a;
    
	one_D(int N){
		n=N;
		T=new double[n];
		delx=L/(n-1);
        a=k/delx;
	}
/********************display array**************************************/
void display(){
	cout<<"using display"<<endl;
	for(int i=0;i<n;i++){
    	cout<<fixed<<setprecision(5)<<i*delx<<"\t\t"<<T[i]<<endl;
	   }
    }
/****************************Initialisation****************************/
void initialise(double T0=0.){
	for(int i=0;i<n;i++){
    	T[i]=T0;
	}	
}	
/******************TDMA for 1D diffusion equation******************/	
void solve(){
	double aE,aW,aP,b,P[n],Q[n];
				aE=a;
			    aW=a;
				aP=aE+aW-Sp*delx;
				b=Sc*delx;
	cout<<Sp<<endl;			
	P[1]=aE/aP;
	Q[1]=(b+aW*T[0])/aP;
	int i;
if(q==0&&h==0){
	for(i=2;i<n-2;i++){
    	P[i]=aE/(aP-aW*P[i-1]);
    	Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
	}
	i=n-2;
	b=b+aE*T[i+1];	
    Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
    T[i]=Q[i];
    for(int i=n-3;i>0;i--){
    	T[i]=P[i]*T[i+1] + Q[i];
	}
}
else if(h==0){
	for(i=2;i<n-1;i++){
    	P[i]=aE/(aP-aW*P[i-1]);
    	Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
	}
	i=n-1;
	aE=0;
	aW=a;
	aP=aE+aW-Sp*delx/2;
	b=Sc*delx/2-q;
	
    Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
    T[i]=Q[i];
    for(int i=n-2;i>0;i--){
    	T[i]=P[i]*T[i+1] + Q[i];
	}	
	
}
else{
	for(i=2;i<n-1;i++){
    	P[i]=aE/(aP-aW*P[i-1]);
    	Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
	}
	i=n-1;
	aE=0;
	aW=a;
	aP=aE+aW+h-Sp*delx/2;
	b=h*T_inf+Sc*delx/2;
	
    Q[i]=(b+aW*Q[i-1])/(aP-aW*P[i-1]);
    T[i]=Q[i];
    for(int i=n-2;i>0;i--){
    	T[i]=P[i]*T[i+1] + Q[i];
	}	
	
}
    cout<<"using TDMA\n";
    display();

}
/********************flux calculation****************************************/
void flux(){
	double q[n];
	int i=0;
	q[i]= - k *(-3*T[i]+4*T[i+1]-T[i+2])/(2*delx);
	
	for(i=1;i<n-1;i++){
		q[i]= -k * (T[i+1]-T[i-1])/(2*delx);
	}
	
	i=n - 1;
	q[i]= - k *(3*T[i]-4*T[i-1]+T[i-2])/(2*delx);
	
	for(i=0;i<n;i++){
    	cout<<fixed<<setprecision(5)<<i*delx<<"\t\t"<<q[i]<<endl;
	}
}		
};

class two_D{
	public:
		
    double **T;
    int n,m;
	double delx,dely;
    double a;
    
	two_D(int N,int M){
		n=N;   // no. of node points in x direction
		m=M;   // no. of node points in y direction
	    T = new double * [m];
     	for(int i=0;i<m;i++){
	       T[i]=new double[n];
 	    }
	}

    /****************************Initialisation****************************/
void initialise(double T0=0.){
	for(int j=0;j<m;j++){                // j's are rows
 	     	for(int i=0;i<n;i++){        //i's s are columns
 		    	T[j][i]=T0;
		   }
	   }	
}
/********************************display*****************************/	
void display(){
	int i,j;
	for(int j=0;j<m;j++){
 		for(int i=0;i<n;i++){
 			cout<<fixed<<setprecision(5)<<T[j][i]<<" ";
		 }
		 cout<<endl;
	 }
}

/****************Gauss siedal solution for all dirichlet bc**********************/
void solve(){
	double aE,aW,aP,aN,aS,b,tempo,error;
			delx=L/(n-1);
			dely=L/(m-1);
			aE=k*dely/delx;
			aW=aE;
			aN=k*delx/dely;
			aS=aN;
			aP=aE+aW+aN+aS-Sp*delx*dely;
			b=Sc*delx*dely;
	int count=0,i,j;
	do{
		count++;
		error=0.;
		for (j=1;j<m-1;j++){
	    	for (i=1;i<n-1;i++){
	    		tempo=T[j][i];
	    		T[j][i]=(aE*T[j][i+1]+aW*T[j][i-1]+aN*T[j+1][i]+aS*T[j-1][i]+b)/aP;
	    		error+=pow((T[j][i]-tempo),2);
			}
	}	
		error=sqrt(error);
		cout<<count<<"\t"<<error<<endl;
	}while(error>0.000001);
    cout<<"using GS\n";
    display();
}
};



