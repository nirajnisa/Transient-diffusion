#include "diffusion.h"

int main(){
	Sc=100.;
    two_D t(10,10);
    t.initialise(0.);
    for(int i=0;i<t.n;i++)
    	t.T[0][i]=10.;
    for(int j=0;j<t.m;j++)
    	t.T[j][0]=10.;
    t.solve();
    FILE *fp =fopen("temp.txt","w");
    delete []t.T;

return 0;
}  

