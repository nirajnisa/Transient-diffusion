#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	int m,n;
	double height,width;
	printf("ENTER THE HEIGHT OF RECTANGLE PLATE:\n");
	scanf("%lf",&height);
	printf("ENTER THE WIDTH OF RECTANGLE PLATE:\n");
	scanf("%lf",&width);
	printf("\n");
	printf("Number of grid along the x(width) direction:\n");
	scanf("%d",&m);
	printf("Number of grid along the y(height) direction:\n");
	scanf("%d",&n);

	int np=m*n; // total number of points
	
	double dx = 1.0/(double)m;
	double dy = 1.0/(double)n;
	
	double T[m+1][n+1],T_prev[m+1][n+1],temp;
	double bb, residual;
	
	int i,j;
	double gamma=1.0; //diffusion coefficient
	double dt=1.0e-5; //time step;
	double t=0.0;
	double error=1.0;
	double errortime=0.0;
	
	double volume= dx*dy;
	double areae = dy, areaw = dy, distancex = dx;
	double arean = dx, areas = dx, distancey = dy;
	
	double as,aw,ap,ae,an,ap0,anb;
	ae=(gamma*areae)/distancex;
	aw=(gamma*areaw)/distancex;
	an=(gamma*arean)/distancey;
	as=(gamma*areas)/distancey;
	
	ap0=volume/dt;
	anb=ae+aw+an+as;
	ap=ap0+anb;
	
	//Initialization and boundary conditions
	
	for(j=0;j<=n;j++)
	{
		T[0][j]=0.0;
		T[m][j]=0.0;
			
	}
	
	for(i=0;i<=m;i++)
	{
		T[i][0]=0.0;
		T[i][n]=1.0;
	}
	
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			T[i][j]=0.0;
		}
	}
	int iterations=0,timestep=0;
	
	//Implicit
	//outer Time Loop
	
	do
	{
		//copy present to previous
		for(i=1;i<m;i++)
		{
			for(j=1;j<n;j++)
			{
				T_prev[i][j]=T[i][j];
			}
		}
	// Inner loop
	iterations =0;
	while(error>0.00001)
	{
		error=0.0;
		for(i=1;i<m;i++)
		{
			for(j=1;j<n;j++)
			{
				residual=ap0*T_prev[i][j] - ap*T[i][j] + ae*T[i+1][j] + aw*T[i-1][j] + an*T[i][j+1] + as*T[i][j-1];
				
				temp=T[i][j];
				T[i][j]=T[i][j]+residual/ap;
				error=error+pow((T[i][j]-temp),2.0);
			}
		}	
		error=sqrt(error/np);
		printf("Iteration-%d\tError=%e\n",iterations,error);
		iterations++;
		
	}
	error=1.0;
	 errortime=0.0;
	 for(i=1;i<m;i++)
	 {
	 	for(j=1;j<n;j++)
	 	{
	 	  errortime=errortime+pow((T[i][j]-T_prev[i][j]),2.0);	
		}
	 }
	errortime=sqrt(errortime/np);
	t=t+dt;
	timestep++;
	printf("Time=%lf\ttimestep-%d\tError=%e\n",t,timestep,errortime);
	 
   }while(errortime>1.0e-8);
   
   FILE *fp;
   fp=fopen("temperature_Implicit.dat","w");
   fp=fopen("temperature_Implicit.txt","w");
   fprintf(fp,"I=%d\t\t  J=%d\t\t ZONE\n",m+1,n+1);
   for(j=0;j<=n;j++)
   {
   	for(i=0;i<=m;i++)
    	{
   			fprintf(fp,"%lf\t%lf\t%lf\n",i*dx,j*dy,T[i][j]);
     	}
     	fprintf(fp,"\n");
   }
   fclose(fp);
   return 0;
	 	
	
}
