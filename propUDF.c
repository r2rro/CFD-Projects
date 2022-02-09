#include "udf.h"
#include "math.h" 

DEFINE_PROPERTY(cell_conductivity,cell,thread)
{
   real k,kbrown,kstatic,Kapa,Beta,f,dp,Cpf;
	real temp = C_T(cell,thread);
	real phi = 1.0-C_VOF(cell, thread);
	real rhof= 993;
	real rhop= 3970;
	real kp= 36;
	real kbf= 0.628;

	Kapa=1.38064852*pow(10,-23);   Cpf= 4178;  dp= 25*pow(10,-9);

   Beta=0.0017*pow(100*phi,-0.0841);
	f=(-0.8467*phi+0.0753)*temp+(237.67*phi-21.998);

	kstatic=kbf*(kp+2*kbf+2*(kp-kbf)*phi)/(kp+2*kbf-(kp-kbf)*phi);
	kbrown=5*pow(10,4)*Beta*phi*rhof*Cpf*sqrt(Kapa*temp/(rhop*dp))*f;

	k=kstatic+kbrown;

	return k;
} 

DEFINE_PROPERTY(cell_viscosity, cell, thread)
{
   real mu_nf; /*Effective viscosity of nanofluid*/
   real phi = 1.0-C_VOF(cell, thread);
   real df= 0.385*pow(10,-9);
   real dp= 25*pow(10,-9);
   real visc_dyn_min=0.0001; /*Define the minimum limit for the value of viscosity*/
   real visc_dyn_max=10000; /*Define the maximum limit for the value of viscosity*/
   real n= 1; /*power law index*/
   real C= 0.0015031; /*consistency index*/
   real T0=0; /*Define the reference temperature T0*/
   real temp=C_T(cell,thread);
   real mu_app; /*Define the variable for the apparent viscosity of fluid*/
   real mu_temp; /*Define a temporary variable to stock the apparent viscosity*/ 

   mu_temp=C*pow(C_STRAIN_RATE_MAG(cell,thread),n-1)*exp(T0/temp);

   if (N_ITER<1)
   {
   mu_app=C; /*Initialize mu_app to the value of k*/
   }

   else if ((N_ITER>=1)&& (mu_temp>visc_dyn_min)&& (mu_temp<visc_dyn_max))
   {
   mu_app=mu_temp;
   }

   else if ((N_ITER>=1)&&(mu_temp>=visc_dyn_max))
   {
   mu_app=visc_dyn_max;
   }

   else if ((N_ITER>=1)&&(mu_temp<=visc_dyn_min))
   {
   mu_app=visc_dyn_min;
   }

   mu_nf=mu_app/(1.0-34.87*pow((dp/df),-0.3)*pow(phi,1.03));
   return mu_nf;

   }