/**********************************************************************
UDF to calculate temperature field function and store in
user-defined memory. Also print min, max, avg temperatures.
***********************************************************************/
#include "udf.h"
DEFINE_ON_DEMAND(entropy)
{
	Domain *d; /* declare domain pointer since it is not passed as an argument to the DEFINE macro  */
	real dtx,dty,dtz,k,kf,mu,ux,uy,vx,vy,dT,T0,R,Ns;
	Thread *t;
	cell_t c;
	d = Get_Domain(1);     /* Get the domain using ANSYS FLUENT utility */
	/* Loop over all cell threads in the domain */
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
			T0 = 310;
			dT = 10;
           R = 0.0018617;
			kf = 0.628;
			dtx = C_T_G(c,t)[0];
			dty = C_T_G(c,t)[1];
			vx = C_V_G(c,t)[0];
			vy = C_V_G(c,t)[1];
			ux = C_U_G(c,t)[0];
			uy = C_U_G(c,t)[1];
			k = C_K_L(c,t);
			mu = C_MU_L(c,t);
			Ns = (kf*dT*dT)/(T0*T0*R*R);

			C_UDMI(c,t,0) = (dtx*dtx+dty*dty)*k/(T0*T0)/Ns;
			C_UDMI(c,t,1) = (mu/T0)*(2*ux*ux+2*vy*vy+(uy+vx)*(uy+vx))/Ns;

           }
		end_c_loop(c,t)
	}
}