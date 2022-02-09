//Arthur Rostami
//
//This Program solves 2-dimensional Navier Stokes equation using Implicit Euler as time advance schemes
//and Second order centered scheme for space discretization

#include<stdio.h>
#include<math.h>
#include<ctime>
#define MAXSIZE (200+2)

const int I=30; // number of Cells in x direction
const int J=60; // number of Cells in y direction
const double dt=0.05;//Time step
const int tmax=200; //maximum number of marching in time
const double Re=250; //Reynolds number
const double p0=1; //Pressure coefficient
const double u0=1; //x velocity coefficient
const double v0=1; //y velocity coefficien
const double ORF=1.0; //Over relaxation factor
const double beta=1.0; 
const double h=2.0; //Height of the cavity
const double w=1.0; //Width of the cavity
const double U_top=1;

double l1_p, l1_u, l1_v, l2_p, l2_u, l2_v, linf_p, linf_u, linf_v;
FILE *f1 = fopen("part 4 h1 convergence history.txt", "w"); // opening the file that we like to save the data in
FILE *f2 = fopen("part 4 h1 Contour.txt", "w"); // opening the file that we like to save the data in

static void SpewMatrix(double Source[3][3])
{
	printf("%10.6f %10.6f %10.6f\n", Source[0][0], Source[1][0], Source[2][0]);
	printf("%10.6f %10.6f %10.6f\n", Source[0][1], Source[1][1], Source[2][1]);
	printf("%10.6f %10.6f %10.6f\n", Source[0][2], Source[1][2], Source[2][2]);
}

static void SpewVector(double Source[3])
{
	printf("%10.6f %10.6f %10.6f\n", Source[0], Source[1], Source[2]);
}

static inline void CopyVec(const double Source[3], double Target[3])
{
	Target[0] = Source[0];
	Target[1] = Source[1];
	Target[2] = Source[2];
}

static inline void Copy3x3(double Source[3][3], double Target[3][3])
{
	Target[0][0] = Source[0][0];
	Target[0][1] = Source[0][1];
	Target[0][2] = Source[0][2];
	Target[1][0] = Source[1][0];Target[1][1] = Source[1][1];
	Target[1][2] = Source[1][2];
	Target[2][0] = Source[2][0];
	Target[2][1] = Source[2][1];
	Target[2][2] = Source[2][2];
}

static inline void Mult3x3(double A[3][3], double B[3][3], double C[3][3])
{
	C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
	C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
	C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
	C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
	C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
	C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
	C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
	C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
	C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
}

static inline void MultVec(double A[3][3], const double Vec[3], double Result[3])
{
	Result[0] = A[0][0]*Vec[0] + A[0][1]*Vec[1] + A[0][2]*Vec[2];
	Result[1] = A[1][0]*Vec[0] + A[1][1]*Vec[1] + A[1][2]*Vec[2];
	Result[2] = A[2][0]*Vec[0] + A[2][1]*Vec[1] + A[2][2]*Vec[2];
}

static inline void Add3x3(double A[3][3], double B[3][3], const double Factor, double C[3][3])
{
	C[0][0] = A[0][0] + Factor * B[0][0];
	C[0][1] = A[0][1] + Factor * B[0][1];
	C[0][2] = A[0][2] + Factor * B[0][2];
	C[1][0] = A[1][0] + Factor * B[1][0];
	C[1][1] = A[1][1] + Factor * B[1][1];
	C[1][2] = A[1][2] + Factor * B[1][2];
	C[2][0] = A[2][0] + Factor * B[2][0];
	C[2][1] = A[2][1] + Factor * B[2][1];
	C[2][2] = A[2][2] + Factor * B[2][2];
}

static inline void AddVec(const double A[3], const double B[3],const double Factor, double C[3])
{
	C[0] = A[0] + Factor * B[0];
	C[1] = A[1] + Factor * B[1];
	C[2] = A[2] + Factor * B[2];
}

static inline void Invert3x3(double Block[3][3], double Inverse[3][3])
{
	double DetInv = 1. / (+ Block[0][0]*Block[1][1]*Block[2][2]
	+ Block[0][1]*Block[1][2]*Block[2][0]
	+ Block[0][2]*Block[1][0]*Block[2][1]
	- Block[0][2]*Block[1][1]*Block[2][0]
	- Block[0][1]*Block[1][0]*Block[2][2]
	- Block[0][0]*Block[1][2]*Block[2][1]);
	/* Expand by minors to compute the inverse */
	Inverse[0][0] = + DetInv * (Block[1][1]*Block[2][2] -
	Block[2][1]*Block[1][2]);
	Inverse[1][0] = - DetInv * (Block[1][0]*Block[2][2] -
	Block[2][0]*Block[1][2]);
	Inverse[2][0] = + DetInv * (Block[1][0]*Block[2][1] -
	Block[2][0]*Block[1][1]);
	Inverse[0][1] = - DetInv * (Block[0][1]*Block[2][2] -
	Block[2][1]*Block[0][2]);
	Inverse[1][1] = + DetInv * (Block[0][0]*Block[2][2] -
	Block[2][0]*Block[0][2]);
	Inverse[2][1] = - DetInv * (Block[0][0]*Block[2][1] -
	Block[2][0]*Block[0][1]);
	Inverse[0][2] = + DetInv * (Block[0][1]*Block[1][2] -
	Block[1][1]*Block[0][2]);
	Inverse[1][2] = - DetInv * (Block[0][0]*Block[1][2] -
	Block[1][0]*Block[0][2]);
	Inverse[2][2] = + DetInv * (Block[0][0]*Block[1][1] -
	Block[1][0]*Block[0][1]);
}

void SolveBlockTri(double LHS[MAXSIZE][3][3][3], double RHS[MAXSIZE][3], int iNRows)
{
	int j;
	double Inv[3][3];
	for (j = 0; j < iNRows-1; j++) 
	{
		/* Compute the inverse of the main block diagonal. */
		Invert3x3(LHS[j][1], Inv);
		/* Scale the right-most block diagonal by the inverse. */
		{
		double Temp[3][3];
		Mult3x3(Inv, LHS[j][2], Temp);
		Copy3x3(Temp, LHS[j][2]);
		}
		/* Scale the right-hand side by the inverse. */
		{
		double Temp[3];
		MultVec(Inv, RHS[j], Temp);
		CopyVec(Temp, RHS[j]);
		}
		/* Left-multiply the jth row by the sub-diagonal on the j+1st row
		and subtract from the j+1st row. This involves the
		super-diagonal term and the RHS of the jth row. */
		{
		/* First the LHS manipulation */
		#define A LHS[j+1][0]
		#define B LHS[j+1][1]
		#define C LHS[j ][2]
		double Temp[3][3], Temp2[3][3];
		double TVec[3], TVec2[3];
		Mult3x3(A, C, Temp);
		Add3x3(B, Temp, -1., Temp2);
		Copy3x3(Temp2, B);
		/* Now the RHS manipulation */
		MultVec(A, RHS[j], TVec);
		AddVec(RHS[j+1], TVec, -1., TVec2);
		CopyVec(TVec2, RHS[j+1]);
		#undef A
		#undef B
		#undef C
		}
	} /* Done with forward elimination loop */
	/* Compute the inverse of the last main block diagonal. */
	j = iNRows-1;
	Invert3x3(LHS[j][1], Inv);
	/* Scale the right-hand side by the inverse. */
	{
	double Temp[3];
	MultVec(Inv, RHS[j], Temp);
	CopyVec(Temp, RHS[j]);
	}
	/* Now do the back-substitution. */
	for (j = iNRows-2; j >= 0; j--) 
	{
		/* Matrix-vector multiply and subtract. */
		#define C LHS[j][2]
		RHS[j][0] -= (C[0][0]*RHS[j+1][0] +
		C[0][1]*RHS[j+1][1] +
		C[0][2]*RHS[j+1][2]);
		RHS[j][1] -= (C[1][0]*RHS[j+1][0] +
		C[1][1]*RHS[j+1][1] +
		C[1][2]*RHS[j+1][2]);RHS[j][2] -= (C[2][0]*RHS[j+1][0] +
		C[2][1]*RHS[j+1][1] +
		C[2][2]*RHS[j+1][2]);
		#undef C
		}
}

/*This subroutine calculates a line in [3x3]x[3x1] Matrix Multiplication*/
double MtrxMlt (int i, int j, int a, int b, int m, double (*A)[J+2][3][3], double (*B)[J+2][3])
{
	return A[i][j][m][0]*B[i+a][j+b][0]+A[i][j][m][1]*B[i+a][j+b][1]+A[i][j][m][2]*B[i+a][j+b][2];
}

/*This subroutine returns the value of the exact flux integral- Part 1*/
void FI_exactval (double (*FI_exact)[J+1][3], const double dx, const double dy)
{
	double a,b,c,d;
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			FI_exact[i][j][0]=-(u0*cos(M_PI*(i-0.5)*dx)*sin(2*M_PI*(j-0.5)*dy)+v0*sin(2*M_PI*(i-0.5)*dx)*cos(M_PI*(j-0.5)*dy))*M_PI/beta;

			a=p0*M_PI*sin(M_PI*(i-0.5)*dx)*cos(M_PI*(j-0.5)*dy);
			b=-pow(u0,2)*M_PI*sin(2*M_PI*(i-0.5)*dx)*pow(sin(2*M_PI*(j-0.5)*dy),2);
			c=-u0*v0*M_PI*sin(M_PI*(i-0.5)*dx)*sin(2*M_PI*(i-0.5)*dx)*(cos(M_PI*(j-0.5)*dy)*sin(2*M_PI*(j-0.5)*dy)+2*cos(2*M_PI*(j-0.5)*dy)*sin(M_PI*(j-0.5)*dy));
			d=-u0*5*pow(M_PI,2)*sin(M_PI*(i-0.5)*dx)*sin(2*M_PI*(j-0.5)*dy)/Re;
			FI_exact[i][j][1]=a+b+c+d;
			
			a=p0*M_PI*cos(M_PI*(i-0.5)*dx)*sin(M_PI*(j-0.5)*dy);
			b=-pow(v0,2)*M_PI*pow(sin(2*M_PI*(i-0.5)*dx),2)*sin(2*M_PI*(j-0.5)*dy);
			c=-u0*v0*M_PI*sin(M_PI*(j-0.5)*dy)*sin(2*M_PI*(j-0.5)*dy)*(cos(M_PI*(i-0.5)*dx)*sin(2*M_PI*(i-0.5)*dx)+2*cos(2*M_PI*(i-0.5)*dx)*sin(M_PI*(i-0.5)*dx));
			d=-v0*5*pow(M_PI,2)*sin(2*M_PI*(i-0.5)*dx)*sin(M_PI*(j-0.5)*dy)/Re;
			FI_exact[i][j][2]=a+b+c+d;
		}
	}
}
	
/*This Subroutine calculates the first initial condition of the project - Except for ghost cell*/
void initial(double (*U)[J+2][3], const double dx, const double dy)
{
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
		U[i][j][0]=p0*cos(M_PI*(i-0.5)*dx)*cos(M_PI*(j-0.5)*dy);
		U[i][j][1]=u0*sin(M_PI*(i-0.5)*dx)*sin(2*M_PI*(j-0.5)*dy);
		U[i][j][2]=v0*sin(2*M_PI*(i-0.5)*dx)*sin(M_PI*(j-0.5)*dy);
		}
	}
}

/*This Subroutine calculates the ghostcell values for T*/
void U_GhostCell (double (*U)[J+2][3], const double dx, const double dy)
{
	int i, j, k;

	for (i=1; i<=I; i++)
	{
		U[i][0][0]=U[i][1][0]; //BL approx
		U[i][0][1]=-U[i][1][1]; //no-slip
		U[i][0][2]=-U[i][1][2]; //no-slip
	}
	
	for (i=1; i<=I; i++)
	{
		U[i][J+1][0]=U[i][J][0]; //BL approx
		U[i][J+1][1]=2*U_top-U[i][J][1]; //no-slip
		U[i][J+1][2]=-U[i][J][2]; //no-slip
	}
	
	for (j=1; j<=J; j++)
	{
		U[0][j][0]=U[1][j][0]; //BL approx
		U[0][j][1]=-U[1][j][1]; //no-slip
		U[0][j][2]=-U[1][j][2]; //no-slip
	}
	
	for (j=1; j<=J; j++)
	{
		U[I+1][j][0]=U[I][j][0]; //BL approx
		U[I+1][j][1]=-U[I][j][1]; //no-slip
		U[I+1][j][2]=-U[I][j][2]; //no-slip
	}
}

/*This subroutine calculates the flux integral*/
void FInt (double (*U)[J+2][3], double (*FI)[J+2][3], double (*E)[J+2][3], const double dx, const double dy, const double dt)
{
	double F1[3], F2[3], G1[3], G2[3], FI_exact[I+1][J+1][3];
	//double src1[I+2][J+2];
	
	FI_exactval(FI_exact, dx, dy);
	
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			F1[0]=(U[i+1][j][1]+U[i][j][1])/(2*beta);
			F1[1]=pow((U[i+1][j][1]+U[i][j][1])/2,2)+(U[i+1][j][0]+U[i][j][0])/2-(U[i+1][j][1]-U[i][j][1])/Re/dx;
			F1[2]=(U[i+1][j][1]+U[i][j][1])*(U[i+1][j][2]+U[i][j][2])/4-(U[i+1][j][2]-U[i][j][2])/Re/dx;
			
			F2[0]=(U[i][j][1]+U[i-1][j][1])/(2*beta);
			F2[1]=pow((U[i][j][1]+U[i-1][j][1])/2,2)+(U[i][j][0]+U[i-1][j][0])/2-(U[i][j][1]-U[i-1][j][1])/Re/dx;
			F2[2]=(U[i][j][1]+U[i-1][j][1])*(U[i][j][2]+U[i-1][j][2])/4-(U[i][j][2]-U[i-1][j][2])/Re/dx;
			
			G1[0]=(U[i][j+1][2]+U[i][j][2])/(2*beta);
			G1[1]=(U[i][j+1][1]+U[i][j][1])*(U[i][j+1][2]+U[i][j][2])/4-(U[i][j+1][1]-U[i][j][1])/Re/dy;
			G1[2]=pow((U[i][j+1][2]+U[i][j][2])/2,2)+(U[i][j+1][0]+U[i][j][0])/2-(U[i][j+1][2]-U[i][j][2])/Re/dy;
			
			G2[0]=(U[i][j][2]+U[i][j-1][2])/(2*beta);
			G2[1]=(U[i][j][1]+U[i][j-1][1])*(U[i][j][2]+U[i][j-1][2])/4-(U[i][j][1]-U[i][j-1][1])/Re/dy;
			G2[2]=pow((U[i][j][2]+U[i][j-1][2])/2,2)+(U[i][j][0]+U[i][j-1][0])/2-(U[i][j][2]-U[i][j-1][2])/Re/dy;
			
			for (int k=0; k<3; k++)
			{
				FI[i][j][k]=-((F1[k]-F2[k])/dx+(G1[k]-G2[k])/dy);
			//	E[i][j][k]=sqrt(pow(FI[i][j][k]-FI_exact[i][j][k],2));
			}
		}
	}
}

/*This subroutine is for calculating jacobians and the coefficients of the equation*/
void Jacobian (double (*U)[J+2][3], double (*A_x)[J+2][3][3], double (*B_x)[J+2][3][3], double (*C_x)[J+2][3][3], double (*A_y)[J+2][3][3], double (*B_y)[J+2][3][3], double (*C_y)[J+2][3][3], const double dx, const double dy)
{
	int i,j,k;
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			A_x[i][j][0][0]=0;
			A_x[i][j][0][1]=-1/(2*beta)/dx;
			A_x[i][j][0][2]=0;
			
			A_x[i][j][1][0]=-0.5;
			A_x[i][j][1][1]=-((U[i][j][1]+U[i-1][j][1])/2+1/(Re*dx))/dx;
			A_x[i][j][1][2]=0;
			
			A_x[i][j][2][0]=0;
			A_x[i][j][2][1]=-(U[i][j][2]+U[i-1][j][2])/4/dx;
			A_x[i][j][2][2]=-((U[i][j][1]+U[i-1][j][1])/4+1/(Re*dx))/dx;
		}
	}
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			B_x[i][j][0][0]=0-0;
			B_x[i][j][0][1]=(1/(2*beta)-1/(2*beta))/dx;
			B_x[i][j][0][2]=0-0;
			
			B_x[i][j][1][0]=(0.5-0.5)/dx;
			B_x[i][j][1][1]=((U[i][j][1]+U[i+1][j][1])/2+1/(Re*dx)-((U[i][j][1]+U[i-1][j][1])/2-1/(Re*dx)))/dx;
			B_x[i][j][1][2]=0-0;
			
			B_x[i][j][2][0]=0-0;
			B_x[i][j][2][1]=((U[i][j][2]+U[i+1][j][2])/4-(U[i][j][2]+U[i-1][j][2])/4)/dx;
			B_x[i][j][2][2]=((U[i][j][1]+U[i+1][j][1])/4+1/(Re*dx)-((U[i][j][1]+U[i-1][j][1])/4-1/(Re*dx)))/dx;
		}
	}
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			C_x[i][j][0][0]=0;
			C_x[i][j][0][1]=1/(2*beta)/dx;
			C_x[i][j][0][2]=0;
			
			C_x[i][j][1][0]=0.5/dx;
			C_x[i][j][1][1]=((U[i][j][1]+U[i+1][j][1])/2-1/(Re*dx))/dx;
			C_x[i][j][1][2]=0;
			
			C_x[i][j][2][0]=0;
			C_x[i][j][2][1]=(U[i][j][2]+U[i+1][j][2])/4/dx;
			C_x[i][j][2][2]=((U[i][j][1]+U[i+1][j][1])/4-1/(Re*dx))/dx;
		}
	}
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			A_y[i][j][0][0]=0;
			A_y[i][j][0][1]=0;
			A_y[i][j][0][2]=-1/(2*beta)/dy;
			
			A_y[i][j][1][0]=0;
			A_y[i][j][1][1]=-((U[i][j][2]+U[i][j-1][2])/4+1/(Re*dy))/dy;
			A_y[i][j][1][2]=-(U[i][j][1]+U[i][j-1][1])/4/dy;
			
			A_y[i][j][2][0]=-0.5/dy;
			A_y[i][j][2][1]=0;
			A_y[i][j][2][2]=-((U[i][j][2]+U[i][j-1][2])/2+1/(Re*dy))/dy;
		}
	}
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			B_y[i][j][0][0]=0-0;
			B_y[i][j][0][1]=0-0;
			B_y[i][j][0][2]=(1/(2*beta)-1/(2*beta))/dy;
			
			B_y[i][j][1][0]=0-0;
			B_y[i][j][1][1]=((U[i][j][2]+U[i][j+1][2])/4+1/(Re*dy)-((U[i][j][2]+U[i][j-1][2])/4-1/(Re*dy)))/dy;
			B_y[i][j][1][2]=((U[i][j][1]+U[i][j+1][1])/4-(U[i][j][1]+U[i][j-1][1])/4)/dy;
			
			B_y[i][j][2][0]=(0.5-0.5)/dy;
			B_y[i][j][2][1]=0-0;
			B_y[i][j][2][2]=((U[i][j][2]+U[i][j+1][2])/2+1/(Re*dy)-((U[i][j][2]+U[i][j-1][2])/2-1/(Re*dy)))/dy;
		}
	}
	
	for (j=1; j<=J; j++)
	{
		for (i=1; i<=I; i++)
		{
			C_y[i][j][0][0]=0;
			C_y[i][j][0][1]=0;
			C_y[i][j][0][2]=1/(2*beta)/dy;
			
			C_y[i][j][1][0]=0;
			C_y[i][j][1][1]=((U[i][j][2]+U[i][j+1][2])/4-1/(Re*dy))/dy;
			C_y[i][j][1][2]=(U[i][j][1]+U[i][j+1][1])/4/dy;
			
			C_y[i][j][2][0]=0.5/dy;
			C_y[i][j][2][1]=0;
			C_y[i][j][2][2]=((U[i][j][2]+U[i][j+1][2])/2-1/(Re*dy))/dy;
		}
	}
}

/*This subroutine uses approximately factored time advance scheme to solve incompressible Navier-Stokes equations*/
void SolveNS (double (*U)[J+2][3], double (*E)[J+2][3], int n, const double dx, const double dy, const double dt)
{
	double U_new[I+2][J+2][3], del_U[I+2][J+2][3], dU_tilde[I+2][J+2][3], FI[I+2][J+2][3], FI_new[I+2][J+2][3];
	double LHS_x[I+2][3][3][3], RHS_x[I+2][3], LHS_y[J+2][3][3][3], RHS_y[J+2][3];
	double A_x[I+2][J+2][3][3], B_x[I+2][J+2][3][3], C_x[I+2][J+2][3][3], A_y[I+2][J+2][3][3], B_y[I+2][J+2][3][3], C_y[I+2][J+2][3][3];
	
	FInt(U, FI, E, dx, dy, dt);
	Jacobian(U, A_x, B_x, C_x, A_y, B_y, C_y, dx, dy);

	//*************x-sweep*****************
	for (int j=1; j<=J; j++)
	{
		//Boundary Conditions:
		for (int l=0; l<3; l++) // l and m are Block indices
		{
			for (int m=0; m<3; m++)
			{
				LHS_x[0][0][l][m]=0;
				LHS_x[I+1][2][l][m]=0;
				
				if (l != m)
				{
					LHS_x[0][1][l][m]=0;
					LHS_x[0][2][l][m]=0;

					LHS_x[I+1][0][l][m]=0;
					LHS_x[I+1][1][l][m]=0;
				}
				
				else if (l == m)
				{
					LHS_x[0][1][l][m]=1;
					LHS_x[0][2][l][m]=1;

					LHS_x[I+1][0][l][m]=1;
					LHS_x[I+1][1][l][m]=1;
				}
			}
			RHS_x[0][l]=0;
			RHS_x[I+1][l]=0;
		}
		LHS_x[0][1][0][0]=-1;
		LHS_x[I+1][1][0][0]=-1;
		
		// Computing LHS and RHS
		for (int i=1; i<=I; i++)
		{
			for (int l=0; l<3; l++) // l and m are Block indices
			{
				for (int m=0; m<3; m++)
				{
					LHS_x[i][0][l][m]=A_x[i][j][l][m]*dt;
					LHS_x[i][2][l][m]=C_x[i][j][l][m]*dt;
				
					if (l == m)
					{
						LHS_x[i][1][l][m]=1+B_x[i][j][l][m]*dt;
					}
					
					else if (l != m)
					{
						LHS_x[i][1][l][m]=B_x[i][j][l][m]*dt;
					}
				}
				
			RHS_x[i][l]=FI[i][j][l]*dt;
			}
		}
		
		SolveBlockTri(LHS_x, RHS_x, I+2); //This gives delta U tilde
		
		for (int i=1; i<=I; i++)
		{
			for (int k=0; k<3; k++)
			{
				dU_tilde[i][j][k]=RHS_x[i][k]*ORF;
			}
		}
	}
			
		
	//*************y-sweep*****************
	for (int i=1; i<=I; i++)
	{
		//Boundary Conditions:
		for (int l=0; l<3; l++) // l and m are Block indices
		{
			for (int m=0; m<3; m++)
			{
				LHS_y[0][0][l][m]=0;
				LHS_y[J+1][2][l][m]=0;
				
				if (l != m)
				{
					LHS_y[0][1][l][m]=0;
					LHS_y[0][2][l][m]=0;

					LHS_y[J+1][0][l][m]=0;
					LHS_y[J+1][1][l][m]=0;
				}
				else
				{
					LHS_y[0][1][l][m]=1;
					LHS_y[0][2][l][m]=1;

					LHS_y[J+1][0][l][m]=1;
					LHS_y[J+1][1][l][m]=1;
				}
			}
			RHS_y[0][l]=0;
			RHS_y[J+1][l]=0;
		}
		
		LHS_y[0][1][0][0]=-1;
		LHS_y[J+1][1][0][0]=-1;
		
		// Computing LHS and RHS
		for (int j=1; j<=J; j++)
		{
			for (int l=0; l<3; l++) // l and m are Block indices
			{
				for (int m=0; m<3; m++)
				{
					LHS_y[j][0][l][m]=A_y[i][j][l][m]*dt;
					LHS_y[j][2][l][m]=C_y[i][j][l][m]*dt;
					
					if (l == m)
					{
						LHS_y[j][1][l][m]=1+B_y[i][j][l][m]*dt;
					}
					
					else if (l != m)
					{
						LHS_y[j][1][l][m]=B_y[i][j][l][m]*dt;
					}

				}
				
			RHS_y[j][l]=dU_tilde[i][j][l];
			}
		}
		
		SolveBlockTri(LHS_y, RHS_y, J+2); // This gives Delta U
		
		for (int j=1; j<=J; j++)
		{
			for (int k=0; k<3; k++)
			{
				U_new[i][j][k]=U[i][j][k]+RHS_y[j][k]*ORF;
			}
		}
	}

	for (int j=1; j<=J; j++)
		{
			for (int i=1; i<=I; i++) 
			{
				for (int k=0; k<3; k++)
				{
					E[i][j][k]=sqrt(pow(U_new[i][j][k]-U[i][j][k],2));
					U[i][j][k]=U_new[i][j][k];
				}
			}
		}

}

/*This subroutine calculates L1 norm*/
void L1Norm(double (*E)[J+2][3])
{
	double sigma[3]; //error sum
	for (int k=0;k<3;k++)
	{
		sigma[k]=0;
	}
	
	
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			for (int k=0;k<3;k++)
			{
			sigma[k]=sigma[k]+E[i][j][k]/(I*J);
			}
		}
	}
	
	l1_p=sigma[0];
	l1_u=sigma[1];
	l1_v=sigma[2];
}

/*This subroutine calculates L2 norm*/
void L2Norm(double (*E)[J+2][3])
{
	double sigma[3]; //error sum
	for (int k=0;k<3;k++)
	{
		sigma[k]=0;
	}
	
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			for (int k=0;k<3;k++)
			{
			sigma[k]=sigma[k]+E[i][j][k]*E[i][j][k]/(I*J);
			}
		}
	}

	l2_p=sqrt(sigma[0]);
	l2_u=sqrt(sigma[1]);
	l2_v=sqrt(sigma[2]);
}

/*This subroutine calculates L_inf norm*/
void L_infNorm(double (*E)[J+2][3])
{
	double max[3]; //maximum error
	for (int k=0;k<3;k++)
	{
		max[k]=0;
	}
	
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			for (int k=0;k<3;k++)
			{
				if (E[i][j][k]>max[k])
				{
					max[k]=E[i][j][k];
				}
			}
		}
	}
	
	linf_p=max[0];
	linf_u=max[1];
	linf_v=max[2];
}

/*This subroutine calculates the vorticity in the domain*/
void vorticity(double (*vrt)[J+1],double (*U)[J+2][3], const double dx, const double dy)
{
	for (int j=1; j<=J; j++)
	{
		for (int i=1; i<=I; i++)
		{
			vrt[i][j]=(U[i+1][j][2]-U[i-1][j][2])/dx-(U[i][j+1][1]-U[i][j-1][1])/dy;
		}
	}
}

/*This subroutine calculates the Circulation. The arguments are vorticity, length of the integration area, i of the vortex, j of the vortex, and mesh sizes respectively*/
double gamma(double (*vrt)[J+1], const double a, const int ii, const int jj, const double dx, const double dy)
{
	double G=0;
	//int k=(a/2)/dx;
	for (int j=jj-2;j<=jj+2;j++)
	{
		for (int i=ii-2;i<=ii+2;i++)
		{
			G=G+vrt[i][j]*dx*dy;
		}
	}
	
	return G;
}

/* main function solves and saves the Energy equation using the above-defined functions*/
int main()
{
	clock_t start, end;
	int i,j,k,n=0;
	l2_p=l2_u=l2_v=1;
	int loc[2];
	double U[I+2][J+2][3], FI[I+2][J+2][3], E[I+2][J+2][3], FI_exact[I+1][J+1][3], vrt[I+1][J+1], max=-50, min=10, cpu_time_used;
	
	const double dx=w/I; //Calculating the mesh size in x direction
	const double dy=h/J; //Calculating the mesh size in y direction

	// Initial Values:
	initial(U, dx, dy);
	U_GhostCell(U, dx, dy);
	
	fprintf(f1, "VARIABLES = \"t\",\"L2_p\",\"L2_u\",\"L2_v\"\n");
	fprintf(f2, "VARIABLES = \"x\",\"y\",\"p\",\"u\",\"v\",\"vrt\",\"v mag\"\n");
	
	start = clock(); //Start reading the system time
	
	//*************SOLVER********************
	while (l2_p>pow(10,-17) || l2_u>pow(10,-17) || l2_v>pow(10,-17))  //Marching in time(n<tmax)
	{
		U_GhostCell(U, dx, dy); //Setting ghost cell values usingGhostCell subroutines
		SolveNS(U, E, n, dx, dy, dt);
		n++;
		L2Norm(E);
		fprintf(f1, "%f %.10f %.10f %.10f %.10f\n", n*dt, l2_p, l2_u, l2_v,vrt); // reporting L2
	}
	//***************************************
	
	U_GhostCell(U, dx,dy);
	end = clock(); //End of reading the system time
	cpu_time_used=((double) (end-start))/CLOCKS_PER_SEC;

//This loop computes the location of the main vortex
	vorticity(vrt, U, dx, dy);
	for (i=0.2*w/dx+0.5; i<=0.8*w/dx+0.5; i++)
	{
	//	for (j=J/2-0.1*h/dy; j<=J-0.1*h/dy; j++)
	for (j=(h-0.6)/dy+0.5; j<=(h-0.2)/dy+0.5; j++)
		{
			if (sqrt(U[i][j][1]*U[i][j][1]+U[i][j][2]*U[i][j][2])<min)
			{
				loc[0]=i; // x-location of the vortex
				loc[1]=j; // y-location of the vortex
				min=sqrt(U[i][j][1]*U[i][j][1]+U[i][j][2]*U[i][j][2]);
			}
		}
	}

	for (i=1; i<=I; i++)
	{
		for (j=1; j<=J; j++)
		{
			fprintf(f2, "%f %f %.10f %.10f %.10f %.10f %.10f\n", (i-0.5)*dx, (j-0.5)*dy, U[i][j][0], U[i][j][1], U[i][j][2], vrt[i][j], sqrt(U[i][j][1]*U[i][j][1]+U[i][j][2]*U[i][j][2]));
		}
	}

//For grid indipendancy test
printf("%f %f %f", (U[I/2][J][0]+U[I/2+1][J][0]+U[I][J/2][0]+U[I][J/2+1][0])/4, (U[I/2][J][1]+U[I/2+1][J][1]+U[I][J/2][1]+U[I][J/2+1][1])/4,(U[I/2][J][2]+U[I/2+1][J][2]+U[I][J/2][2]+U[I][J/2+1][2])/4);


	/*for (j=1; j<=J; j++)
	{
	for (i=1; i<=I; i++)
	{
			printf("%f,",  vrt[i][j]);//, (j-0.5)*dy,U[i][j][2]);
		}
		printf("\n");
	}
			
	printf("\n\n");
	
	for (j=1; j<=J; j++)
	{
		printf("u(:,%d)=[", j);
	for (i=1; i<=I; i++)
	{
			printf("%f,",  U[i][j][1]);//, (j-0.5)*dy,U[i][j][2]);
		}
		printf("];\n");
	}

	for (j=1; j<=J; j++)
	{
		printf("v(:,%d)=[", j);
	for (i=1; i<=I; i++)
	{
			printf("%f,",  U[i][j][2]);//, (j-0.5)*dy,U[i][j][2]);
		}
		printf("];\n");
	}*/

	L1Norm(E);
	L2Norm(E);
	L_infNorm(E);
	
	printf("\np L1 norm is: %.16f", l1_p); // reporting L1
	printf("\nu L1 norm is: %.16f", l1_u); // reporting L1
	printf("\nv L1 norm is: %.16f\n", l1_v); // reporting L1
	
	printf("\np L2 norm is: %.16f", l2_p); // reporting L2
	printf("\nu L2 norm is: %.16f", l2_u); // reporting L2
	printf("\nv L2 norm is: %.16f\n", l2_v); // reporting L2
	
	printf("\np L_inf norm is: %.16f", linf_p); // reporting L inf
	printf("\nu L_inf norm is: %.16f", linf_u); // reporting L inf
	printf("\nv L_inf norm is: %.16f\n", linf_v); // reporting L inf
	
	printf("\n Main vortex location is: %f, %f, %d, %d\n", (loc[0]-0.5)*dx,-h+(loc[1]-0.5)*dy,loc[0], loc[1]);
	printf("\n Gamma: %f", gamma(vrt, 0.12, 15, 24, dx, dy));
	
	fclose(f1); // closing the data storage file
	fclose(f2); // closing the data storage file

}
