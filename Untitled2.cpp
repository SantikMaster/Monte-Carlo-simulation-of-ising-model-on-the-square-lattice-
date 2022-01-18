#include <iostream>
#include <math.h>
#include <stdio.h>

#define M_PI       3.14159265358979323846
#include <cmath>

#define SIZE 30
#define a_lattice 10.0

#define LEAPS 5000
#define CYCLES 100
 
// this program provides calculation for 2D ising system using Monte Carlo alorhythm
// You will get the file which is M(h) Magnetization as function of exterlan fiels
// Or you can get a spin configuration for any field or Temperature that you want.
// Then you get file in the format XYAM - XY - are the coordinates of the spin and A- angle M - magnitude


// The framework is Ising model on square lattice. The random interaction (RandCof) can be included and long range interaction (J_d) ~1/r^2



using namespace std;

	double h = 0.3;
	double T=0.001;
	double J_d = 9; //500
	double RandCof = 1.0;
	
	long double Hamiltonian[8];
	long double fn();
	
	struct Spin
	{
		double x;
		double y;
		int Sz;
	};
	
	Spin S[SIZE][SIZE];
	double J[SIZE][SIZE][2];
	
void randHamiltonianFill();
void TextOut(string Text, int Number);
double distance (Spin, Spin);
double EnergyCalc(int, int);
double Magnetization();
void MonteCarlo();
void PrintSpins(bool ToFile);
void PrintSpins_XYAM( bool ToFile);

int main(int argc, char** argv) 
{
	
	int i, j;
	  double M;
	randHamiltonianFill();

//	while (h<2)
	{
	randHamiltonianFill();	

	for (i = 0; i<CYCLES; i++)
	{	
		for (j = 0; j<LEAPS; j++)
	{	
	     MonteCarlo();
    }
	
    }cout << "!!!!!!!!!!\n\n\n";

	char buf[256];

		  M = Magnetization();
			sprintf(buf, "%4.4f    %4.4f  \n", h, M);

			
			string ST(buf);
			TextOut(ST, 0);
			cout << ST;

     h+=0.05;
   }
	PrintSpins_XYAM(1);
//	PrintSpins(1);


//    getchar();	
	return 0;
}
void PrintSpins( bool ToFile)
{

	char buf[256];


	int i, j;
	double Val;
	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{ 
			Val = S[i][j].Sz;
			sprintf(buf, "%10d    %10d  %1.1f  \n  ", i, j, Val);
			string ST(buf);
			if(ToFile == 0)
			{
				cout << ST;
			}
			else TextOut(ST, 0);
		}

		string ST(buf);
		if(ToFile == 0)
		{
				cout << ST;
		}
		else TextOut(ST, 0);
	}
}
void PrintSpins_XYAM( bool ToFile)
{
// angle magnitude
	char buf[256];


	int i, j;
	double Val;
	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{ 
			Val = M_PI/2*S[i][j].Sz;
			sprintf(buf, "%10d    %10d  %1.1f  1  \n  ", i, j, Val);
			string ST(buf);
			if(ToFile == 0)
			{
				cout << ST;
			}
			else TextOut(ST, 0);
		}
	//	sprintf(buf, "\n");
		string ST(buf);
		if(ToFile == 0)
		{
				cout << ST;
		}
		else TextOut(ST, 0);
	}
}
void TextOut(string Text, int Number)
{

	// #include <stdio.h>

	FILE	 	*stream;
	char str[10];


	char Stri[25];


	stream =
		fopen("DataEnh.txt",  "a");
fopen(Stri,  "a");


	//  "Data.txt"
	char * p;
	p = (char *)Text.c_str();
	if (stream != NULL)
	{
		fprintf(stream, p, "DataEnh.txt");

		fclose(stream);
	}
	return;
}

double Magnetization()
{ 
	int i, j;
	double M = 0;
   for (i = 0; i<SIZE; i++)
	{
		for (j = 0; j<SIZE; j++)
		{
			M+= S[i][j].Sz;
 
		}

	}
	M/=SIZE*SIZE;
return M;
}
void randHamiltonianFill()
{
	int i, j;
	for (i = 0; i<SIZE; i++)
	{
		for (j = 0; j<SIZE; j++)
		{
		
			S[i][j].x = a_lattice*i; 
			S[i][j].y = a_lattice*j; 	
			S[i][j].Sz = 1-2*(rand()%2);
			
			J[i][j][0] =-1 +RandCof*(1-0.01*(rand()%201)); // -1 + рандом от  ot 1 RandCof до 1 RandCof
			J[i][j][1] =-1 +RandCof*(1-0.01*(rand()%201)); // 	

		}

	}
	
}
double distance (Spin S1, Spin S2)
{
	double x1, y1, R;
	
	x1 = S1.x- S2.x;
	y1 = S1.y- S2.y;
	
	R = x1*x1+y1*y1;
	if (R!=0) R= sqrt(R);
	return R;
	
}
double EnergyCalc(int i, int j)
{
double Energy = 0, dist;
	
	if (i!= SIZE)
    {
		Energy+= J[i][j][0]*S[i][j].Sz*S[i+1][j].Sz;
	}
	if (j!= SIZE)
    {
		Energy+= J[i][j][1]*S[i][j].Sz*S[i][j+1].Sz;
	}
	if (i!= 0)
    {
		Energy+= J[i-1][j][0]*S[i][j].Sz*S[i-1][j].Sz;
	}
	if (j!= 0)
    {
		Energy+= J[i][j-1][1]*S[i][j].Sz*S[i][j-1].Sz;
	}
	
		int ii, jj;
	for (ii = 0; ii<SIZE; ii++)
	{
		for (jj = 0; jj<SIZE; jj++)
		{
        	dist = distance (S[i][j], S[ii][jj]);
			if (dist!=0) Energy += J_d/dist/dist*S[i][j].Sz*S[ii][jj].Sz;
		}
   }
Energy -= h*S[i][j].Sz;		
			
return Energy;	
		
}
void MonteCarlo()
{
	int i, j;
	double OldEnergy, NewEnergy;
	long double W;
	i = rand()%SIZE;
	j = rand()%SIZE;
	OldEnergy =  EnergyCalc(i, j);
	S[i][j].Sz = -1*S[i][j].Sz;
	NewEnergy = EnergyCalc(i, j);
	if (OldEnergy < NewEnergy)
	{
		S[i][j].Sz = -1*S[i][j].Sz;	
		
		W = exp (-(NewEnergy-OldEnergy )/T);

		if ((rand()%1001)/1000 < W)	
		{
				S[i][j].Sz = -1*S[i][j].Sz;		
		}
	}
	
	
	
}

