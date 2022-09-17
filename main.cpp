#pragma warning (disable : 4996)

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>      


#define M_PI       3.14159265358979323846
#include <cmath>

#define a_lattice 10.0

#define LEAPS 5000
#define CYCLES 10
 
// this program provides calculation for 2D ising system using Monte Carlo alorhythm
// You will get the file which is M(h) Magnetization as function of exterlan fiels
// Or you can get a spin configuration for any field or Temperature that you want.
// Then you get file in the format XYAM - XY - are the coordinates of the spin and A- angle M - magnitude
	
struct Spin
{
	double x;
	double y;
	int Sz;
};

void TextOut(const std::string &Text, int Number);
double Distance (Spin, Spin);

template <const int x, const int y>
	void PrintSpins(bool ToFile, Spin (&S)[x][y]);
template <const int x, const int y>
	void PrintSpins_XYAM(bool ToFile, Spin (&S)[x][y]);
template <const int x, const int y>
	double Magnetization(Spin (&S)[x][y]);
template <const int x, const int y>
	void RandSpinFill(Spin (&S)[x][y], double (&J)[x][y][2], double RandCof);
template <const int x, const int y>
	double EnergyCalc(int i, int j, Spin (&S)[x][y], double (&J)[x][y][2], double J_d, double h);
template <const int x, const int y>
	void MonteCarlo(Spin (&S)[x][y], double (&J)[x][y][2], double T,  double J_d, double h);

int main(int argc, char** argv) 
{
	double h = 0.3;
	double T=0.001;
	double J_d = 9; //500
	double RandCof = 1.0;
	
	const int x = 30, y = 20;
	
	Spin S[x][y];
	double J[x][y][2];
	
	int i, j;
	double M;
	RandSpinFill(S, J, RandCof);

	while (h<2)
	{
		RandSpinFill(S, J, RandCof);

		for (i = 0; i<CYCLES; i++)
		{	
			for (j = 0; j<LEAPS; j++)
			{	
	  	  	 	MonteCarlo(S, J, T, J_d, h);
    		}
   		}
		std::cout << "\n";

		M = Magnetization(S);

		std::stringstream ss;
		ss << h << '\t' << M;
	
			
		std::string ST = ss.str();
		TextOut(ST, 0);
		std::cout << ST;

     	h+=0.05;
    }
	PrintSpins_XYAM(true, S);
	
	return 0;
}

template <const int x, const int y>
void PrintSpins(bool ToFile, Spin (&S)[x][y])
{
	char buf[256];

	int i, j;
	double Val;
	for (i = 0; i < x; i++)
	{
		for (j = 0; j < y; j++)
		{ 
			Val = S[i][j].Sz;
			sprintf(buf, "%10d    %10d  %1.1f  \n  ", i, j, Val);
			std::string ST(buf);
			if(ToFile == 0)
			{
				std::cout << ST;
			}
			else TextOut(ST, 0);
		}

		std::string ST(buf);
		if(ToFile == 0)
		{
				std::cout << ST;
		}
		else TextOut(ST, 0);
	}
}

template <const int x, const int y>
void PrintSpins_XYAM(bool ToFile, Spin (&S)[x][y])
{
	int i, j;
	double Val;
	for (i = 0; i < x; i++)
	{
		for (j = 0; j < y; j++)
		{ 
			Val = M_PI/2*S[i][j].Sz;

			std::stringstream ss;
			ss << i << '\t' << j << '\t' << Val << '\n';
			
			std::string ST = ss.str();
			if(ToFile == 0)
			{
				std::cout << ST;
			}
			else TextOut(ST, 0);
		}
	}
}
void TextOut(const std::string &Text, int Number)
{
	// #include <stdio.h>
	FILE	 	*stream;

	stream =
		fopen("DataEnh.txt",  "a");

	char * p;
	p = (char *)Text.c_str();
	if (stream != NULL)
	{
		fprintf(stream, p, "DataEnh.txt");

		fclose(stream);
	}
	return;
}

template <const int x, const int y>
double Magnetization(Spin (&S)[x][y])
{ 
	int i, j;
	double M = 0;
   for (i = 0; i<x; i++)
	{
		for (j = 0; j<y; j++)
		{
			M+= S[i][j].Sz;
		}
	}
	M/=x*y;
return M;
}

template <const int x, const int y>
void RandSpinFill(Spin (&S)[x][y], double (&J)[x][y][2], double RandCof)
{
	int i, j;
	for (i = 0; i<x; i++)
	{
		for (j = 0; j<y; j++)
		{
			S[i][j].x = a_lattice*i; 
			S[i][j].y = a_lattice*j; 	
			S[i][j].Sz = 1-2*(rand()%2);
			
			J[i][j][0] =-1 + RandCof*(1-0.01*(rand()%201)); // -1 + рандом от  ot 1 RandCof до 1 RandCof
			J[i][j][1] =-1 + RandCof*(1-0.01*(rand()%201)); // 	
		}
	}
}

double Distance (Spin S1, Spin S2)
{
	double x1, y1, R;
	
	x1 = S1.x- S2.x;
	y1 = S1.y- S2.y;
	
	R = x1*x1+y1*y1;
	if (R!=0) R = sqrt(R);
	return R;
	
}

template <const int x, const int y>
double EnergyCalc(int i, int j, Spin (&S)[x][y], double (&J)[x][y][2], double J_d, double h)
{
double Energy = 0, dist;
	
	if (i!= x)
    {
		Energy += J[i][j][0]*S[i][j].Sz*S[i+1][j].Sz;
	}
	if (j!= y)
    {
		Energy += J[i][j][1]*S[i][j].Sz*S[i][j+1].Sz;
	}
	if (i!= 0)
    {
		Energy += J[i-1][j][0]*S[i][j].Sz*S[i-1][j].Sz;
	}
	if (j!= 0)
    {
		Energy += J[i][j-1][1]*S[i][j].Sz*S[i][j-1].Sz;
	}
	
	int ii, jj;
	for (ii = 0; ii<x; ii++)
	{
		for (jj = 0; jj<y; jj++)
		{
        	dist = Distance (S[i][j], S[ii][jj]);
			if (dist!=0) Energy += J_d/dist/dist*S[i][j].Sz*S[ii][jj].Sz;
		}
    }
	Energy -= h*S[i][j].Sz;		
			
	return Energy;	
}

template <const int x, const int y>
void MonteCarlo(Spin (&S)[x][y], double (&J)[x][y][2], double T, double J_d, double h)
{
	int i, j;
	double OldEnergy, NewEnergy;
	long double W;
	i = rand()%x;
	j = rand()%y;
	OldEnergy =  EnergyCalc(i, j, S, J, J_d, h);
	S[i][j].Sz = -1*S[i][j].Sz;
	NewEnergy = EnergyCalc(i, j, S, J, J_d, h);
	if (OldEnergy < NewEnergy)
	{
		S[i][j].Sz = -1*S[i][j].Sz;	
		
		W = exp (-(NewEnergy-OldEnergy )/T);

		if (rand() / (RAND_MAX + 1.0) < W)
		{
				S[i][j].Sz = -1*S[i][j].Sz;		
		}
	}
}

