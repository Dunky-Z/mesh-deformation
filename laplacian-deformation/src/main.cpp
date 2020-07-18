#include "Laplacian.h"

int main(int argc, char **argv)
{
	double start = GetTickCount();

	LaplaceDeformation Deform;
	Deform.AllpyLaplaceDeformation(argv);


	
	double finish = GetTickCount();
	cout << finish - start << endl; 
	return 0;
}