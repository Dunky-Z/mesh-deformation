#include "stdafx.h"
#include "Laplacian1.h"

int main(int argc, char **argv)
{
	double start = GetTickCount();  //开始时间
	LaplaceDeformation Deform;
	Deform.AllpyLaplaceDeformation(argv);
	double finish = GetTickCount();   //结束时间
	double t = finish - start;
	cout << t << endl; //输出时间
	return 0;
}