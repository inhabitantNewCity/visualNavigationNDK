#ifndef data_HH
#define data_HH
#pragma once
#include <stdio.h>

//extern  double *Xe,*Ye,*Se;
//extern double StepEtalX,StepEtalY,StepObjX,StepObjY,StepRestX,StepRestY;
//extern int TypeScene,Ny,ThrEtal,ThrObj;
//extern int NumEtal,KetImage;
//extern double X[1024],Y[1024],S[1024];
//extern long LenFile;
//extern long XYEtal[16384],XYObj[16384];

// объединение для точки
 union uXY {
	 unsigned int XY;
	 unsigned short Ptr[2];
};
union ul { 
	long ll; 
	float fl; 
	unsigned short Sh[2];
} ;
// заголовок эталона
struct hEtalon {
	unsigned short N;           // количество точек
	unsigned short Class;       // класс
};

// точка
struct sPoint {
	unsigned short X;
	unsigned short Y;
};

struct CapEtalons {
	unsigned short N;
        unsigned short Class;
};
struct ParmEtal {
	long NbP;
	unsigned short n;
	double Xn;
	double Yn;
        double Xk;
	double Yk;
	double DE;
	double Sq;
	double Qm;
        double EM;
};
struct DataRecObj {
	unsigned short Xn:15, SgnDx:1;
	unsigned short Yn:15, SgnDy:1;
	unsigned short NE:15, Inet:1;
	unsigned short dx:8, dy:8;
};
#endif
