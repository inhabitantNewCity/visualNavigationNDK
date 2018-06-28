#include "AutoRecognition.h"
#include "TPoint.h"
#include "TContour.h"
#include "data.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>




FILE* HederFileKmp,HederFileEtalons,HederFileCompressKmp,HederFileStandard;
char NameFileKmp[256],NameFileEtalons[256],NameFileCompressKmp[256],*Ptr,
     Path[256],CurrentDirectory[256],NameFileStandard[256];
int Pmax,Pmin,Ny=20,Nb,NbPk,NumParm,OutMin,InMin,ThrDO,ThrBO;
struct CapEtalons CEC,*CEM;
bool SelectData;
double Xb[10],Yb[10],Sb[10],DWb,SE,SimPor;
unsigned int *BigObjects,*ZME,*NumSmall5;
char text[256];
long ANE,APE,XDmax,XDmin,YDmax,YDmin,N_Dmax,AdrBig;
long SizeFileKmp,NumWordKmp,*AllKmp,*CapKmp,*PaspKmp,*MetrKmp,*XYObj;
int KSmObj,SMSmObj,KBObj,SMBObj,KDObj,SMDObj,KPObj,SMPObj,SMEtal;
int KSmObj3,KSmObj4,KSmObj5,NObj;
long SizeBigObjects,SizeDiscreteObjects,SizeEtalons,SizeSmallObjects,
     MoveCompressKmp;
int NumStop;
//AnsiString AnsiStr;
int ErrorCode;
int i, j;
struct DataRecObj *DiscreteObjects;

  union {long ll; float fl;} ul;

// сжатие целочисленного описания объекта
// xr[] - координаты по x
// yr[] - координаты по y
// n - количество точек
// d - порог сжатия кусочно-линейного
// x[] - выходные координаты
// y[] -
// возвращает количество точек
int CompressInteger(long xr[],long yr[],int n,int d,long x[],long y[])
  { int xn,yn,xk,yk;
  int m11,j11,k11,i11,jt,jj11;
  long xt,yt,d11,e11,s11,st11,p11,rn11,rk11;
  d11=(long)d*d;   xn=xr[0]; yn=yr[0];   x[0]=xn; y[0]=yn;
  m11=1;j11=1; k11=2;
  if(n==1) return(n);
  while((i11=j11+k11-1)<n)
    { xk=xr[i11]; yk=yr[i11]; xt=(long)xk-xn; yt=(long)yk-yn;
    s11=(xt*xt+yt*yt); p11=d11*s11;
    for(jt=1;jt<k11;jt++)
      { jj11=j11+jt-1;
      e11=((long)xr[jj11]-xn)*yt-((long)yr[jj11]-yn)*xt;
      if(e11*e11>p11) break;
      else
        { st11=((long)xr[jj11]-xn)*xt+((long)yr[jj11]-yn)*yt;
        if(st11<=s11 && st11>0) continue;
        else
          { if(st11<=0)
            { rn11=((long)xr[jj11]-xn)*((long)xr[jj11]-xn)+
                   ((long)yr[jj11]-yn)*((long)yr[jj11]-yn);
            if(rn11<=d11) continue;  else break;
            }
          else
            { rk11=((long)xr[jj11]-xk)*((long)xr[jj11]-xk)+
                   ((long)yr[jj11]-yk)*((long)yr[jj11]-yk);
            if(rk11<=d11) continue;  else break;
            }
          }
        }
      }
    if(jt>=k11) { k11++; continue;}
    else
      { j11+=k11-1; xn=xr[j11-1]; yn=yr[j11-1];
      x[m11]=xn; y[m11]=yn; k11=2; m11++;
      }
    }
  x[m11]=xr[n-1]; y[m11]=yr[n-1];  m11++;
  return (m11);
  }

//---------------------------------------------------------------------
// сжатие точечного описания объекта
// Zi[] - исходные координаты точек
//int N - количество точек
//int d - порог
// long Ze - выходные координаты
// возвращает количетсво точек
int CompressPoints(long Zi[],int N,int d,long Ze[])
  { long xn,yn,xk,yk,xz,yz,xp,yp;
  int m,j,k,i,jt,jp,j2,jp2,mp;
  long xt,yt,d2,e,s,st,rn,rk;
  __int64 e2,p;
  bool Ibreak;
  xn=Zi[0]; yn=Zi[1];
  Ze[0]=xn;   Ze[1]=yn;
  if(N==1) return(N);

  m=1; i=1; k=0; j=1;
  d2=(long)d*d;
  for(;j<N;)
    {
    j2=j*2;
    xk=Zi[j2]; yk=Zi[j2+1];  Ibreak=false;
    if(k>0)
      {
      xp=xk-xn; yp=yk-yn;
      s=(xp*xp+yp*yp);
      p=(__int64)s; p*=d2;
      for(jt=0; jt<k; jt++)
        { jp=i+jt; jp2=jp*2;
        xt=Zi[jp2]-xn; yt=Zi[jp2+1]-yn;
        e=xt*yp-yt*xp;
        e2=(__int64)e*(__int64)e;
        if(e2>p) Ibreak=true;
        else
          { st=xp*xt+yp*yt;
          if(st<=s && st>0) continue;
          else
            {
            if(st<=0)
              { rn=xt*xt+yt*yt;
              if(rn<=d2) continue;
              else Ibreak=true;
              }
            else
              { xt-=xp; yt-=yp;
              rk=xt*xt+yt*yt;
              if(rk<=d2) continue;
              else Ibreak=true;
              }
            }
          }
        if(Ibreak==true) break;
        }
      }
    if(Ibreak==false) { k++; j++; xz=xk; yz=yk; continue;}
    else
      { i+=k; xn=xz; yn=yz;    mp=m*2;
      Ze[mp]=xz; Ze[mp+1]=yz; k=0; m++; j=i;
      }
    }
  mp=m*2; Ze[mp]=xz; Ze[mp+1]=yz;  m++;

  if(m==2)
    {
    if(Ze[0]==Ze[2]&&Ze[1]==Ze[3]) m=1;
    }

  return (m);
  }

//---------------------------------------------------------------------
// X, - исходные координаты
// Y, 
// S, - расстояние между точками в описании конутров
// Nt - количество точек в описании контура
// Np - номер точки, определяющей начальную точку описания контура
// Sp - смещение начальной точки нового описания от точки с номером Np
// Se - длина базового контура
// выходные параметры
// Xr - текущее описание по x
// Yr,
// Sr - расстояние между точками (новое)
// возвращает количетсво точек
int CurrentDecribeObject( double *X, double *Y, double *S,
                   int Nt, int Np, double Sp, double Se,
                   double *Xr, double *Yr, double *Sr)
  {
  int i,n,i1,i2;
  double sg,Sfr,st;

  if(S[Np-1]<Sp) return(-1);
  i=Np;
  i2=i%Nt; i1=(i-1)%Nt;
  sg=Sp/S[i1];
  Xr[0]=X[i1]*(1.-sg)+X[i2]*sg;
  Yr[0]=Y[i1]*(1.-sg)+Y[i2]*sg;
  Sfr=S[i1]-Sp;
  st=Sfr;
  n=1;
  for(;;)
    { i2=i%Nt; i1=(i-1)%Nt;
    if(Sfr>=Se)
      {
      sg=Sfr-Se;
      Sr[n-1]=st-sg;
      sg/=st;
      Xr[n]=X[i1]*sg+X[i2]*(1.-sg);
      Yr[n]=Y[i1]*sg+Y[i2]*(1.-sg);
      n++;
      break;
      }

    Xr[n]=X[i2]; Yr[n]=Y[i2]; Sr[n-1]=st;
    st=S[i2]; Sfr+=st;
    i++; n++;
    }

  return(n);
  }

//---------------------------------------------------------------------
#define EPS 1.e-5
// Ne
// Xo - координаты исходного описания эталона
// Yo
// So - расстояние между точками эталона
// X - координаты исходного описания объекта
// Y,
// S - расстояние между точками объекта
// выходные параметры
// XXo - вспомогательное описание эталона по x
// YYo
// XX - - вспомогательное описание объекта по x
// YY
// C[] - расстояния между точками эталона и объекта
int AuxiliaryDescriptionShapes(int Ne,
        double *Xo, double *Yo, double *So,
        double *X, double *Y, double *S,
        double *XXo, double *YYo, double *XX, double *YY, double C[])
  {
  int in,jn,NM,nt,jt,N,ii;
  double sp,spo,sg;
  if(Xo[0]!=Xo[Ne-1]||Yo[0]!=Yo[Ne-1]) N=Ne+1;
  else N=Ne;
  in=1; jn=1; NM=1;
  sp=S[0];
  XXo[0]=Xo[0]; YYo[0]=Yo[0];
  XX[0]=X[0];   YY[0]=Y[0];
  spo=So[0];
  while (in<N)
    { ii=in%Ne;
    if(sp<EPS) { sp+=S[jn]; jn++; }
    if(spo>=EPS)
      {
      if(spo<=sp)
        { nt=NM-1;
        XXo[NM]=Xo[ii]; YYo[NM]=Yo[ii];
        jt=jn; sg=spo/sp;
        XX[NM]=XX[nt]*(1.-sg)+X[jt]*sg;
        YY[NM]=YY[nt]*(1.-sg)+Y[jt]*sg;
        C[nt]=spo; sp-=spo; spo=So[ii];
        NM++; in++;
        }
      else
        { nt=NM-1; jt=jn;
        XX[NM]=X[jt]; YY[NM]=Y[jt];
        sg=sp/spo;
        XXo[NM]=XXo[nt]*(1.-sg)+Xo[ii]*sg;
        YYo[NM]=YYo[nt]*(1.-sg)+Yo[ii]*sg;
        C[nt]=sp; spo-=sp; sp=S[jt];
        NM++; jn++;
          }
        }
      else { spo+=So[ii]; in++; }
      }

  return(NM);
  }

//---------------------------------------------------------------------
// Xo - координаты вспомогательного опиcания эталона по x 
// Yo
// X - координаты вспомогательного опиcания эталона по x
// Y
// C - расстояние между точками
// Se - длина базового контура (20 точек)
// N - количество точек
// Sn - ковареация
// Cs - ковареация
// Re - сумма квадратов (Sn^2 + Re^2)
void  ErrorMatching( double *Xo, double *Yo, double *X, double *Y,
                double *C, double Se, int N,
                double *Sn, double *Cs, double *Re)
  {
  double RXYo,RYYo,RXXo,RYXo,px,ppx,py,ppy,xp,yp,zx,zy,zzx,zzy,
         xz,yz,sd,sn,cs,re,Se6,st;
  int j;
  Se6=Se*6.;
  px=Xo[0]; py=Yo[0];
//  zx=X[0]; zy=Y[0];
  zx=0.; zy=0.;
  RXYo=0.; RYYo=0.; RXXo=0.; RYXo=0.;

  for(j=1; j<N; j++)
    { ppx=Xo[j]; ppy=Yo[j]; st=C[j-1];
    xp=px+ppx; yp=py+ppy;
    zzx=X[j]-X[0]; zzy=Y[j]-Y[0];
    xz=zx+zzx; yz=zy+zzy;
    sd=xp*xz; sd=sd+sd;
    RXXo+=(sd-px*zzx-ppx*zx)*st;
    sd=xp*yz; sd=sd+sd;
    RYXo+=(sd-px*zzy-ppx*zy)*st;
    sd=yp*xz; sd=sd+sd;
    RXYo+=(sd-py*zzx-ppy*zx)*st;
    sd=yp*yz; sd=sd+sd;
    RYYo+=(sd-py*zzy-ppy*zy)*st;
    px=ppx; py=ppy;
    zx=zzx; zy=zzy;
    }
  sn=(RYXo-RXYo)/Se6; *Sn=sn;
  cs=(RXXo+RYYo)/Se6; *Cs=cs;

  re= sn*sn+cs*cs;
  *Re=re;

  return;
  }

//---------------------------------------------------------------------
#include <math.h>
#define EPS 1.e-5
// NumEtal - количество эталонов
// ParmEtals - структура, описывающая эталон
// StepX - dpi
// StepY - dpi
// SE - длина базового контура
// выходные параметры:
// Xe - описание эталона
// Ye
// Se - длина контура
void ReadEtalons(int NumEtal, struct ParmEtal *ParmEtals,
                 double StepX, double StepY, double SE,
                 double *Xe, double *Ye, double *Se)

				 // StepX, double StepY - размеры сетки (зависит от dpi) (посмотреть dpi!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! перевод)
				 // длина эталона приводится к размеру = 20
  {
  int n;
  long NbP,NbPi,NbPj,                         // размер смещения при чтении данных из файла 
	 Xmax,Xmin, Ymax,Ymin,ix,iy;
  double Xt,Yt,Xp,Yp,St,So,Sq,MX,MY,MX2,MY2,DX,DY,DW,DWo,Xa,Ya;
  union {unsigned int XY; unsigned short Ptr[2];} CR;
  double Xr[1025],Yr[1025],
	  Sr[1025],                                  // расстояния между точками
	  Xbv[2050],Ybv[2050],                       // для записи вспомогательного описания базового эталона (треугольник эталона масштабированный на 20)
         Xv[2050],Yv[2050],Sv[2050];             // для записи вспомогательного описания эталона
                                                 // !!!!!!!!!!!!!!метод Кифера для нахождения min (глава 10)
  double Em,Ds,Sp,Re,Cs,Sn,Emin[21],Qmin[21],Qmax[21];
  int Nt,Np,Nmin,Nmax,Class;

  NbPk=0;  NbP=0;  NbPi=0;
  CEM=(struct CapEtalons *)ZME;
  for(i=0; i<NumEtal; i++)
    {
    n=CEM[NbPi].N;
    Class=CEM[NbPi].Class;
    NbPi++;

    So=0.;
    for(j=0; j<n; j++)
      { NbPj=NbP+j;
      CR.XY=ZME[NbPi+j];
      ix=CR.Ptr[0]; iy=CR.Ptr[1];
	  // получение эталонного значения в миллиметрах
      Xe[NbPj]=Xt=ix*StepX;
      Ye[NbPj]=Yt=iy*StepY;
      if(j==0)
        { Xmax=ix; Xmin=ix;
        Ymax=iy; Ymin=iy;
        }
      else
        {
        if(Xmax<ix) Xmax=ix;
        if(Xmin>ix) Xmin=ix;
        if(Ymax<iy) Ymax=iy;
        if(Ymin>iy) Ymin=iy;
		// вычисление расстояния между точками
        Se[NbPj-1]=sqrt((Xt-Xp)*(Xt-Xp)+(Yt-Yp)*(Yt-Yp));
		// длина эт в мил
        So+=Se[NbPj-1];
        }
      Xp=Xt; Yp=Yt;
      }
    Se[NbPj]=sqrt((Xt-Xe[NbP])*(Xt-Xe[NbP])+(Yt-Ye[NbP])*(Yt-Ye[NbP]));
    So+=Se[NbPj];

//Вычисление статистических характеристик эталона
    MX=0.;  MY=0.;
    MX2=0.; MY2=0.;
    Xp=Xe[NbP]; Yp=Ye[NbP];
	// вычисление корреляционных моментов
    for(j=1; j<n; j++)
      { NbPj=NbP+j;
      Xt=Xe[NbPj]; Yt=Ye[NbPj]; St=Se[NbPj-1];
      Xa=Xt+Xp; Ya=Yt+Yp;
      MX+=Xa*St; MY+=Ya*St;
      MX2+=(Xa*Xa-Xt*Xp)*St; MY2+=(Ya*Ya-Yt*Yp)*St; // начальный момент
      Xp=Xt; Yp=Yt;
      }
    Xt=Xe[NbP]; Yt=Ye[NbP];  St=Se[NbPj];
    Xa=Xt+Xp; Ya=Yt+Yp;
    MX+=Xa*St; MY+=Ya*St;
    MX/=(2.*So); MY/=(2.*So);
    MX2+=(Xa*Xa-Xt*Xp)*St; MY2+=(Ya*Ya-Yt*Yp)*St;
    MX2/=(3.*So); MY2/=(3.*So);
    // дисперчия по x
	DX=MX2-MX*MX;
	DY=MY2-MY*MY;
    DW=DX+DY;
    Sq=SE/So;                        // коэффициент масштабирования эталона относительно базового
    ParmEtals[i].Sq=Sq;
    DWo=DW*Sq*Sq;                    // величина дисперсии привиденного эталона
    ParmEtals[i].DE=DWo;
    ParmEtals[i].Xn=-MX;
    ParmEtals[i].Yn=-MY;
    ParmEtals[i].Xk=(Xmax-Xmin)*StepX-MX;       // конечное значение
    ParmEtals[i].Yk=(Ymax-Ymin)*StepY-MY;
//центрирование и масштабирование (на SE) описаний эталонов
    for(j=0; j<n; j++)
      {
	  // масштабирование
	  NbPj=NbP+j;
      Xe[NbPj]-=MX; Ye[NbPj]-=MY;
      Xe[NbPj]*=Sq; Ye[NbPj]*=Sq; Se[NbPj]*=Sq;
      }
    ParmEtals[i].n=(unsigned short)n;
    ParmEtals[i].NbP=NbP;
//======== Нахождение экстремумов =========
    Em=DWo; 
	Ds=SE/Ny;                        // шаг по контуру
	Sp=-Ds;
	St=Se[NbP];
	Np=1;
    // расчет ошибок от начальной до конечной с шагом
	for(j=0; j<Ny; j++)
      { 
		  Sp+=Ds;
      for(;St<Sp;) {Sp-=St; St=Se[Np+NbP]; Np++; }
	  // текущее описание объекта
      CurrentDecribeObject(&Xe[NbP],&Ye[NbP],&Se[NbP],n,Np,Sp,SE,Xr,Yr,Sr);
	  //Np,Sp -номер точки, смещение

	  // вспомогательное описание форм
      Nt=AuxiliaryDescriptionShapes(Nb,Xb,Yb,Sb,Xr,Yr,Sr,Xbv,Ybv,
                                   Xv,Yv,Sv);
	  // ошибка сходства
      ErrorMatching(Xbv,Ybv,Xv,Yv,Sv,SE,Nt,&Sn,&Cs,&Re);
	  // Xbv,Ybv,  вспомогательное описание объекта
 	  // Xv,Yv
      Emin[j]=DWb-Re/DWo;
      }
    Nmin=0;  Nmax=0;
    if(Emin[0]<Emin[Ny-1]&&Emin[0]<=Emin[1])
      { Qmin[Nmin]=0.; Nmin++;
      }
    if(Emin[0]>Emin[Ny-1]&&Emin[0]>=Emin[1])
      { Qmax[Nmax]=0.; Nmax++;
      }
    for(j=1; j<Ny-1; j++)
      {
      if(Emin[j]<Emin[j-1]&&Emin[j]<=Emin[j+1])
        { Qmin[Nmin]=j*Ds; Nmin++;
        }
      if(Emin[j]>Emin[j-1]&&Emin[j]>=Emin[j+1])
        { Qmax[Nmax]=j*Ds; Nmax++;
        }
      }
    if(Emin[Ny-1]<Emin[Ny-2]&&Emin[Ny-1]<=Emin[0])
      { Qmin[Nmin]=(Ny-1)*Ds; Nmin++;
      }
    if(Emin[Ny-1]>Emin[Ny-2]&&Emin[Ny-1]>=Emin[0])
      { Qmax[Nmax]=(Ny-1)*Ds; Nmax++;
      }
//========= поиск наилучшего минимума ===========
    Nt=Qmin[0]/Ds+0.5;
    Em=Emin[Nt];
    ParmEtals[i].Qm=Qmin[0];
    for(j=1; j<Nmin; j++)
      {
      Nt=Qmin[j]/Ds+0.5;
      if(Emin[Nt]<Em)
        { Em=Emin[Nt];
        ParmEtals[i].Qm=Qmin[j];
        }
      }
//======= Определение порога для поиска согласованных описаний =======
    Np=0;  Sp=Ds;
    for(j=0; j<n; j++)
      { Np=j+1;
      if(Se[NbP+j]>Sp) break;
      Sp-=Se[NbP+j];
      }
    CurrentDecribeObject(&Xe[NbP],&Ye[NbP],&Se[NbP],n,Np,Sp,SE,
                          Xr,Yr,Sr);
    Nt=AuxiliaryDescriptionShapes(n,&Xe[NbP],&Ye[NbP],&Se[NbP],Xr,Yr,Sr,
                                  Xbv,Ybv,Xv,Yv,Sv);
    ErrorMatching(Xbv,Ybv,Xv,Yv,Sv,SE,Nt,&Sn,&Cs,&Re);
    Em=DWo-Re/DWo;
    ParmEtals[i].EM=Em;
//====================================================================
    NbP+=n;
    NbPi+=n;
    }
  NbPk=NbP;

  return;
  }

//---------------------------------------------------------------------
void RecognitionObject(TScene *SceneIn, long KObj, int *NumEtal,  struct ParmEtal *ParmEtals,
                 double StepX, double StepY, double SE,
                 double *Xe, double *Ye, double *Se)
  { // int FileEtalons, int FileKMP,
  int n,k,jt,N,Nt,Nsm,Nln,Np,Ne,Thr,LP=16;
  long NbP,Xmax,Xmin, Ymax,Ymin,ZMo[2048];
  long Xi[1025],Yi[1025];
  double Xd[1025],Yd[1025],Sd[1025];
  double Xt,Yt,Xp,Yp,St,So,Sq,MX,MY,MX2,MY2,DX,DY,DW,DWo,DE,Xa,Ya;
  long ii,NumObj,Lobj;
  double Em,Rm,Eq,Rq,Qmin[21],Emin[21],Qmax[21],Ds,Sp,Sn,Cs,Re,Qt;
  double Snq,Csq,Snm,Csm,Xc,Yc,Hg,Wd;
  int Nmin,Nmax,NEtal;
  double Xr[1025],Yr[1025],Sr[1025], Xbv[2050],Ybv[2050],
         Xv[2050],Yv[2050],Sv[2050];
//======= переменные уточнения минимума невязки ==========
  double s0,s1,s2,sc,e0,e1,e2,Er,emax=0.5,Emax,step=0.5,w1,w2;
  int ol;

  TContour *wContour;
  T_Point maxXY;
  T_Point minXY;
  T_Point wPoint;

  KPObj=0; SMPObj=0;  KBObj=0; SMBObj=0;
  KDObj=0; SMDObj=0;  KSmObj=0; SMSmObj=0;
  KSmObj3=0; KSmObj4=0; KSmObj5=0;

  for(ii=0; ii<KObj; ii++)
    { Nsm=0; Nln=0;

    itoa(ii+1, text, 10);
    
	//AP=ii*LP;

    wContour = SceneIn->GetContour(ii);

    //NumObj=PaspKmp[AP+4l];
	// номер объекта
	NumObj = ii;
	maxXY = wContour->FindMaxXY();
	minXY = wContour->FindMinXY();
    //Xmin=PaspKmp[AP+5l];
	Xmin = minXY.GetX();
    //Ymin=PaspKmp[AP+7l];
	Ymin = minXY.GetY();
    //Xmax=PaspKmp[AP+6l];
	Xmax = maxXY.GetX();
    //Ymax=PaspKmp[AP+8l];
	Ymax = maxXY.GetY();
    Lobj=wContour->GetNumPoints()*2;        // количество признаков
	
    n=Lobj/2;                   // количество точек

//===== отсечение по не попаданию в прямоугольник =====
    if(Xmin<XDmin) continue;
    if(Xmax>XDmax) continue;
    if(Ymin<YDmin) continue;
    if(Ymax>YDmax) continue;
//===== отсечение контуров  по местоположению =====
    if(maxXY.GetX()>32767||maxXY.GetY()>32767) continue;
//===== отсечение помеченных контуров =====
    if(NumObj<0) continue;
//===== отсечение по количеству узлов в метрике объекта =====
    if(Lobj>2048)    Nln=1;
//================ отсечение по габаритам ===============
//============= обработка большеформатных объектов =============
    if(Xmax-Xmin>Pmax||Ymax-Ymin>Pmax||Nln==1)
      {  
		  continue;
      }

    if(Xmax-Xmin<Pmin&&Ymax-Ymin<Pmin)
      { KPObj++; SMPObj+=Lobj;
      //if(PaspKmp[AP+11l]>0)
	    if( wContour->GetTypeContour() == 0)
        { if(Xmax-Xmin<OutMin&&Ymax-Ymin<OutMin) continue;
        }
      else
        { if(Xmax-Xmin<InMin&&Ymax-Ymin<InMin) continue;
        }
      Nsm=1;
      KSmObj++; SMSmObj+=Lobj;

      }
    else { SMDObj+=Lobj;}

    Thr=(Xmax-Xmin)/50;
    j=(Ymax-Ymin)/50;
    if(Thr<j) Thr=j;
    if(Thr>ThrDO)  Thr=ThrDO;
    Xc=(Xmax+Xmin)*StepX/2.;
    Yc=(Ymax+Ymin)*StepY/2.;


	//MoveObj=(PaspKmp[AP+1l]-1)*4096l+(PaspKmp[AP+2l]-1);
    //for(j=0; j<Lobj; j++) ZMo[j]=MetrKmp[MoveObj+j];
	for(j=0; j<n; j++)
	{
		wPoint = wContour->GeT_Point(j);
		ZMo[j*2] = wPoint.GetX();
		ZMo[j*2 + 1] = wPoint.GetY();
	}

//===== формирование описания объекта ========
	// < 0 - внутренний контур
    //if(PaspKmp[AP+11l]<0)
	if(wContour->GetTypeContour() == 1)
      {
      for(j=0; j<n; j++)
        { i=j+j;
	    // zmo - точки текущего контура
        Xi[n-j-1]=ZMo[i];
        Yi[n-j-1]=ZMo[i+1];
        }
      }
    else
      {
      for(j=0; j<n; j++)
        { i=j+j;
        Xi[j]=ZMo[i];
        Yi[j]=ZMo[i+1];
        }
      }
    Xi[n]=Xi[0];
    Yi[n]=Yi[0];
    n++;
	// целочисленное сжатие
    N=CompressInteger(Xi,Yi,n,Thr,Xi,Yi);
    Xmin=Xi[0]; Ymin=Yi[0];
    Xmax=Xi[0]; Ymax=Yi[0];
    for(j=1; j<N; j++)
      {
      if(Xmin>Xi[j]) Xmin=Xi[j];
      if(Xmax<Xi[j]) Xmax=Xi[j];
      if(Ymin>Yi[j]) Ymin=Yi[j];
      if(Ymax<Yi[j]) Ymax=Yi[j];
      }
    for(j=0; j<N; j++) {Xi[j]-=Xmin;  Yi[j]-=Ymin;}
//=================== обработка мелких объектов ===================
    if(Nsm==1)
      {
		  continue;
      }
//===================================================================
    So=0.;
    Xd[0]=Xp=Xi[0]*StepX;
    Yd[0]=Yp=Yi[0]*StepY;
    for(j=1; j<N; j++)
      {
      Xd[j]=Xt=Xi[j]*StepX;
      Yd[j]=Yt=Yi[j]*StepY;
      Sd[j-1]=St=sqrt((Xt-Xp)*(Xt-Xp)+(Yt-Yp)*(Yt-Yp));
      So+=St;
      Xp=Xt; Yp=Yt;
      }
    Sd[N-1]=0.;
//=== Вычисление статистических характеристик объекта ===
    MX=0.;  MY=0.;
    MX2=0.; MY2=0.;
    Xp=Xd[0]; Yp=Yd[0];
    for(j=1; j<N; j++)
      {
      Xt=Xd[j]; Yt=Yd[j]; St=Sd[j-1];
      Xa=Xt+Xp; Ya=Yt+Yp;
      MX+=Xa*St; MY+=Ya*St;
      MX2+=(Xa*Xa-Xt*Xp)*St; MY2+=(Ya*Ya-Yt*Yp)*St;
      Xp=Xt; Yp=Yt;
      }
    MX/=(2.*So); MY/=(2.*So);
    MX2/=(3.*So); MY2/=(3.*So);
    DX=MX2-MX*MX; DY=MY2-MY*MY;
    DW=DX+DY;
//==== Центрирование и масштабирование описания объекта ====
    Sq=SE/So;
    for(j=0; j<N; j++)
      {
      Xd[j]-=MX; Yd[j]-=MY;
      Xd[j]*=Sq; Yd[j]*=Sq; Sd[j]*=Sq;
      }
    DWo=DW*Sq*Sq;
//============ вычисление центральной точки ===============
    Xc=Xmin*StepX+MX;
    Yc=Ymin*StepY+MY;
    Wd=(Xmax-Xmin)*StepX;
    Hg=(Ymax-Ymin)*StepY;
//==== Вычисление минимальной невязки ====
    Em=DWb;  Ds=SE/Ny;  Sp=-Ds; St=Sd[0]; Np=1;
    for(i=0; i<Ny; i++)
      { Sp+=Ds;
      for(;St<Sp;) {Sp-=St; St=Sd[Np]; Np++; }
      k=CurrentDecribeObject(Xd,Yd,Sd,N-1,Np,Sp,SE,Xr,Yr,Sr);
      Nt=AuxiliaryDescriptionShapes(Nb,Xb,Yb,Sb,Xr,Yr,Sr,Xbv,Ybv,
                                   Xv,Yv,Sv);
      ErrorMatching(Xbv,Ybv,Xv,Yv,Sv,SE,Nt,&Sn,&Cs,&Re);
      Emin[i]=DWb-Re/DWo;
      }
    j=0;  k=0;
    if(Emin[0]<Emin[Ny-1]&&Emin[0]<=Emin[1])
      { Qmin[j]=0.; j++;
      }
    if(Emin[0]>Emin[Ny-1]&&Emin[0]>=Emin[1])
      { Qmax[k]=0.; k++;
      }
    for(i=1; i<Ny-1; i++)
      {
      if(Emin[i]<Emin[i-1]&&Emin[i]<=Emin[i+1])
        { Qmin[j]=i*Ds; j++;
        }
      if(Emin[i]>Emin[i-1]&&Emin[i]>=Emin[i+1])
        { Qmax[k]=i*Ds; k++;
        }
      }
    if(Emin[Ny-1]<Emin[Ny-2]&&Emin[Ny-1]<=Emin[0])
      { Qmin[j]=(Ny-1)*Ds; j++;
      }
    if(Emin[Ny-1]>Emin[Ny-2]&&Emin[Ny-1]>=Emin[0])
      { Qmax[k]=(Ny-1)*Ds; k++;
      }
    Nmin=j; Nmax=k;


    Rm=0.;
    for(i=0; i<*NumEtal; i++)
      {
      Ne=ParmEtals[i].n;
      NbP=ParmEtals[i].NbP;
      DE=ParmEtals[i].DE;
      if(emax<ParmEtals[i].EM) Emax=ParmEtals[i].EM;
      else Emax=emax;
      Eq=DE;
      Rq=1./(1.+Eq/2.5);


	  // расчет номера точки и величины сдвига
      for(j=0; j<Nmin; j++)
        { 
			Qt=Qmin[j]-ParmEtals[i].Qm;
        if(Qt<0.) Qt+=SE;
        Sp=Qt;
        ol=0;
        for(k=0; k<5; k++)
          {
          if(Sp<0.) Sp+=SE;
          else  if(Sp>SE) Sp-=SE;
          Np=0;                    // номер точки
          for(jt=0; jt<N; jt++)
            { Np=jt+1;
            if(Sd[jt]>Sp) break;
            Sp-=Sd[jt];            // расстояние от точки до текущего описания
            }
          CurrentDecribeObject(Xd,Yd,Sd,N-1,Np,Sp,SE,Xr,Yr,Sr);
          Nt=AuxiliaryDescriptionShapes(Ne,&Xe[NbP],&Ye[NbP],&Se[NbP],
                                    Xr,Yr,Sr,Xbv,Ybv,Xv,Yv,Sv);
          ErrorMatching(Xbv,Ybv,Xv,Yv,Sv,SE,Nt,&Sn,&Cs,&Re);
          // ошибка
		  Er=DE-Re/DWo;
          if(Eq>Er)  { Eq=Er; Snq=Sn; Csq=Cs; Rq=1./(1.+Eq/2.5);}
          // метод парабол
		  switch(k)
            {
            case 0: e1=Er; s1=sc=Qt; s0=Sp=sc-step; if(e1>Emax) ol=1;
                    break;
            case 1: e0=Er; s2=Sp=sc+step;
                    break;
            case 2: e2=Er;
                    if(e2+e0<=2.*e1)
                      {
                      if(e2>e0) Sp=s0-1.5;
                      else Sp=s2+1.5;
                      }
                    else
                      {
                      Sp=s1+(e0-e2)/(e0+e2-2.*e1)*step/2.;
                      if(e0>e2) { if(Sp>s2+1.5) Sp=s2+1.5;}
                      else {if(Sp<s0-1.5) Sp=s0-1.5;}
                      }
                    sc=Sp;
                    break;
            case 3: if(sc>s1)
                      { e0=e1; s0=s1;
                      if(sc<s2) { e1=Er; s1=sc;}
                      else { e1=e2; s1=s2; e2=Er; s2=sc;}
                      }
                    else
                      { e2=e1; s2=s1;
                      if(sc>s0) { e1=Er; s1=sc;}
                      else { e1=e0; s1=s0; e0=Er; s0=sc;}
                      }
                    if(e2*(s1-s0)+e0*(s2-s1)<=e1*(s2-s0))
                      {
                      if(e2>e0) Sp=s0-1.5;
                      else Sp=s2+1.5;
                      }
                    else
                      {
                      w1=(e0-e1)*(s2-s1); w2=(e2-e1)*(s1-s0);
                      Sp=s1+(w1*(s2-s1)-w2*(s1-s0))/(w1+w2)/2.;
                      if(e0>e2) { if(Sp>s2+1.5) Sp=s2+1.5;}
                      else {if(Sp<s0-1.5) Sp=s0-1.5;}
                      }
                    break;
            case 4: break;
            default: break;
            }
          if(ol==1) break;
          }

        if(Er<Eq)
          { Eq=Er; Snq=Sn; Csq=Cs;
          Rq=1./(1.+Eq/2.5);
          }
        }
      if(Rq>Rm)
        { Rm=Rq;  Snm=Snq; Csm=Csq;
        NEtal=i+1;
        }
      }

      { struct DataRecObj  RecObj;
      double Rat,Km,Xn,Yn,Xk,Yk,z1,z2;

	 
      if(Rm>SimPor)
        {
//======= формирование данных для передачи ========
        RecObj.NE=NEtal;
        Rat=sqrt(Snm*Snm+Csm*Csm);
        Snm/=Rat; Csm/=Rat;
        Km=Rat/DWo;
        Xn=ParmEtals[NEtal-1].Xn*ParmEtals[NEtal-1].Sq/(Sq*Km);
        Yn=ParmEtals[NEtal-1].Yn*ParmEtals[NEtal-1].Sq/(Sq*Km);
        Xk=ParmEtals[NEtal-1].Xk*ParmEtals[NEtal-1].Sq/(Sq*Km);
        Yk=ParmEtals[NEtal-1].Yk*ParmEtals[NEtal-1].Sq/(Sq*Km);
        if(wContour->GetTypeContour() == 1) RecObj.Inet=1;
        else RecObj.Inet=0;

        z1=Xn*Csm-Yn*Snm+Xc;
        z2=Xk*Csm-Yk*Snm+Xc;
        RecObj.Xn=z1/StepX+0.5;
        z2-=z1;
        if(z2<0) { RecObj.SgnDx=1; z2=-z2;}
        else RecObj.SgnDx=0;
        RecObj.dx=z2/StepX+0.5;

        z1=Xn*Snm+Yn*Csm+Yc;
        z2=Xk*Snm+Yk*Csm+Yc;
        RecObj.Yn=z1/StepY+0.5;
        z2-=z1;
        if(z2<0) { RecObj.SgnDy=1; z2=-z2;}
        else RecObj.SgnDy=0;
        RecObj.dy=z2/StepY+0.5;
        }
      else
        { //======= запись нового эталона и объекта ==========
			continue;
        }

      DiscreteObjects[KDObj]=RecObj;
      KDObj++;
      NObj+=(2*ParmEtals[RecObj.NE-1].n);

      }
    }
  }



//---------------------------------------------------------------------------
// главная функция распознавания
int AutoRecognition (TScene *SceneIn, double stepxin, double stepyin)
  {
  int Lp=16,NumEtal;
  long Mn,Mk,KObj;
  long NumPoints;
  double *Xe,*Ye,*Se;
  union { long ll; float fl; unsigned short Sh[2];} ul;
  struct ParmEtal *ParmEtals;
  long KPMetr,SNE,SPE,NbPi;
  T_Point maxXY;
  T_Point minXY;
  TContour *wContour;
  int NumberFileEtalon;
  FILE *pFileEtalons;
  double StepX,StepY;

  

  // описание базового контура (треугольник Пифагора)
  ul.ll=0xBFD55555;  Xb[0]=Xb[1]=Xb[3]=ul.fl;
  ul.ll=0x40555555;  Xb[2]=ul.fl;
  ul.ll=0xC0855555;  Yb[0]=Yb[3]=ul.fl;
  ul.ll=0x40200000;  Yb[1]=Yb[2]=ul.fl;
  ul.ll=0x40D55555;  Sb[0]=ul.fl;
  ul.ll=0x40A00000;  Sb[1]=ul.fl;
  ul.ll=0x41055555;  Sb[2]=ul.fl;
  ul.ll=0x40F471C6;  DWb=ul.fl;
  

  Nb=4; SE=20.;
  SimPor=0.95;
  Pmax=1000; Pmin=5;    // габариты
  StepX=25.4/stepxin;
  StepY=25.4/stepyin;
  //StepX=stepxin;
  //StepY=stepyin;
  //StepX=25.4/300.; StepY=25.4/300.;
  OutMin=InMin=Pmin;
  SelectData=false;
  ThrDO=1;
  ThrBO=2;
  AdrBig=1;
  NObj=0;

  pFileEtalons = fopen("E:\\PROJECT\\NavigationSystem\\3D-Navigation\\A\\table\\mark_5.etl", "rb");


  KObj = SceneIn->GetNumObjects();
  Mn=1; Mk=KObj;
  // Xmin, Xmax, Ymin, Ymax растра
  //XDmin=CapKmp[15]; XDmax=CapKmp[16];
  //YDmin=CapKmp[17]; YDmax=CapKmp[18];
  maxXY = SceneIn->FindMaxXY();
  minXY = SceneIn->FindMinXY();
  XDmin=minXY.GetX(); XDmax = maxXY.GetX();
  YDmin=minXY.GetY(); YDmax = maxXY.GetY();

  KPObj=0; SMPObj=0;  KBObj=0; SMBObj=0;
  KDObj=0; SMDObj=0;  KSmObj=0; SMSmObj=0;
  KSmObj3=0; KSmObj4=0; KSmObj5=0; N_Dmax=0;
  
  long i,AP,Dx,Dy,N,G;
  for(i=0; i<Mk; i++)
    {
    AP=i*Lp;
    //Dx=PaspKmp[AP+6l]-PaspKmp[AP+5l];
    //Dy=PaspKmp[AP+8l]-PaspKmp[AP+7l];
	wContour = SceneIn->GetContour(i);
	maxXY = wContour->FindMaxXY();
	minXY = wContour->FindMinXY();
	Dx=maxXY.GetX() - minXY.GetX();
	Dy=maxXY.GetY() - minXY.GetY();
    //N=PaspKmp[AP+3l];
	N=wContour->GetNumPoints();
	// выбираем контур с макс числом точек
    if(N_Dmax<N) N_Dmax=N;
	// G - габарит
    if(Dx>Dy) G=Dx;
    else G=Dy;
    if(G<Pmin)
      { KSmObj++;  SMSmObj+=N;
      if(G<3) KSmObj3++;
      else
        {
        if(G<4) KSmObj4++;
        else KSmObj5++;
        }
      }
    else
      {
      if(G>Pmax||N>2048) { KBObj++; SMBObj+=N; }
      else  { KDObj++; SMDObj+=N; }
      }
    }
  // KDObj - количество дискретных объектов
  // KBObj - количество больших объектов

  NumberFileEtalon = fileno(pFileEtalons);
  SizeEtalons=filelength(NumberFileEtalon);
  NumPoints=SizeEtalons/sizeof(unsigned int);
  SMEtal=NumPoints;
  KPMetr=NumPoints; //+SMDObj/2+KDObj;
  ZME= new unsigned int[KPMetr];
  if(ZME==NULL)
    {
    if(NumParm>0) ErrorCode = 1;
    else
      {
		  strcpy(text,"Нет памяти под эталоны");
      }
    goto END;
    }

  NumEtal=0;  NbPk=0; NbPi=0;
  KPMetr=NumPoints;
  if(NumPoints>0)
    {

    if(fread(ZME, sizeof(int), KPMetr, pFileEtalons) != SizeEtalons / sizeof(unsigned int))
      {
      if(NumParm>0) ErrorCode = 1;
      else
        {
        strcpy(text,"Ошибка чтения метрики эталонов");
        }
      goto END;
      }

    CEM=(struct CapEtalons *)ZME;
    for(;;)
      {
      NbPk+=CEM[NbPi].N;
      NbPi+=CEM[NbPi].N+1;
      NumEtal++;
      if(NbPi>=NumPoints) break;
      }
    }


//======== отведение памяти под объекты и эталоны ===========
  SNE=NumEtal+KDObj;
  SPE=KPMetr; //+SMDObj+KDObj+KDObj;
  ParmEtals=new struct ParmEtal[SNE];
  if(ParmEtals==NULL)
    {
    if(NumParm>0) ErrorCode = 1;
    else
      {
      strcpy(text,"Нет памяти под параметры эталонов");
      }
    goto END;
    }
  Xe=new double[SPE]; Ye=new double[SPE]; Se=new double[SPE];
  if(Xe==NULL||Ye==NULL||Se==NULL)
    {
    if(NumParm>0) ErrorCode = 1;
    else
      {
      strcpy(text,"Нет памяти под метрику эталонов");
      }
    goto END;
    }
  DiscreteObjects=new struct DataRecObj[KDObj+1];
  //BigObjects=new unsigned int[SMBObj*5/4+4*KBObj+1];
  //XYObj=new long[N_Dmax];

  ReadEtalons(NumEtal,ParmEtals,StepX,StepY,SE,Xe,Ye,Se);

  RecognitionObject(SceneIn, KObj,&NumEtal,ParmEtals,StepX,StepY,SE,Xe,Ye,Se);



END:
  if(BigObjects!=NULL) delete BigObjects;
  if(DiscreteObjects!=NULL) delete DiscreteObjects;
  fclose(pFileEtalons);  if(SMEtal<1) unlink(NameFileEtalons);
  if(ZME!=NULL) delete ZME;
  if(NumSmall5!=NULL) delete NumSmall5;
  if(XYObj!=NULL) delete XYObj;
  if(ParmEtals!=NULL) delete ParmEtals;
  if(Xe!=NULL) delete Xe;
  if(Ye!=NULL) delete Ye;
  if(Se!=NULL) delete Se;
  if(AllKmp!=NULL) delete AllKmp;


  return 0;
  }
//---------------------------------------------------------------------------
