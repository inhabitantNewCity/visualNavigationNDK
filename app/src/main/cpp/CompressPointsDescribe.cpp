//======= сжатие целочисленного точечного описания объекта =======
//------------ входные данные ---------
// N - исходное количество точек в описании объекта
// Pin - исходное описание объекта (массив точек)
//------------ выходные данные ---------
// Pout - описание объекта после сжатия (массив точек)
// m - количество точек в описании объекта после сжатия
//---------------------------------------------------------------------------
#pragma warn -pck
#pragma option -a1
//---------------------------------------------------------------------------
int CompressPointsDescribe(long  *Pin, int N, int d, long *Pout)
  { 
  long xn,yn,xk,yk,xz,yz,xp,yp;
  int m,j,k,i,jt,jp,jk,j2,jp2,mp;
  long xt,yt,d2,e,s,st,rn,rk;
  __int64 e2,p;
  bool Ibreak;
  xn=Pin[0]; yn=Pin[1];
  Pout[0]=xn;   Pout[1]=yn;
  if(N==1) return(N);

  m=1; i=1; k=0; j=1;
  d2=(long)d*d;
  for(;j<N;)
    {
    j2=j+j;
    xk=Pin[j2++]; yk=Pin[j2];  Ibreak=false;
    if(k>0)
      {
      xp=xk-xn; yp=yk-yn;
      s=(xp*xp+yp*yp);
      p=(__int64)s; p*=d2;
      for(jt=0; jt<k; jt++)
        { jp=i+jt; jp2=jp+jp;
        xt=Pin[jp2++]-xn; yt=Pin[jp2]-yn;
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
      { i+=k; xn=xz; yn=yz;    mp=m+m;
      Pout[mp]=xz; Pout[mp+1]=yz; k=0; m++; j=i;
      }
    }
  mp=m+m; Pout[mp]=xz; Pout[mp+1]=yz;  m++;

  if(m==2)
    {
    if(Pout[0]==Pout[2]&&Pout[1]==Pout[3]) m=1;
    }

  return (m);
  }
