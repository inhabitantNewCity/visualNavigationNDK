//========= обращение матрицы   ===========
//---------------------------------------------------------------------------
#pragma warn -pck
#pragma option -a1
//---------------------------------------------------------------------------
#include <math.h>
void longMinv(long double *a, int n, long double *det, int *l, int *m)
  {
//------------ входные данные ---------
// a - массив элементов матрицы рамером nхn , построчный. ћатрица разрушаетс€,
//      на ее месте будут записаны элементы обратной матрицы.
//  n - размерность матрицы
//  l, m - рабочие массивы размерами   n
//------------ выходные данные ---------
// det -определитель матрицы
// a - массив элементов обратной матрицы
 
  int nk,k,i,j,kk,iz,ij,ki,ji,jp,jk,ik,kj,jq,jr;
  long double biga,hold,d;
  d=1.0; nk=-n; k=-1;
  while(++k<n)
    { nk+=n; l[k]=k+1; m[k]=k+1; kk=nk+k; biga=a[kk]; j=k-1;
    while(++j<n)
      { iz=n*j; i=k-1;
      while(++i<n)
        { ij=iz+i;
        if(fabsl(biga)<fabsl(a[ij]))
          { biga=a[ij]; l[k]=i+1; m[k]=j+1;
          }
        }
      }
      j=l[k];
      if(j>k+1)
        { ki=k-n; i=-1;
        while(++i<n)
          { ki+=n; hold=-a[ki]; ji=ki-k+j-1; a[ki]=a[ji]; a[ji]=hold;
          }
        }
      i=m[k];
      if(i>k+1)
        { jp=n*(i-1); j=-1;
        while(++j<n)
          { jk=nk+j; ji=jp+j; hold=-a[jk]; a[jk]=a[ji]; a[ji]=hold;
          }
        }
      if(biga==0.) {d=0.0; *det=d; return;}
      i=-1;
      while(++i<n) { if(i!=k) {ik=nk+i; a[ik]/=(-biga);}}
      i=-1;
      while(++i<n)
        {ik=nk+i; hold=a[ik]; ij=i-n; j=-1;
        while(++j<n)
          { ij+=n; if(i!=k) { if(j!=k) { kj=ij-i+k; a[ij]+=hold*a[kj];}}
          }
        }
      kj=k-n; j=-1;
      while(++j<n)
        { kj+=n; if(j!=k) a[kj]/=biga;
        }
      d*=biga; a[kk]=1.0/biga;
    }
  k=n;
  while(--k>0)
    { i=l[k-1];
    if(i>k)
      {jq=n*(k-1); jr=n*(i-1); j=-1;
      while(++j<n)
        { jk=jq+j; hold=a[jk]; ji=jr+j; a[jk]=-a[ji]; a[ji]=hold;
        }
      }
    j=m[k-1];
    if(j>k)
      { ki=k-n; i=-1;
      while(++i<n)
        { ki+=n; hold=a[ki-1]; ji=ki-k+j; a[ki-1]=-a[ji-1]; a[ji-1]=hold;
        }
      }
    }
  *det=d;
  return;
  }