#include "TScene.h"
#include "AutoRecognition.h"
#include <stdio.h>
#include <stdlib.h>
#include <io.h>


TScene::TScene ()
{
	numObj = 0;
	Contours = new vector<TContour>;
}

TScene::~TScene ()
{
	delete Contours;
}
// конструктор копирования
TScene::TScene (const TScene &obj)
{
	numObj = obj.numObj;
	Contours = new vector<TContour>;

	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours->push_back(obj.Contours[0][vi]);
	}

	Contours = obj.Contours;
}

// добавить в сцену контур
int TScene::AddContour (TContour *addContour)
{
	if (!addContour->isEmpty())
	{
		Contours->push_back(*addContour);
		numObj ++;
	}
	else 
	{
		return 1;
	}
	return 0;
}
// отобразить сцену
void TScene::DrawScene (Graphics ^gr, int Height)
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].DrawContour(gr, Height);
	}
}
// записать контуры как эталоны
int TScene::WriteContours (char *pFileName)
{
	FILE *pHeaderEtalons;

	pHeaderEtalons = fopen(pFileName, "wb");
	if (pHeaderEtalons == NULL) return -1;

	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].WriteContour(pHeaderEtalons);
	}

	fclose(pHeaderEtalons);

	return 0;
}
// считать контуры из файла эталонов
int TScene::ReadContours(char *pFileName)
{
	FILE *pHeaderEtalons;
	int filenoetl;
	unsigned int filelengthelt;                           // количество точек
	int NumPoints = 0;
	unsigned int *pMemEtalons;                            // временная память для хранения эталонов
	int NumEtalons;                                       // количество эталонов
	hEtalon *phEtalons;                                   // заголовки эталонов
	int index = 1;
	int CurNumPoints = 0;                                 // текущее количество точек
	sPoint *pPointEtalons;                                 // точки эталонов 
	sPoint PointEtl;
	TContour *nContour;                                   

	
	pHeaderEtalons = fopen(pFileName, "rb");
	if (pHeaderEtalons == NULL) return -1;

	filenoetl = fileno(pHeaderEtalons);
	filelengthelt = filelength(filenoetl);

	
	NumPoints = filelengthelt/sizeof(int);

	// выделение памяти под эталоны
	pMemEtalons = new unsigned int[NumPoints];
	// чтение эталонов
	fread(pMemEtalons, sizeof(int), NumPoints, pHeaderEtalons);

	// извлечение эталонов
	phEtalons = (hEtalon *) pMemEtalons;
	pPointEtalons = (sPoint *) pMemEtalons;

	while (CurNumPoints < NumPoints)
	{
		nContour = new TContour ('1');
		for (int vi = 0; vi < phEtalons[CurNumPoints].N; vi ++)
		{
			PointEtl = pPointEtalons [index + vi];
			nContour->AddPoint(PointEtl.X, PointEtl.Y);
		}

		AddContour(nContour);
		CurNumPoints += phEtalons[CurNumPoints].N + 1;
		index = CurNumPoints + 1;
		NumEtalons ++;
	}


	delete pMemEtalons;

	return 0;
}

// получение контура
TContour *TScene::GetContour (int nContour)
{
	if (nContour < 0 || nContour >= numObj || isEmpty ())
		return NULL;
	return &Contours[0][nContour];
}


// проверка на пустоту сцены
bool TScene::isEmpty ()
{
	if(!numObj) return true;
	else return false;
}

// возвращает количетсво конутров в сцене
int TScene::GetNumObjects ()
{
	return numObj;
}

// приведение объектов сцены к математической системе координат
void TScene::toMathCoord (int height)
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].ToMathCoord(height);
	}
}
// приведение каждого объекта к локальной системе координат
void TScene::toLocalSysCoord ()
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].ToLocalSysCoord();
	}
}
// найти max по X и Y
T_Point TScene::FindMaxXY ()
{
	T_Point maxXY;
	T_Point curPoint;

	if (!isEmpty())
	{
		maxXY = Contours[0][0].FindMaxXY();
		for (int vi = 1; vi < numObj; vi ++)
		{
			curPoint = Contours[0][vi].FindMaxXY();
			if (curPoint.GetX() >= maxXY.GetX()) maxXY.SetX(curPoint.GetX());
			if (curPoint.GetY() >= maxXY.GetY()) maxXY.SetY(curPoint.GetY());
		}
	}

	return maxXY;
}
// найти min по X и Y
T_Point TScene::FindMinXY ()
{
	T_Point minXY;
	T_Point curPoint;

	if (!isEmpty())
	{
		minXY = Contours[0][0].FindMinXY();
		for (int vi = 1; vi < numObj; vi ++)
		{
			curPoint = Contours[0][vi].FindMinXY();
			if (curPoint.GetX() <= minXY.GetX()) minXY.SetX(curPoint.GetX());
			if (curPoint.GetY() <= minXY.GetY()) minXY.SetY(curPoint.GetY());
		}
	}

	return minXY;
}


