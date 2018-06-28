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
// ����������� �����������
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

// �������� � ����� ������
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
// ���������� �����
void TScene::DrawScene (Graphics ^gr, int Height)
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].DrawContour(gr, Height);
	}
}
// �������� ������� ��� �������
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
// ������� ������� �� ����� ��������
int TScene::ReadContours(char *pFileName)
{
	FILE *pHeaderEtalons;
	int filenoetl;
	unsigned int filelengthelt;                           // ���������� �����
	int NumPoints = 0;
	unsigned int *pMemEtalons;                            // ��������� ������ ��� �������� ��������
	int NumEtalons;                                       // ���������� ��������
	hEtalon *phEtalons;                                   // ��������� ��������
	int index = 1;
	int CurNumPoints = 0;                                 // ������� ���������� �����
	sPoint *pPointEtalons;                                 // ����� �������� 
	sPoint PointEtl;
	TContour *nContour;                                   

	
	pHeaderEtalons = fopen(pFileName, "rb");
	if (pHeaderEtalons == NULL) return -1;

	filenoetl = fileno(pHeaderEtalons);
	filelengthelt = filelength(filenoetl);

	
	NumPoints = filelengthelt/sizeof(int);

	// ��������� ������ ��� �������
	pMemEtalons = new unsigned int[NumPoints];
	// ������ ��������
	fread(pMemEtalons, sizeof(int), NumPoints, pHeaderEtalons);

	// ���������� ��������
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

// ��������� �������
TContour *TScene::GetContour (int nContour)
{
	if (nContour < 0 || nContour >= numObj || isEmpty ())
		return NULL;
	return &Contours[0][nContour];
}


// �������� �� ������� �����
bool TScene::isEmpty ()
{
	if(!numObj) return true;
	else return false;
}

// ���������� ���������� �������� � �����
int TScene::GetNumObjects ()
{
	return numObj;
}

// ���������� �������� ����� � �������������� ������� ���������
void TScene::toMathCoord (int height)
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].ToMathCoord(height);
	}
}
// ���������� ������� ������� � ��������� ������� ���������
void TScene::toLocalSysCoord ()
{
	for (int vi = 0; vi < numObj; vi ++)
	{
		Contours[0][vi].ToLocalSysCoord();
	}
}
// ����� max �� X � Y
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
// ����� min �� X � Y
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


