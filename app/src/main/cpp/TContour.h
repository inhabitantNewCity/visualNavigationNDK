#ifndef TContour_HH
#define TContour_HH
#pragma once
#include <stdio.h>
#include <vector>
#include "TPoint.h"
#include "data.h"

using namespace std;


class TContour {

private:
	vector <T_Point> *VectorPoints;              // ������ �����
	int NumPoints;                              // ���������� ����� � �������
	char TypeContour;                           // ��� ������� (0 - ������� ������, 1 - ���������� ������)
	T_Point minXY;                               // ����������� �����
	T_Point maxXY;                               // ������������ �����
	short Class;                                // �����

public:
	TContour ();
	TContour (char typecontour);
	~TContour ();          
	TContour (const TContour &bj);
	// ���������� ����� � ������
	void AddPoint (T_Point *point);
	void AddPoint (int x, int y);
	// �������� �� ���������� ����� � �������
	bool CheckPointInContour (T_Point *point);
	// �������� �����
	int DeletePoint (int x, int y);
	// �������� �� ������� �������
	bool isEmpty ();
	// ����������, �������� �� ������ ������� 
	bool isExternalContour ();
	// ���������� ���������� ����� � �������
	int GetNumPoints ();
	// ���������� ����� ������� �� ���������� x
	int GeT_PointsX (int *pX);
	// ���������� ����� ������� �� ���������� y
	int GeT_PointsY (int *pY);
	// ���������� ����� ������� �� ���������� x � y
	int GeT_PointsXY (int *pXY);
	// ���������� �����
	T_Point GeT_Point (int nPoint);
	// ���������� ��� �������
	char GetTypeContour ();
	// ������� � �������������� ����������
	void ToMathCoord (int height);
	// ���������� min �� ���������� x � y
	T_Point FindMinXY ();
	// ���������� max �� ���������� x � y
	T_Point FindMaxXY ();
	// ���������� � ������� ������� ���������
	void ToLocalSysCoord ();
	// ������ ������� � ���� (��������)
	int WriteContour (FILE *pHeaderFile);
};
#endif
