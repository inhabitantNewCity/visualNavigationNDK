#ifndef TBlob_HH
#define TBlob_HH
#pragma once
#include <stdio.h>
#include <list>
#include "TShtrih.h"
#include "TPRO.h"
#include "TShtrihPicture.h"

// using namespace System::Drawing;

using namespace std;

class TBlob {

public:
	list<TShtrih> *ListShtrih_Blob;         // ������ ���
    int NumShtrih_Blob;                    // ���������� ��� � �����

public:
    TBlob();
	TBlob (const TBlob &Obj);
    ~TBlob();
	// ��������� ����� �� ��������� ����� 
    bool CheckingToAddShtrih(TShtrih &addshtrih);
    // ���������� ������ ������� � ����� �����
    void AddListShtrih(list<TShtrih> *addlistshtrih);
    // ���������� � ����� ������ �������
    void AddShtrih(TShtrih &addshtrih);
    // �������� ������ ������
    list<TBlob>* CreateListBLOB(TShtrihPicture &ShPicture);
	// ��������� �����
//    void DrawBlobOnPicture(Bitmap ^imagein, list <TShtrih> *ListShtrih, int zoom, Color color);
	// ��������� ����� �� ����������� imagein � ������ nstr
//	void DrawLine(Bitmap ^imagein, int nstr, int start, int end, int width, Color color);   
	// ����� ������ ��� � ������, ������� �������� shtrih; ���������� -1, ���� �� ����������� �� ������ ��� �� ������
    int SearchInList_PROWithShtrih(list<TPRO> *listpro, TShtrih shtrih);
	// �������� ������ ��� � ������ �������������� ������� ��� �����������
    void CreateListPRO_and_ConnectShtr(list<TShtrih> *ListConnectShtrihs, list<TPRO> *ListPro, TShtrihPicture ShPicture);
	// ������� ����������� �����
    void Clear();          
	// �������� �� �������������� ������ �����
    bool Check_BlobHasShtrih(TShtrih &shtrih);

};
#endif
