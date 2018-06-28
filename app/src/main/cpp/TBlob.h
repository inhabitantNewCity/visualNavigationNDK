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
	list<TShtrih> *ListShtrih_Blob;         // список ПРО
    int NumShtrih_Blob;                    // количество ПРО в блобе

public:
    TBlob();
	TBlob (const TBlob &Obj);
    ~TBlob();
	// проверяет нужно ли добавлять штрих 
    bool CheckingToAddShtrih(TShtrih &addshtrih);
    // добавление списка штрихов в конец блоба
    void AddListShtrih(list<TShtrih> *addlistshtrih);
    // добавление в конец списка штрихов
    void AddShtrih(TShtrih &addshtrih);
    // создание списка блобов
    list<TBlob>* CreateListBLOB(TShtrihPicture &ShPicture);
	// отрисовка блоба
//    void DrawBlobOnPicture(Bitmap ^imagein, list <TShtrih> *ListShtrih, int zoom, Color color);
	// отрисовка линии на изображении imagein в строке nstr
//	void DrawLine(Bitmap ^imagein, int nstr, int start, int end, int width, Color color);   
	// поиск номера ПРО в списке, который содержит shtrih; возвращает -1, если не принадлежит ни одному ПРО из списка
    int SearchInList_PROWithShtrih(list<TPRO> *listpro, TShtrih shtrih);
	// создание списка ПРО и списка соединительных штрихов для изображения
    void CreateListPRO_and_ConnectShtr(list<TShtrih> *ListConnectShtrihs, list<TPRO> *ListPro, TShtrihPicture ShPicture);
	// очистка содержимого блоба
    void Clear();          
	// проверка на принадлежность штриха блобу
    bool Check_BlobHasShtrih(TShtrih &shtrih);

};
#endif
