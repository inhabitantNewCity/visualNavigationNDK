#pragma once
#include <list>
#include <fstream>
#include "TShtrih.h"

using namespace std;

class TShtrihPicture {

public:
    list<TShtrih> *MasShtrih;
    int Width;
    int Height;

public:
    TShtrihPicture();
	TShtrihPicture(int width, int height);
	TShtrihPicture (const TShtrihPicture &Obj);

    ~TShtrihPicture();
    void ClearShtihPicture();
	// формирование штрихов
    int ConvertToShtrih (unsigned char* mas);
	// определение коэффициента связи штриха
    void Get_SWp_SWs_Shtrih(TShtrih &shtrih, int &SWs, int &SWp);
	// запись штрихов в файл
	void SaveShtrihToFile (char *name);
	// получение штриха из списка
    TShtrih * GetShtrFromList (int NumStr, int Index);
	// получение самого левого штриха в предыдущей строке
	TShtrih* GetTopLeftShtrih (TShtrih *pShtr);
    // получение самого левого штриха в следующей строке
	TShtrih* GetBottomLeftShtrih (TShtrih *pShtr);
	// получение самого правого штриха в предыдущей строке
	TShtrih* GetTopRightShtrih (TShtrih *pShtr);
	// получение самого правого штриха в следующей строке
	TShtrih* GetBottomRightShtrih (TShtrih *pShtr);
	// определение направления обхода во внешнем контуре
    int GetDirectionInExtContour (TShtrih *pShtr, int prev, int &next);
	// получение предыдущего штриха в строке
	TShtrih* GetPrevShtr (TShtrih *pShtr);
	// получение следующего штриха в строке
	TShtrih* GetNextShtr (TShtrih *pShtr);
	// определение направления обхода во внутреннем контуре
	int GetDirectionInterContour (TShtrih *pShtr, int prev, int &next);
};