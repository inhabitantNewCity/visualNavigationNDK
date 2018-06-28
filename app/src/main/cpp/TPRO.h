#ifndef TPRO_HH
#define TPRO_HH
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "TShtrih.h"
#include <list>

using namespace std;

class TPRO {

	
public:
	list<TShtrih> *ListShtrih_PRO;          // список штрихов ПРО
    int NumShtrih;                         // количество штрихов в ПРО

public:
    TPRO();
	TPRO (const TPRO &obj);
    ~TPRO();
    // добавление списка штрихов в конец ПРО
	void AddListShtrih(list<TShtrih> *addlistshtrih);   
    // добавление в конец списка штрихов
	void AddShtrih(TShtrih &addshtrih);       
	// отрисовка ПРО
    //void DrawPROOnPicture(Bitmap imagein, int zoom, Color color;  
	// проверка принадлежит ли штрих ПРО
    bool CheckPRO_Shtrih(TShtrih &shtrih);  
    // очистка содержимого ПРО
	void Clear();

};
#endif
