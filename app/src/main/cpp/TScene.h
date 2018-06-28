#ifndef TScene_HH
#define TScene_HH
#pragma once
#include <stdio.h>
#include <vector>
#include "TContour.h"
// using namespace System::Drawing;
using namespace std;


// класс сцена
class TScene {

private:
	int numObj;                             // количество контуров
	vector <TContour> *Contours;              // список контуров

public:
	TScene ();
	// конструктор копирования
	TScene (const TScene &obj);
	~TScene();

	// добавить в сцену контур
	int AddContour (TContour *addContour);
	// отобразить сцену
  	void DrawScene (Graphics ^gr, int Height);
	// записать контуры как эталоны
	int WriteContours (char *pFileName);
	// считать контуры из файла эталонов
	int ReadContours (char *pFileName);
	// получение контура
	TContour *GetContour (int nContour);
	// проверка на пустоту сцены
	bool isEmpty ();
	// возвращает количетсво конутров в сцене
	int GetNumObjects ();
	// приведение объектов сцены к математической системе координат
	void toMathCoord (int height);
	// приведение каждого объекта к локальной системе координат
	void toLocalSysCoord ();
	// найти max по X и Y
	T_Point FindMaxXY ();
	// найти min по X и Y
	T_Point FindMinXY ();

};
#endif
