#include "TContour.h"


TContour::TContour ()
{
	VectorPoints = new vector <T_Point> ();
	NumPoints = 0;
	Class = 0;
}

TContour::TContour (char typecontour)
{
	VectorPoints = new vector <T_Point> ();
	NumPoints = 0;
	TypeContour = typecontour;
	Class = 0;
}
TContour::~TContour ()
{
	delete VectorPoints;
}                   
TContour::TContour (const TContour &obj)
{

	NumPoints = obj.NumPoints;
    T_Point minXY0 = obj.minXY;                               // минимальная точка
	T_Point maxXY0 = obj.maxXY;                               // максимальная точка = obj.minXY.GetX(); 
   //     minXY = obj.minXY;
   //     maxXY = obj.maxXY;
        minXY.SetX(minXY0.GetX());
        minXY.SetY(minXY0.GetY());
        maxXY.SetX(maxXY0.GetX());
        maxXY.SetY(maxXY0.GetY());
	VectorPoints = new vector <T_Point> ();

	for (int vi = 0; vi < NumPoints; vi ++)
	{
		VectorPoints->push_back(obj.VectorPoints[0][vi]);
	}
	TypeContour = obj.TypeContour;
	Class = obj.Class;

}
// добавление точки в контур
void TContour::AddPoint (T_Point *point)
{
	VectorPoints->push_back (*point);
	NumPoints ++;
}

void TContour::AddPoint (int x, int y)
{
	T_Point point (x, y);

	VectorPoints->push_back (point);
	NumPoints ++;
}

// проверка на нахождении точки в контуре
bool TContour::CheckPointInContour (T_Point *point)
{
	vector <T_Point>::iterator iterpoint;

	for (iterpoint = VectorPoints->begin (); iterpoint != VectorPoints->end (); iterpoint++)
	{
		if  (*point == *iterpoint)
			return true;
	}

	return false;
}
// удаление точки
int TContour::DeletePoint (int x, int y)
{
	if (!isEmpty ())
	{
		vector <T_Point>::iterator iterpoint;

		for (iterpoint = VectorPoints->begin (); iterpoint != VectorPoints->end (); iterpoint++)
		{
			if  ((x == (*iterpoint).GetX ()) && (y == (*iterpoint).GetY ()))
			{
			    iterpoint = VectorPoints->erase (iterpoint);
				NumPoints --;
				return 0;
			}		
		}

	}
	return -1;
}

// проверка на пустоту контура
bool TContour::isEmpty ()
{
	if (!NumPoints) return true;
	else return false;
}

// определяет, является ли контур внешним 
bool TContour::isExternalContour ()
{
	if (TypeContour == '0') return 1;
	else return 0; 
}

// возврвщает количество точек в контуре
int TContour::GetNumPoints ()
{
	return NumPoints;
}


// возврвщает точки контура по координате x
int TContour::GeT_PointsX (int *pX)
{
	vector <T_Point>::iterator iterpoint;
	if (NumPoints == 0 || pX == NULL) return -1;

	int vi = 0;
	for (iterpoint = VectorPoints->begin (); iterpoint != VectorPoints->end (); iterpoint++)
	{
		pX[vi] = iterpoint->GetX ();
		vi ++;
	}
	return 0;
}
// возврвщает точки контура по координате y
int TContour::GeT_PointsY (int *pY)
{
	vector <T_Point>::iterator iterpoint;
	if (NumPoints == 0 || pY == NULL) return -1;
	
	int vi = 0;
	for (iterpoint = VectorPoints->begin (); iterpoint != VectorPoints->end (); iterpoint++)
	{
		pY[vi] = iterpoint->GetY ();
		vi ++;
	}
	return 0;
}

// возвращает точки контура по координате x и y
int TContour::GeT_PointsXY (int *pXY)
{
	vector <T_Point>::iterator iterpoint;
	if (NumPoints == 0 || pXY == NULL) return -1;

	int vi = 0;
	for (iterpoint = VectorPoints->begin (); iterpoint != VectorPoints->end (); iterpoint++)
	{
		pXY[vi] = iterpoint->GetX ();
		pXY[vi+1] = iterpoint->GetY ();
		vi ++;
	}

	return 0;


}

// возвращает точку
T_Point TContour::GeT_Point (int nPoint)
{
	if (nPoint < 0 || nPoint >= NumPoints)
		return -1;

	return VectorPoints[0][nPoint];
}

// возвращает тип контура
char TContour::GetTypeContour ()
{
	return TypeContour;
}


// перевод в математические координаты
void TContour::ToMathCoord (int height)
{
	int newY = 0;
	if (!isEmpty ())
	{
		for (int vi = 0; vi < NumPoints; vi ++)
		{
			newY = height - VectorPoints[0][vi].GetY();
			VectorPoints[0][vi].SetY(newY);
		}
	}
}

// нахождение min по координате x и y
T_Point TContour::FindMinXY ()
{
	T_Point minxy;
	int minx = VectorPoints[0][0].GetX();
	int miny = VectorPoints[0][0].GetY();

	if (!isEmpty ())
	{
		for (int vi = 0; vi < NumPoints; vi ++)
		{
			// нахождение минимума
			if (VectorPoints[0][vi].GetX() <= minx)
				minx = VectorPoints[0][vi].GetX();
			if (VectorPoints[0][vi].GetY() <= miny)
				miny = VectorPoints[0][vi].GetY();
		}
	}
	minxy.SetX(minx);
	minxy.SetY(miny);
	minXY = minxy;

	return minxy;
}
// нахождение max по координате x и y
T_Point TContour::FindMaxXY ()
{
	T_Point maxxy;
	
	int maxx = VectorPoints[0][0].GetX();
	int maxy = VectorPoints[0][0].GetY();

	if (!isEmpty ())
	{
		for (int vi = 0; vi < NumPoints; vi ++)
		{
			// нахождение максимума
			if (VectorPoints[0][vi].GetX() >= maxx)
				maxx = VectorPoints[0][vi].GetX();
			if (VectorPoints[0][vi].GetY() >= maxy)
				maxy = VectorPoints[0][vi].GetY();
		}
	}
	maxxy.SetX(maxx);
	maxxy.SetY(maxy);
	maxXY = maxxy;

	return maxxy;
}

// приведение к местной системе координат
void TContour::ToLocalSysCoord ()
{
	T_Point minxy;

	if (!isEmpty ())
	{
		minxy = FindMinXY ();

		for (int vi = 0; vi < NumPoints; vi ++)
		{
			VectorPoints[0][vi] = VectorPoints[0][vi] - minxy;
		}
	}
}


// запись контура в файл (эталонов)
int TContour::WriteContour (FILE *pHeaderFile)
{
	unsigned int *pBlock;
	uXY headEtal;
	uXY pointXY;

	if (pHeaderFile == NULL || isEmpty())
		return -1;

	// выделение блока памяти для размещения заголовочной информации и данных
	// заголовок состоит:
	// 1) количество точек
	// 2) номер класса
	// данные состоят из точек (x,y): short
	pBlock = new unsigned int[NumPoints + 1];   // +1 для заголовчной информации

	// упаковка информации
	// блок данных
	headEtal.Ptr[0] = NumPoints;
	headEtal.Ptr[1] = Class;

	pBlock[0] = headEtal.XY;

	for (int vi = 0; vi < NumPoints; vi ++)
	{
		pointXY.Ptr[0] = VectorPoints[0][vi].GetX();
		pointXY.Ptr[1] = VectorPoints[0][vi].GetY();

		pBlock[vi + 1] = pointXY.XY;
	}

	fwrite (pBlock, sizeof(int) * (NumPoints + 1), 1, pHeaderFile);

	delete []pBlock;

	return 0;
}
// отрисовка контура
void TContour::DrawContour (Graphics ^gr, int Height)
{
	vector <T_Point>::iterator iterprev;
	vector <T_Point>::iterator iternext;
	vector <T_Point>::iterator iterend;

	int zoom = 1;
	Color color;

	if (TypeContour == '0') color = Color::Blue;
	else color = Color::Green;

	Pen ^WorkPen = gcnew Pen (color);

	iterend = VectorPoints->end ();
    iterend --;
	iternext = VectorPoints->begin ();

	for (iterprev = VectorPoints->begin (); iterprev != iterend; iterprev++)
	{
		iternext ++;
		gr->DrawLine (WorkPen, iterprev->GetX (), Height - iterprev->GetY (), iternext->GetX (), Height - iternext->GetY());
    }


	gr->DrawLine (WorkPen, iterend->GetX (), Height - iterend->GetY (), (VectorPoints->begin ())->GetX (), Height -(VectorPoints->begin ())->GetY());

	delete WorkPen;
}


