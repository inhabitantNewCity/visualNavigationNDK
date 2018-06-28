#include <stdlib.h>
#include "TContour.h"
#include "TShtrihPicture.h"
#include "TShtrih.h"
#include "ContourModel.h"
#include "TScene.h"

// главная функция прослеживания контуров
void FindListContour (TShtrihPicture *shtrpic, TScene *SceneIn)
{
	//list <TContour> *ListContour;
	list<TShtrih>::iterator itercurshtr;

    //ListContour = new list <TContour> ();

	for (int vi = 0; vi < shtrpic->Height; vi ++)
	{
		int vk = 0;
		for (itercurshtr = shtrpic->MasShtrih[vi].begin (); itercurshtr != shtrpic->MasShtrih[vi].end (); itercurshtr ++)
		{
			// если снизу или сверху есть штрихи
			if ((itercurshtr->SWs > 0) || (itercurshtr->SWp > 0))
			{   
				// если координаты начала штриха не встречаются ни в одном из контуров
				if (!(CheckPointInListContour (SceneIn, itercurshtr->begin, vi)))
				{
					// создание контура с пометкой "внешний"
					TContour *AddContour = new TContour ('0');
					// прослеживание внешних контуров
					FindExternalContour (shtrpic, AddContour, vi, vk);
					// добавление найденного контура в список контуров
					SceneIn->AddContour(AddContour);
				}
				// если координаты конца штриха не встречаются ни в одном из контуров
				if (!(CheckPointInListContour (SceneIn, itercurshtr->end, vi)))
				{
					// создание контура с пометкой "внутренний"
					TContour *AddContour = new TContour ('1');
				    // прослеживание внутренних контуров
				    FindInternalContour (shtrpic, AddContour, vi, vk);
					// добавление найденного контура в список контуров
				    SceneIn->AddContour(AddContour);
				}
			}
			vk ++;
		}


	}
}

// проверка на нахождение точки в списке контуров
bool CheckPointInListContour (TScene *wSceneIn, int x, int y)
{
	list <TContour>::iterator itercontour;
	T_Point checkpoint (x, y);
    int NumContours = wSceneIn->GetNumObjects ();
	TContour *wContour;

	for (int vi = 0; vi < NumContours; vi ++)
	{
		wContour = wSceneIn->GetContour(vi);
		if (wContour->CheckPointInContour (&checkpoint))
			return true;
	}

	return false;
}
// прослеживание внешнего контура
void FindExternalContour (TShtrihPicture *shtrpic, TContour *extcontour, int NumStr, int NumShtr)
{
	// -----!обход внешнего контура осуществляется против часовой стрелки!----- 

	TShtrih *CurShtr = shtrpic->GetShtrFromList (NumStr, NumShtr);    // текущий штрих
	T_Point StartPoint (CurShtr->begin, CurShtr->Nstr);                // стартовая точка контура
	T_Point AddNextPoint (-1, -1);                                     // следующая добавляемая точка

    // добавление точки в контур
	extcontour->AddPoint (CurShtr->begin, CurShtr->Nstr);
	int prev, next;            // куда был, будет совершен переход

	prev = 1;
	next = 1;
    
	// определение направления во внешнем контуре
	int direction = shtrpic->GetDirectionInExtContour (CurShtr, prev, next);

		switch (direction)
		{
			case 0:
				{
					CurShtr = shtrpic->GetBottomLeftShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 1:
				{
					extcontour->AddPoint (&AddNextPoint);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 2:
				{
					CurShtr = shtrpic->GetTopRightShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 3:
				{
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 4:
				{
					CurShtr = shtrpic->GetNextShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 5:
				{
					CurShtr = shtrpic->GetPrevShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
		}

		prev = next;

	while (StartPoint != AddNextPoint)
	{
		int direction = shtrpic->GetDirectionInExtContour (CurShtr, prev, next); 

		switch (direction)
		{
			//   0-------------0
            //   |
			//  \|/
			//   0-------------0
			case 0:
				{
                    extcontour->AddPoint (&AddNextPoint);
					// текущий штрих становится нижним левым
					CurShtr = shtrpic->GetBottomLeftShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   0-------------0
			//   0------------->
			case 1:
				{
					extcontour->AddPoint (&AddNextPoint);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   0-------------0
            //                /|\
			//                 |
			//   0-------------0
			case 2:
				{
					extcontour->AddPoint (&AddNextPoint);
					// текущий штрих становится верхним правым
					CurShtr = shtrpic->GetTopRightShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
            //   <-------------0
			//   0-------------0
			case 3:
				{
					extcontour->AddPoint (&AddNextPoint);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			// 0--------------------0
			// 0-------0===>0-------0 
			case 4:
				{
					extcontour->AddPoint (&AddNextPoint);
                    // текущий штрих становится следующим
					CurShtr = shtrpic->GetNextShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   0-----0<===0-----0  
		    // 0--------------------0
			case 5:
				{
					extcontour->AddPoint (&AddNextPoint);
					// текущий штрих становится предыдущим
					CurShtr = shtrpic->GetPrevShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
		}
		prev = next;
	}

}
// прослеживание внутреннего контура
void FindInternalContour (TShtrihPicture *shtrpic, TContour *intercontour, int NumStr, int NumShtr)
{
	// -----!обход внутреннего контура осуществляется по часовой стрелке!----- 

	TShtrih *CurShtr = shtrpic->GetShtrFromList (NumStr, NumShtr);   // текущий штрих
	T_Point StartPoint (CurShtr->end, CurShtr->Nstr);                 // стартовая точка контура
	T_Point AddNextPoint (-1, -1);                                    // следующая добавляемая точка 


	intercontour->AddPoint (CurShtr->end, CurShtr->Nstr);
	int prev, next;

	prev = 0;
	next = 0;
    // получение направления во внутреннем контуре
	int direction = shtrpic->GetDirectionInterContour (CurShtr, prev, next);

	switch (direction)
		{
			case 0:
				{
					CurShtr = shtrpic->GetBottomLeftShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 1:
				{	
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 2:
				{
					// текущий штрих становится верхним правым
					CurShtr = shtrpic->GetTopRightShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 3:
				{
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 4:
				{
					// текущий штрих становится следующим
					CurShtr = shtrpic->GetNextShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 5:
				{
					// текущий штрих становится предыдущим
					CurShtr = shtrpic->GetPrevShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
		}
	prev = next;

	while (StartPoint != AddNextPoint)
	{
		int direction = shtrpic->GetDirectionInterContour (CurShtr, prev, next); 

		switch (direction)
		{
			//   0-------------0
            //                 |
			//                \|/
			//   0-------------0
			case 0:
				{
					intercontour->AddPoint (&AddNextPoint);
					// текущий штрих становится нижним левым
					CurShtr = shtrpic->GetBottomLeftShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   0-------------0
			//   0------------->
			case 1:
				{
					intercontour->AddPoint (&AddNextPoint);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
            //   0-------------0
            //                /|\
			//                 |
			//   0-------------0
			case 2:
				{
					intercontour->AddPoint (&AddNextPoint);
					CurShtr = shtrpic->GetTopRightShtrih (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   <-------------0
			//   0-------------0
			case 3:
				{
					intercontour->AddPoint (&AddNextPoint);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			// 0--------------------0
			// 0-------0===>0-------0 
			case 4:
				{
					intercontour->AddPoint (&AddNextPoint);
					CurShtr = shtrpic->GetNextShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			//   0-----0<===0-----0  
		    // 0--------------------0
			case 5:
				{
					intercontour->AddPoint (&AddNextPoint);
					CurShtr = shtrpic->GetPrevShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
		}
		prev = next;
	}
}