#include <stdlib.h>
#include "TContour.h"
#include "TShtrihPicture.h"
#include "TShtrih.h"
#include "ContourModel.h"
#include "TScene.h"

// ������� ������� ������������� ��������
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
			// ���� ����� ��� ������ ���� ������
			if ((itercurshtr->SWs > 0) || (itercurshtr->SWp > 0))
			{   
				// ���� ���������� ������ ������ �� ����������� �� � ����� �� ��������
				if (!(CheckPointInListContour (SceneIn, itercurshtr->begin, vi)))
				{
					// �������� ������� � �������� "�������"
					TContour *AddContour = new TContour ('0');
					// ������������� ������� ��������
					FindExternalContour (shtrpic, AddContour, vi, vk);
					// ���������� ���������� ������� � ������ ��������
					SceneIn->AddContour(AddContour);
				}
				// ���� ���������� ����� ������ �� ����������� �� � ����� �� ��������
				if (!(CheckPointInListContour (SceneIn, itercurshtr->end, vi)))
				{
					// �������� ������� � �������� "����������"
					TContour *AddContour = new TContour ('1');
				    // ������������� ���������� ��������
				    FindInternalContour (shtrpic, AddContour, vi, vk);
					// ���������� ���������� ������� � ������ ��������
				    SceneIn->AddContour(AddContour);
				}
			}
			vk ++;
		}


	}
}

// �������� �� ���������� ����� � ������ ��������
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
// ������������� �������� �������
void FindExternalContour (TShtrihPicture *shtrpic, TContour *extcontour, int NumStr, int NumShtr)
{
	// -----!����� �������� ������� �������������� ������ ������� �������!----- 

	TShtrih *CurShtr = shtrpic->GetShtrFromList (NumStr, NumShtr);    // ������� �����
	T_Point StartPoint (CurShtr->begin, CurShtr->Nstr);                // ��������� ����� �������
	T_Point AddNextPoint (-1, -1);                                     // ��������� ����������� �����

    // ���������� ����� � ������
	extcontour->AddPoint (CurShtr->begin, CurShtr->Nstr);
	int prev, next;            // ���� ���, ����� �������� �������

	prev = 1;
	next = 1;
    
	// ����������� ����������� �� ������� �������
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
					// ������� ����� ���������� ������ �����
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
					// ������� ����� ���������� ������� ������
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
                    // ������� ����� ���������� ���������
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
					// ������� ����� ���������� ����������
					CurShtr = shtrpic->GetPrevShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->end);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
		}
		prev = next;
	}

}
// ������������� ����������� �������
void FindInternalContour (TShtrihPicture *shtrpic, TContour *intercontour, int NumStr, int NumShtr)
{
	// -----!����� ����������� ������� �������������� �� ������� �������!----- 

	TShtrih *CurShtr = shtrpic->GetShtrFromList (NumStr, NumShtr);   // ������� �����
	T_Point StartPoint (CurShtr->end, CurShtr->Nstr);                 // ��������� ����� �������
	T_Point AddNextPoint (-1, -1);                                    // ��������� ����������� ����� 


	intercontour->AddPoint (CurShtr->end, CurShtr->Nstr);
	int prev, next;

	prev = 0;
	next = 0;
    // ��������� ����������� �� ���������� �������
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
					// ������� ����� ���������� ������� ������
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
					// ������� ����� ���������� ���������
					CurShtr = shtrpic->GetNextShtr (CurShtr);
					AddNextPoint.SetX (CurShtr->begin);
					AddNextPoint.SetY (CurShtr->Nstr);
					break;
				}
			case 5:
				{
					// ������� ����� ���������� ����������
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
					// ������� ����� ���������� ������ �����
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