#include "TShtrihPicture.h"

TShtrihPicture::TShtrihPicture()
{
    Width = 0;
    Height = 0;
}
TShtrihPicture::TShtrihPicture (const TShtrihPicture &Obj)
{
	
	Width = Obj.Width;
	Height = Obj.Height;
	MasShtrih = new list<TShtrih>[Height];

	for (int vi = 0; vi < Height; vi ++)
	{
		MasShtrih[vi] = Obj.MasShtrih[vi];
		
	}
    
}


TShtrihPicture::~TShtrihPicture()
{
	delete []MasShtrih;
}
TShtrihPicture::TShtrihPicture(int width, int height): Width (width), Height (height)
{
	//��������� ������
    MasShtrih = new list<TShtrih>[Height];
}
void TShtrihPicture::ClearShtihPicture()
{
	for (int i = 0; i < Height; i++)
    {
		MasShtrih[i].clear ();
    }
}
// ������������ �������
int TShtrihPicture::ConvertToShtrih (unsigned char* mas)        
{
    // ��������� ���� ��� ������� ������
    int SWp = 0;
    int SWs = 0;
	bool flag = 0; 
	TShtrih tmp;

	for(int y = 0; y < Height; y++)

		for(int x = 0; x < Width; x++)
		{
			// ��������� ���������� ������ ������
			if(( mas[x + y * Width] == 1) && (flag == 0)) 
			{
				tmp.Nstr = y;
				tmp.begin = x;
				flag = 1;
			}
			// ��������� �������� ����� ����� �������� ��������� � ����� ������
			if(( mas[x + y * Width] == 1) && (x == Width -1)) 
			{
				flag = 0;
				tmp.end = x;
				MasShtrih[y].push_back(tmp);
			}
			// ��������� ���������� ����� ������
			if(( mas[x + y * Width] == 0) && (flag == 1)) 
			{
				flag = 0;
				tmp.end = x-1;
				MasShtrih[y].push_back(tmp);
			}
		}


	list <TShtrih>::iterator itershtr;

	// ����������� ������������ ����� ������ � �������
    for (int i = 0; i < Height; i++)
    {
		for(itershtr = MasShtrih[i].begin(); itershtr != MasShtrih[i].end (); itershtr++)
        {
            Get_SWp_SWs_Shtrih(*itershtr, SWs, SWp);
            (*itershtr).SWp = SWp;
			(*itershtr).SWs = SWs;
        }
    } 


    return 0;
}
 // ����������� ������������ ����� ������
void TShtrihPicture::Get_SWp_SWs_Shtrih(TShtrih &shtrih, int &SWs, int &SWp)   
{
	list <TShtrih>::iterator itershtr;

	SWs = 0;
    SWp = 0;
	
    // ���� ����� ��������� � ������ ������
    if (shtrih.Nstr == 0)
    {
        SWp = 0;
    }
    else
    {
		for (itershtr = MasShtrih[shtrih.Nstr - 1].begin (); itershtr != MasShtrih[shtrih.Nstr - 1].end (); itershtr++)
        {
            if ((((*itershtr).begin <= shtrih.begin) && ((*itershtr).end >= shtrih.begin)) ||
                (((*itershtr).end >= shtrih.end) && ((*itershtr).begin <= shtrih.end)) ||
                (((*itershtr).begin >= shtrih.begin) && ((*itershtr).end <= shtrih.end)) ||
                (((*itershtr).begin <= shtrih.begin) && ((*itershtr).end >= shtrih.end)))
            {
                SWp++;
            }
            if ((*itershtr).begin > shtrih.end)
                break;
        }
    }
    // ���� ����� ��������� � ��������� ������
    if (shtrih.Nstr == (Height - 1))
    {
        SWs = 0;
    }
    else
    {
        for (itershtr = MasShtrih[shtrih.Nstr + 1].begin (); itershtr != MasShtrih[shtrih.Nstr + 1].end (); itershtr++)
        {
            if ((((*itershtr).begin <= shtrih.begin) && ((*itershtr).end >= shtrih.begin)) ||
                (((*itershtr).end >= shtrih.end) && ((*itershtr).begin <= shtrih.end)) ||
                (((*itershtr).begin >= shtrih.begin) && ((*itershtr).end <= shtrih.end)) ||
                (((*itershtr).begin <= shtrih.begin) && ((*itershtr).end >= shtrih.end)))
            {
                SWs++;
            }
            if ((*itershtr).begin > shtrih.end)
                break;
        }
    }
}

void TShtrihPicture::SaveShtrihToFile (char *name)
{

	TShtrih tmpLine;
	ofstream F;
	F.open(name,ios::out | ios::binary);//������� �������� ����
	F.write((char*)&Width,sizeof(int));//������ ������ ��������� �����������, sizeof - ���������� ������ ���� �������� ������������� � ����
	F.write((char*)&Height,sizeof(int));//������ ������ ��������� �����������
	for (int vi = 0; vi < Height; vi ++)
	{
		while( !MasShtrih[vi].empty())//���� �� �������� ������
		{
			tmpLine = MasShtrih[vi].front();
			MasShtrih[vi].pop_front();
			F.write ((char*)&tmpLine.Nstr, sizeof (tmpLine.Nstr));
			F.write ((char*)&tmpLine.begin, sizeof (tmpLine.Nstr));
			F.write ((char*)&tmpLine.end, sizeof (tmpLine.Nstr));	
		}
	}
	F.close();	
}

// ��������� ������ �� ������
TShtrih* TShtrihPicture::GetShtrFromList (int NumStr, int Index)
{
	
	if (Index >= MasShtrih[NumStr].size ()) return NULL;
	
	list <TShtrih>::iterator itershtr;

	itershtr = MasShtrih[NumStr].begin ();

	for (int vi = 0; vi < Index; vi ++)
	{
		itershtr ++;
	}

	return &(*itershtr);

}
// ��������� ������ ������ ������ � ���������� ������
TShtrih* TShtrihPicture::GetTopLeftShtrih (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;

	for (itershtr = MasShtrih[pShtr->Nstr - 1].begin (); itershtr != MasShtrih[pShtr->Nstr - 1].end (); itershtr++)
    {
        if ((((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->begin)) ||
            (((*itershtr).end >= pShtr->end) && ((*itershtr).begin <= pShtr->end)) ||
            (((*itershtr).begin >= pShtr->begin) && ((*itershtr).end <= pShtr->end)) ||
            (((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->end)))
        {
           return &(*itershtr);
        }
    }

	return NULL;
}
// ��������� ������ ������ ������ � ��������� ������
TShtrih* TShtrihPicture::GetBottomLeftShtrih (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;

	for (itershtr = MasShtrih[pShtr->Nstr + 1].begin (); itershtr != MasShtrih[pShtr->Nstr + 1].end (); itershtr++)
    {
        if ((((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->begin)) ||
            (((*itershtr).end >= pShtr->end) && ((*itershtr).begin <= pShtr->end)) ||
            (((*itershtr).begin >= pShtr->begin) && ((*itershtr).end <= pShtr->end)) ||
            (((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->end)))
        {
           return &(*itershtr);
        }
    }

	return NULL;
}

// ��������� ������ ������� ������ � ���������� ������
TShtrih* TShtrihPicture::GetTopRightShtrih (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;
	bool flag = false;

	for (itershtr = MasShtrih[pShtr->Nstr - 1].begin (); itershtr != MasShtrih[pShtr->Nstr - 1].end (); itershtr++)
    {
		if ((((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->begin)) ||
            (((*itershtr).end >= pShtr->end) && ((*itershtr).begin <= pShtr->end)) ||
            (((*itershtr).begin >= pShtr->begin) && ((*itershtr).end <= pShtr->end)) ||
            (((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->end)))
        {
			flag = true;
		}

		if ((*itershtr).begin > pShtr->end)
		{
			break;
		}
    }
	if (flag)
	{
		itershtr --;
		return &(*itershtr);
	}
	else 
	{
		return NULL;
	}
}

// ��������� ������ ������� ������ � ��������� ������
TShtrih* TShtrihPicture::GetBottomRightShtrih (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;
	bool flag = false;

	for (itershtr = MasShtrih[pShtr->Nstr + 1].begin (); itershtr != MasShtrih[pShtr->Nstr + 1].end (); itershtr++)
    {
		if ((((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->begin)) ||
            (((*itershtr).end >= pShtr->end) && ((*itershtr).begin <= pShtr->end)) ||
            (((*itershtr).begin >= pShtr->begin) && ((*itershtr).end <= pShtr->end)) ||
            (((*itershtr).begin <= pShtr->begin) && ((*itershtr).end >= pShtr->end)))
        {
			flag = true;
		}

		if ((*itershtr).begin > pShtr->end)
		{
			break;
		}
    }
	if (flag)
	{
		itershtr --;
		return &(*itershtr);
	}
	else 
	{
		return NULL;
	}
}
// ��������� ����������� ������ � ������
TShtrih* TShtrihPicture::GetPrevShtr (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;
	int count = 1;

	itershtr = MasShtrih[pShtr->Nstr].begin ();

	if (*pShtr == *itershtr) return NULL;

	for (itershtr = MasShtrih[pShtr->Nstr].begin (); itershtr != MasShtrih[pShtr->Nstr].end (); itershtr++)
	{

		if (*pShtr == *itershtr)
		{
			itershtr --;
			return &(*itershtr);
		}


	}


	return NULL;
}

// ��������� ���������� ������ � ������
TShtrih* TShtrihPicture::GetNextShtr (TShtrih *pShtr)
{
	list <TShtrih>::iterator itershtr;
	int count = 1;
	size_t size = MasShtrih[pShtr->Nstr].size ();

	for (itershtr = MasShtrih[pShtr->Nstr].begin (); itershtr != MasShtrih[pShtr->Nstr].end (); itershtr++)
	{
		if (count == size) return NULL;
        if (*pShtr == *itershtr)
		{
			itershtr ++;
			return &(*itershtr);
		}
		count ++;
	}

	return NULL;
}

int TShtrihPicture::GetDirectionInExtContour (TShtrih *pShtr, int prev, int &next)
{
	// prev, next:  1 - ����, 0 - �����.
    // ���� �������� ���� ����
	if (prev == 1)
	{
		// ���� ����� ���� �����
		if (pShtr->SWs > 0)
		{
			    // ���� � ������ ��� ����������� ������ 
				if (GetPrevShtr (pShtr) == NULL)
				{
					// ��������� ������� ����
					next = 1;
					// ������� ����
					return 0;
				}
				else 
				{   // ���� ������ ����������� ��������� �������:
					//   0-----0<===0-----0  
					// 0--------------------0
					if ((pShtr->begin <= GetBottomLeftShtrih (pShtr)->end) && (GetPrevShtr (pShtr)->end >= GetBottomLeftShtrih (pShtr)->begin))
					{
						// ��������� ������� �����
						next = 0;
						// ������� �����
						return 5;
					}
				   else
				   {  
					   // ��������� ������� ����
					   next = 1;
                       // ������� ����
					   return 0;
				    }
				}
			
		}
		else    // ����� SWp > 0
		{
			// ��������� ������� �����
			next = 0;
			// ������� �������
			return 1;		
		}
	}
	else
	{
		// ���� ������ ���� ������
		if (pShtr->SWp > 0)
		{
			    // ���� ��� ���������� ������
				if (GetNextShtr (pShtr) == NULL)
				{
					// ��������� ������� �����
					next = 0;
					// ������� �����
					return 2;
				}
				else
				{
					// ���� ������ ����������� ��������� �������:
					// 0--------------------0
					// 0-------0===>0-------0  
					if ((pShtr->end >= GetTopRightShtrih (pShtr)->begin) && (GetNextShtr (pShtr)->begin <= GetTopRightShtrih (pShtr)->end))
					{
						// ��������� ������� ����
						next = 1;
						// ������� ������
						return 4;
					}
					else
					{
						// ��������� ������� �����
						next = 0;
						// ������� �����
						return 2;
				     }
				}
		}
		else    // ����� SWp > 0
		{
			    // ������� ����
				next = 1;
				// ������� �����
				return 3;
		}
	}
}

int TShtrihPicture::GetDirectionInterContour (TShtrih *pShtr, int prev, int &next)
{
	// prev, next:  1 - ����, 0 - �����.
    
	// ���� ��� ������� ����
	if (prev == 1)
	{
		// ���� ����� ���� ������
		if (pShtr->SWs > 0)
		{
			// ���� ��� ����������� ������
			if (GetPrevShtr (pShtr) == NULL) 
				{
					// ��������� ������� ����
					next = 1;
					// ������� ����
					return 0;
				}
				else 
				{
					// ���� ������ ����������� ��������� �������:
					//   0-----0<===0-----0  
					// 0--------------------0
					if ((pShtr->begin <= GetBottomLeftShtrih (pShtr)->end) && (GetPrevShtr (pShtr)->end >= GetBottomLeftShtrih (pShtr)->begin))
					{
						// ��������� ������� �����
						next = 0;
						// ������� �����
						return 5;
					}
				   else
				   {   // ��������� ������� ���� 
					   next = 1;
					   // ������� ����
					   return 0;
					}
				}
		}
		else    // ����� SWp > 0
		{
			    // ��������� ������� �����
				next = 0;
				// ������� ������
				return 1;
		}
	}
	else
	{
		// ���� ������ ���� ������
		if (pShtr->SWp > 0)
		{
			    // ���� ��� ���������� ������
				if (GetNextShtr (pShtr) == NULL) 
				{
					// ��������� ������� �����
					next = 0;
					// ������� �����
					return 2;
				}
				else
				{
					// ���� ������ ����������� ��������� �������:
					// 0--------------------0
					// 0-------0===>0-------0  
					if ((pShtr->end >= GetTopRightShtrih (pShtr)->begin) && (GetNextShtr (pShtr)->begin <= GetTopRightShtrih (pShtr)->end))
					{
						// ��������� ������� ���� 
						next = 1;
						// ������� ������
						return 4;
					}
					else
					{
						// ��������� ������� ����� 
						next = 0;
						// ������� �����
						return 2;
					}
				}
		}
		else    // ����� SWs > 0
		{
			    // ��������� ������� ����
				next = 1;
				// ������� �����
				return 3;
		}
	}
}