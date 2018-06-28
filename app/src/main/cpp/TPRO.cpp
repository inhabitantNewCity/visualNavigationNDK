#include "TPRO.h"

TPRO::TPRO()
{
    NumShtrih = 0;
}

TPRO::TPRO(const TPRO &obj) 
{
	NumShtrih = obj.NumShtrih;
	list <TShtrih>::iterator itershtr;

	ListShtrih_PRO = new list <TShtrih> ();
	ListShtrih_PRO = obj.ListShtrih_PRO;

}


TPRO::~TPRO()
{
	delete ListShtrih_PRO;
}
// ���������� ������ ������� � ����� ���
void TPRO::AddListShtrih(list<TShtrih> *addlistshtrih)      
{
	list <TShtrih>::iterator itershtr;

	for (itershtr = addlistshtrih->begin (); itershtr != addlistshtrih->end (); itershtr++)
        AddShtrih(*itershtr);
}
// ���������� � ����� ������ �������
void TPRO::AddShtrih(TShtrih &addshtrih)          
{
	ListShtrih_PRO->push_back(addshtrih);
    NumShtrih++;
}
/*
void TPRO::DrawPROOnPicture(Bitmap imagein, int zoom, Color color)  // ��������� ���
{
    for(int i = 0; i < NumShtrih; i++)
    {
        TShtrih.DrawLine(imagein, ListShtrih_PRO[i].Nstr * zoom, ListShtrih_PRO[i].begin * zoom, (ListShtrih_PRO[i].end + 1) * zoom, 1, color);
    }
}
*/
// �������� ����������� �� ����� ���
bool TPRO::CheckPRO_Shtrih(TShtrih &shtrih)       
{
	list <TShtrih>::iterator itershtr;

    bool res = false;
    for (itershtr = ListShtrih_PRO->begin (); itershtr != ListShtrih_PRO->end (); itershtr++)
    {
        if (*itershtr == shtrih)
        {
            res = true;
            break;
        }
		
    }
    return res;
}
// ������� ����������� ���
void TPRO::Clear()            
{
    ListShtrih_PRO->clear();
    NumShtrih = 0;
}
