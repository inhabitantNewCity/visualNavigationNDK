#include "TBlob.h"

TBlob::TBlob()
{
    NumShtrih_Blob = 0;
}

TBlob::TBlob (const TBlob &Obj)
{
	NumShtrih_Blob = Obj.NumShtrih_Blob;
	ListShtrih_Blob = new list<TShtrih>();
	ListShtrih_Blob = Obj.ListShtrih_Blob;
    
}

TBlob::~TBlob()
{
    delete ListShtrih_Blob;
}
// проверяет нужно ли добавлять штрих 
bool TBlob::CheckingToAddShtrih(TShtrih &addshtrih)  
{
	list <TShtrih>::iterator itershtr;

    bool res = false;
	for (itershtr = ListShtrih_Blob->begin (); itershtr != ListShtrih_Blob->end (); itershtr++)
    {
        if (((((*itershtr).begin <= addshtrih.begin) && ((*itershtr).end >= addshtrih.begin)) ||
            (((*itershtr).end >= addshtrih.end) && ((*itershtr).begin <= addshtrih.end)) ||
            (((*itershtr).begin >= addshtrih.begin) && ((*itershtr).end <= addshtrih.end)) ||
            (((*itershtr).begin <= addshtrih.begin) && ((*itershtr).end >= addshtrih.end))) &&
            ((((*itershtr).Nstr - addshtrih.Nstr) == 1) || ((addshtrih.Nstr - (*itershtr).Nstr) == 1)))
        {
            res = true;
            break;
        }
    }
    return res;
}
// добавление списка штрихов в конец блоба
void TBlob::AddListShtrih(list<TShtrih> *addlistshtrih)      
{
	list <TShtrih>::iterator itershtr;

	for (itershtr = addlistshtrih->begin (); itershtr != addlistshtrih->end (); itershtr++)
        AddShtrih(*itershtr);
}
// добавление в конец списка штрихов
void TBlob::AddShtrih(TShtrih &addshtrih)          
{
	ListShtrih_Blob->push_back(addshtrih);
    NumShtrih_Blob++;
}

// создание списка блобов
list<TBlob>* TBlob::CreateListBLOB(TShtrihPicture &ShPicture)    
{
    list<TBlob> *ResultListBlob = new list<TBlob>();
    list<TPRO> *ListPRO = new list<TPRO>();
    list<TShtrih> *ListConnectShtrih = new list<TShtrih>();

    CreateListPRO_and_ConnectShtr(ListConnectShtrih, ListPRO, ShPicture);
    
	list <TPRO>::iterator iterpro;
    // копируем список ПРО в список блобов
	for (iterpro = ListPRO->begin (); iterpro != ListPRO->end (); iterpro++)
    {
        TBlob *ADD_Blob = new TBlob();
        ADD_Blob->ListShtrih_Blob = new list<TShtrih>();
        ADD_Blob->AddListShtrih((*iterpro).ListShtrih_PRO);
		ResultListBlob->push_back((*ADD_Blob));
    }
    ListPRO->clear();
    int Pos = -1;    // позиция в списке ПРО первого ПРО, к к которму был добавлен штрих
    // добавляем соединительные штрихи в ПРО и одновременно соединяем ПРО, связанные штрихом
	list<TShtrih>::iterator iterlcshtr;
	list<TBlob>::iterator iterblob;
	list<TBlob>::iterator itertemp;
	list<TBlob>::iterator iterprev;

	for (iterlcshtr = ListConnectShtrih->begin (); iterlcshtr != ListConnectShtrih->end (); iterlcshtr++)
    {
		for (iterblob = ResultListBlob->begin (); iterblob != ResultListBlob->end (); iterblob++)
        {
            if ((*iterblob).CheckingToAddShtrih(*iterlcshtr))
            {
                if (Pos == -1)
                {
                    //Pos = j;
					itertemp = iterblob;
                    (*iterblob).AddShtrih(*iterlcshtr);
					Pos = 1;
                }
                else
                {
					iterprev = iterblob;
                    iterprev --;
					
                    // соединяем ПРО под номером Pos и ПРО под номером j; j-ый ПРО удаляем из списка 
                    itertemp->AddListShtrih(iterblob->ListShtrih_Blob);
					ResultListBlob->erase (iterblob);
					iterblob = iterprev;
                }
            }
			
        }
        Pos = -1;
    }
    delete ListConnectShtrih;

    return ResultListBlob;

}

/*
// отрисовка блоба
void TBlob::DrawBlobOnPicture(Bitmap ^imagein, list <TShtrih> * ListShtrih, int zoom, Color color)  
{
	list <TShtrih>::iterator itershtr;

	for (itershtr = ListShtrih->begin (); itershtr != ListShtrih->end (); itershtr++)
    {
        DrawLine(imagein, itershtr->Nstr * zoom, itershtr->begin * zoom, (itershtr->end + 1) * zoom, 1, color);
    }
	
}
*/

/*
void TBlob::DrawLine(Bitmap ^imagein, int nstr, int start, int end, int width, Color color)   // отрисовка линии на изображении imagein в строке nstr
{
    for (int i = 0; i < width; i++)
    {
        for (int k = start; k < end; k++)
        {
            imagein->SetPixel(k, nstr + i, color);
        }
    }
}
*/
// поиск номера ПРО в списке, который содержит shtrih; возвращает -1, если не принадлежит ни одному ПРО из списка
int TBlob::SearchInList_PROWithShtrih(list<TPRO> *listpro, TShtrih shtrih) 
{
    int res = -1;
	list <TPRO>::iterator iterpro;
	int i = 0;

	for (iterpro = listpro->begin (); iterpro != listpro->end (); iterpro ++)
    {
		if (iterpro->CheckPRO_Shtrih(shtrih))
        {
            res = i;
            break;
        }
		i++;
    }
    return res;
}

// создание списка ПРО и списка соеденительных штрихов для изображения
void TBlob::CreateListPRO_and_ConnectShtr(list<TShtrih> *ListConnectShtrihs, list<TPRO> *ListPro, TShtrihPicture ShPicture) 
{
    //ListPro = new list<TPRO>();
    int PosPRO = 0;                         // номер ПРО в списке ПРО
	list <TShtrih>::iterator itershtr;

	int j = 0;


    for (int i = 0; i < ShPicture.Height; i++)
    {
		j = 0;
		
		for (itershtr = ShPicture.MasShtrih[i].begin (); itershtr != ShPicture.MasShtrih[i].end (); itershtr++)
		{
			

            if ((itershtr->SWp == 0) && (itershtr->SWs == 1))
            {
                // создаем новый ПРО, который будет входить в блоб
                TPRO *ADD_PRO = new TPRO();
                ADD_PRO->ListShtrih_PRO = new list<TShtrih>();
                TShtrih *ADDshtrihPRO = new TShtrih(
                                itershtr->Nstr,
                                itershtr->begin,
                                itershtr->end,
                                itershtr->SWp,
                                itershtr->SWs);
				ADD_PRO->ListShtrih_PRO->push_back(*ADDshtrihPRO);
                ADD_PRO->NumShtrih++;
                // добавляем в блоб новый ПРО
				ListPro->push_back(*ADD_PRO);
            }
            if (((itershtr->SWp == 1) && (itershtr->SWs == 1)) ||
                ((itershtr->SWp == 1) && (itershtr->SWs == 0)))
            {
				list <TShtrih>::iterator itershtrp;
                // поиск штриха в предыдущей строке, связанный с текущим штрихом
				for (itershtrp = ShPicture.MasShtrih[i - 1].begin(); itershtrp != ShPicture.MasShtrih[i - 1].end (); itershtrp++)
                {

                    if (((itershtrp->begin <= itershtr->begin) && (itershtrp->end >= itershtr->begin)) ||
                        ((itershtrp->end >= itershtr->end) && (itershtrp->begin <= itershtr->end)) ||
                        ((itershtrp->begin >= itershtr->begin) && (itershtrp->end <= itershtr->end)) ||
                        ((itershtrp->begin <= itershtr->begin) && (itershtrp->end >= itershtr->end)))
                    {
                        PosPRO = SearchInList_PROWithShtrih(ListPro, *itershtrp);
                        if (PosPRO == -1)
                        {
                            // верхний штрих не подходит
                            // создаем новый ПРО
                            TPRO *ADD_PRO = new TPRO();
                            ADD_PRO->ListShtrih_PRO = new list<TShtrih>();
                            TShtrih *ADDshtrihPRO = new TShtrih(
                                            itershtr->Nstr,
											itershtr->begin,
											itershtr->end,
											itershtr->SWp,
											itershtr->SWs);
                            ADD_PRO->ListShtrih_PRO->push_back(*ADDshtrihPRO);
                            ADD_PRO->NumShtrih++;
                            // добавляем в список ПРО новый ПРО
                            ListPro->push_back(*ADD_PRO);
                        }
                        else
                        {
						    list <TPRO>::iterator iterpro = ListPro->begin ();

                            // добавляем еще один штрих в ПРО с номером PosPRO
                            TShtrih *ADDshtrihPRO = new TShtrih(
                                            itershtr->Nstr,
											itershtr->begin,
											itershtr->end,
											itershtr->SWp,
											itershtr->SWs);
							
							for (int vi = 0; vi < PosPRO; vi++)
								iterpro++;

							iterpro->ListShtrih_PRO->push_back(*ADDshtrihPRO);
							iterpro->NumShtrih++;
                        }
                        break;
                    }
                }
            }
            if ((itershtr->SWp > 1) || (itershtr->SWs > 1))
            {
                TShtrih *ADDshtrihPRO = new TShtrih(
					                        itershtr->Nstr,
											itershtr->begin,
											itershtr->end,
											itershtr->SWp,
											itershtr->SWs);
				ListConnectShtrihs->push_back(*ADDshtrihPRO);
            }
			j++;
        }
    }
	int f = 0;
}
// очистка содержимого блоба
void TBlob::Clear()             
{
	ListShtrih_Blob->clear();
    NumShtrih_Blob = 0;
}
// проверка на принадлежность штриха блобу
bool TBlob::Check_BlobHasShtrih(TShtrih &shtrih)     
{
    bool res = false;
	list <TShtrih>::iterator itershtr;

	
	for (itershtr = ListShtrih_Blob->begin (); itershtr != ListShtrih_Blob->end (); itershtr++)
    {
        if ((*itershtr) == shtrih)
        {
            res = true;
            break;
        }
    }

    return res;
}