#ifndef TShtrih_HH
#define TShtrih_HH
#pragma once

// класс штрих
class TShtrih {

public:
	 int Nstr;       // номер строки
     int begin;      // начало штриха
     int end;        // конец штриха
     int SWp;        // коэффициент связи со штрихом предыдущей строки
     int SWs;        // коэффициент связи со штрихом следующей строки

public:
	~TShtrih ()
	{

	}
	TShtrih (const TShtrih &obj)
	{
		Nstr = obj.Nstr;
		begin = obj.begin;
		end = obj.end;
		SWp = obj.SWp;
		SWs = obj.SWs;
	}
	TShtrih(int _Nstr, int _begin, int _end, int _SWp, int _SWs):Nstr (_Nstr), begin (_begin), end (_end), SWp (_SWp), SWs (_SWs) {}
	TShtrih ()
	{
		Nstr = 0;
		begin = 0;
		end = 0;
		SWp = 0;
		SWs = 0;
	}
	friend bool operator == (TShtrih shtr1, TShtrih shtr2)
	{
		 return ((shtr1.Nstr == shtr2.Nstr) && (shtr1.begin == shtr2.begin) && (shtr1.end == shtr2.end) && (shtr1.SWp == shtr2.SWp) && (shtr1.SWs == shtr2.SWs));
	}
	friend bool operator != (TShtrih shtr1, TShtrih shtr2)
	{
		return ((shtr1.Nstr != shtr2.Nstr) || (shtr1.begin != shtr2.begin) || (shtr1.end != shtr2.end) || (shtr1.SWp != shtr2.SWp) || (shtr1.SWs != shtr2.SWs));
	}

};
#endif
