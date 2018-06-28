#ifndef TPoint_HH
#define TPoint_HH
#pragma once
#include <stdio.h>
// класс точка
class T_Point {

private:
	short x;
	short y;

public:
	T_Point (short _x = 0, short _y = 0):x(_x), y(_y){}

	T_Point (const T_Point &obj)
	{
		x = obj.x;
		y = obj.y;
	}

	friend bool operator == (T_Point point1, T_Point point2)
	{
		return ((point1.x == point2.x) && (point1.y == point2.y));
	}
	friend bool operator != (T_Point point1, T_Point point2)
	{
		return ((point1.x != point2.x) || (point1.y != point2.y));
	}

	friend bool operator < (T_Point point1, T_Point point2)
	{
		return ((point1.x < point2.x) && (point1.y < point2.y));
	}

	friend bool operator <= (T_Point point1, T_Point point2)
	{
		return ((point1.x <= point2.x) && (point1.y <= point2.y));
	}

	T_Point operator - (T_Point point)
	{
		T_Point temp;
		short vx = 0;
		short vy = 0;

		vx = x - point.GetX();
		vy = y - point.GetY();
		temp.SetX(vx);
		temp.SetY(vy);

		return temp;
	}

    T_Point operator = (T_Point point)
	{
		x = point.x;
		y = point.y;

		return *this;
	}

	int GetX ()
	{
		return x;
	}
	int GetY ()
	{
		return y;
	}
	void SetX (short _x)
	{
		x = _x;
	}
	void SetY (short _y)
	{
		y = _y;
	}

};
#endif
