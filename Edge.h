#pragma once
#include "Point.h"

class Edge
{
private: 
	Point pPoint1;
	Point pPoint2;

public:
	Edge() : pPoint1(0.0, 0.0, 0.0), pPoint2(0.0, 0.0, 0.0) {};

	Edge(Point _pPoint1, Point _pPoint2) : pPoint1(_pPoint1), pPoint2(_pPoint2) {};

	void SetX1(double dSetX1) { pPoint1.SetX(dSetX1); }
	double GetX1() { return pPoint1.GetX(); }

	void SetY1(double dSetY1) { pPoint1.SetX(dSetY1); }
	double GetY1() { return pPoint1.GetY(); }

	void SetZ1(double dSetZ1) { pPoint1.SetX(dSetZ1); }
	double GetZ1() { return pPoint1.GetZ(); }

	void SetX2(double dSetX2) { pPoint2.SetX(dSetX2); }
	double GetX2() { return pPoint2.GetX(); }

	void SetY2(double dSetY2) { pPoint2.SetX(dSetY2); }
	double GetY2() { return pPoint2.GetY(); }

	void SetZ2(double dSetZ2) { pPoint2.SetX(dSetZ2); }
	double GetZ2() { return pPoint2.GetZ(); }
};

