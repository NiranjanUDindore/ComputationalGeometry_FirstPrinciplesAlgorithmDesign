#pragma once
#include "Point.h"
#include "Shape.h"
class Circle : public Shape
{
private:
	Point pCenter;
	double dRad = 0.0;

public:
	Circle() : pCenter(0.0, 0.0, 0.0) {};

	Circle(Point _pCenter, double _dRad) : pCenter(_pCenter), dRad(_dRad) {};

	void SetRadius(double dSetRad) { dRad = dSetRad; }
	double GetRadius() { return dRad; }

	void SetCenter(Point pSetCenter) { pCenter = pSetCenter; }
	Point GetCenter() { return pCenter; }

	bool IsPointInside(Point& pIsPointInside, bool& bOnEdge);
};

