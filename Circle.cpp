#include <iostream>
#include "Circle.h"

bool Circle::IsPointInside(Point& pIsPointInside, bool& bOnCircum)
{
	double dXp = pIsPointInside.GetX();
	double dYp = pIsPointInside.GetY();

	Point pCircCenter = GetCenter();
	double dXCircCenter = pCircCenter.GetX();
	double dYCircCenter = pCircCenter.GetY();

	double dCircRad = GetRadius();

	double dDistFromCircCenter = sqrt((dXp - dXCircCenter) * (dXp - dXCircCenter) + (dYp - dYCircCenter) * (dYp - dYCircCenter));

	if (dDistFromCircCenter == dCircRad)
	{
		bOnCircum = true;
		return false;
	}
	else if (dDistFromCircCenter > dCircRad)
		return false;
	else if (dDistFromCircCenter < dCircRad)
		return true;
}
