#include <iostream>
#include <vector>
#include "TriangularFace.h"

void TriangularFace::SetPoint(int iIndex, Point pPoint)
{
	if (iIndex > 3)
		return;

	if (0 == iIndex)
		pP0 = pPoint;
	else if (1 == iIndex)
		pP1 = pPoint;
	else if (2 == iIndex)
		pP2 = pPoint;
}

Point TriangularFace::GetPoint(int iIndex)
{
	if (iIndex < 3)
	{
		if (0 == iIndex)
			return this->pP0;
		else if (1 == iIndex)
			return this->pP1;
		else if (2 == iIndex)
			return this->pP2;
	}

	return Point(0,0,0);
}

void TriangularFace::SetPointsSet(std::vector<Point>& vOtherSetOfPoints)
{
	vSetOfPoints = vOtherSetOfPoints;
}

bool TriangularFace::HasSetOfPoints()
{
	return bHasSetOfPoints;
}

void TriangularFace::SetHasSetOfPoints(bool bSetHasSetOfPoints)
{
	bHasSetOfPoints = bSetHasSetOfPoints;
}
