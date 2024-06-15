#pragma once
#include "Point.h"
class TriangularFace
{
private:
	Point pP0;
	Point pP1;
	Point pP2;

	std::vector<Point> vSetOfPoints;
	bool bHasSetOfPoints = false;

public:
	TriangularFace() : pP0(0.0, 0.0, 0.0), pP1(0.0, 0.0, 0.0), pP2(0.0, 0.0, 0.0) {};

	TriangularFace(Point _pP0, Point _pP1, Point _pP2) : pP0(_pP0), pP1(_pP1), pP2(_pP2) {};

	void SetPoint(int iIndex, Point pPoint);
	Point GetPoint(int iIndex);

	// For 3D Convex Hull
	std::vector<Point> GetSetOfPoints() { return vSetOfPoints; }
	void AddPointInSet(Point& pPoint) { vSetOfPoints.push_back(pPoint); }
	void SetPointsSet(std::vector<Point>& vOtherSetOfPoints);
	void SetHasSetOfPoints(bool bSetHasSetOfPoints);
	bool HasSetOfPoints();
};

