#pragma once
#include "Point.h"
#include "Shape.h"
#include <vector>
class Polygon : public Shape
{
private:
	Point* pHeadPoint;
	int dEdgesCount = 0;

public:
	Polygon() : pHeadPoint(nullptr) {};

	Polygon(Point& _pHeadPoint) : pHeadPoint(&_pHeadPoint) {};

	void SetHeadPoint(Point& pSetHeadPoint) { pHeadPoint = &pSetHeadPoint; }
	Point GetHeadPoint() { return *pHeadPoint; }

	void AddPoint(Point pNextPoint, bool bLastPoint = false);
	int GetEdgesCount() { return dEdgesCount; }

	enum MethodType
	{
		RAY_CASTING = 0,
		VECTOR_ALGEBRA = 1
	};

	Point GetEdgeStart(int iEdgeCount);
	void GetEdgeEndsFromPolygon(int iEdgeCount, double& dX1, double& dY1, double& dX2, double& dY2);
	bool IsPointInside(Point& pIsPointInside, MethodType eMethodType, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge);
	bool IsPointInsideUsingRayCasting(Point& pIsPointInside, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge);
	bool IsPointInsideUsingVectorAlgebra(Point& pIsPointInside, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge);
	bool IsConvexType();
};

