#include "Polygon.h"
#include "Vector.h"
#include <Windows.h>

void Polygon::AddPoint(Point pPoint, bool bLastPoint /*= false*/)
{
	if (nullptr == pHeadPoint)
	{
		pHeadPoint = &pPoint;
	}
	else
	{
		Point* pCurrPoint = pHeadPoint;
		while (nullptr != pCurrPoint->pNext)// Gets Point having pNext = nullptr
		{
			pCurrPoint = pCurrPoint->pNext;
		}
		pCurrPoint->pNext = &pPoint;
		pPoint.pPrev = pCurrPoint;//
		dEdgesCount++;
		if (true == bLastPoint)
		{
			pPoint.pNext = pHeadPoint;
			pHeadPoint->pPrev = &pPoint;//
			dEdgesCount++;
		}
	}
}

Point Polygon::GetEdgeStart(int iEdgeCount)
{
	Point pHeadPoint = GetHeadPoint();
	if (0 == iEdgeCount)
		return pHeadPoint;
	else
	{
		Point* pStartPoint = &pHeadPoint;
		for (int iEdgeIndex = 0; iEdgeIndex < iEdgeCount; iEdgeIndex++)
		{
			pStartPoint = pStartPoint->pNext;
		}
		return *pStartPoint;
	}
}

void Polygon::GetEdgeEndsFromPolygon(int iEdgeCount, double& dX1, double& dY1, double& dX2, double& dY2)
{
	Point pStartPoint = GetEdgeStart(iEdgeCount);
	dX1 = pStartPoint.GetX();
	dY1 = pStartPoint.GetY();

	Point* pEndPoint = pStartPoint.pNext;
	dX2 = pEndPoint->GetX();
	dY2 = pEndPoint->GetY();
}

bool Polygon::IsPointInsideUsingRayCasting(Point& pIsPointInside, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge)
{
	if (pIsPointInside.GetZ() < 0.0)
		return false;

	double dXp = pIsPointInside.GetX();
	double dYp = pIsPointInside.GetY();

	int dTotalEdges = GetEdgesCount();

	int iIntersectionCount = 0;
	bool bFirstIntersectionAtVertex = false;
	for (int iEdgeCount = 0; iEdgeCount < dTotalEdges; iEdgeCount++)
	{
		double dX1 = 0.0;
		double dY1 = 0.0;
		double dX2 = 0.0;
		double dY2 = 0.0;
		GetEdgeEndsFromPolygon(iEdgeCount, dX1, dY1, dX2, dY2);

		if (dY1 == dY2)
		{
			// The edge is hori
			if (pIsPointInside.GetY() == dY1)
			{
				// Point is collinear w.r.t. the edge
				if (pIsPointInside.GetX() == dX1 || pIsPointInside.GetX() == dX2)
				{
					// Point satisfies one of the end of the edge
					bOnVertex = true;
				}
				else if ((pIsPointInside.GetX() < dX1) != (pIsPointInside.GetX() < dX2))
				{
					// Point is on the edge
					bOnEdge = true;
				}
				else
				{
					bOnLineOfEdge = true;
				}
			}
		}
		else
		{
			double pIntersectionPointX = ((dX2 - dX1) * (pIsPointInside.GetY() - dY1) / (dY2 - dY1)) + dX1;
			Point pIntersectionPoint(pIntersectionPointX, pIsPointInside.GetY(), 0.0);
			if ((pIsPointInside.GetX() == dX1 && pIsPointInside.GetY() == dY1) || (pIsPointInside.GetX() == dX2 && pIsPointInside.GetY() == dY2))
			{
				// Point is on the vertex
				bOnVertex = true;
				return false;
			}
			else if (((pIsPointInside.GetY() < dY1) != (pIsPointInside.GetY() < dY2)) && (pIsPointInside.GetX() == pIntersectionPointX))
			{
				// Point is on the edge
				bOnEdge = true;
				return false;
			}
			else if ((pIsPointInside.GetY() == dY1) && (pIsPointInside.GetX() < pIntersectionPointX))
			{
				// Point lies in the line of Polygon vertex and inside the Polygon
				if (false == bFirstIntersectionAtVertex)
				{
					iIntersectionCount++;
					bFirstIntersectionAtVertex = true;
				}
			}
			else if (((pIsPointInside.GetY() < dY1) != (pIsPointInside.GetY() < dY2)) && (pIsPointInside.GetX() < pIntersectionPointX) 
				&& (pIsPointInside.GetY() != dY1) && (pIsPointInside.GetY() != dY2))
				// Point lies inside the Polygon
				iIntersectionCount++;
		}
	}

	if (iIntersectionCount % 2 == 1)
		return true;
	else
		return false;
}

bool Polygon::IsPointInsideUsingVectorAlgebra(Point& pIsPointInside, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge)
{
	// This may fail in some conditions like
	/*________________    _________
	  |	             |    |       |
	  |      . point |____|       |
	  |___________________________|
	*/
	if (pIsPointInside.GetZ() != 0.0 || (false == IsConvexType()))
	{
		MessageBox(NULL, L"Vector Algebra method requires convex polygon.", L"Invalidity", MB_OK | MB_ICONERROR);
		return false;
	}

	double dXp = pIsPointInside.GetX();
	double dYp = pIsPointInside.GetY();

	int dTotalEdges = GetEdgesCount();

	int iIntersectionCount = 0;
	double dPreviousCrossProd = NULL;
	for (int iEdgeCount = 0; iEdgeCount < dTotalEdges; iEdgeCount++)
	{
		double dX1 = 0.0;
		double dY1 = 0.0;
		double dX2 = 0.0;
		double dY2 = 0.0;
		GetEdgeEndsFromPolygon(iEdgeCount, dX1, dY1, dX2, dY2);

		// Vectors formation
		double vP1X = (dXp - dX1);
		double vP1Y = (dYp - dY1);
		double v21X = (dX2 - dX1);
		double v21Y = (dY2 - dY1);

		double dCrossProduct = (vP1X * v21Y) - (vP1Y * v21X);
		if ((dCrossProduct > 0 && dPreviousCrossProd > 0) || (dCrossProduct < 0 && dPreviousCrossProd < 0))
			continue;
		else if (0 == iEdgeCount)
		{
			dPreviousCrossProd = dCrossProduct;
			continue;
		}
		else
			return false;
		dPreviousCrossProd = dCrossProduct;
	}

	return true;
}

bool Polygon::IsPointInside(Point& pIsPointInside, MethodType eMethodType, bool& bOnEdge, bool& bOnVertex, bool& bOnLineOfEdge)
{
	bool bIsPointInside = false;
	if (RAY_CASTING == eMethodType)
		bIsPointInside = IsPointInsideUsingRayCasting(pIsPointInside, bOnEdge, bOnVertex, bOnLineOfEdge);
	else if (VECTOR_ALGEBRA == eMethodType)
		bIsPointInside = IsPointInsideUsingVectorAlgebra(pIsPointInside, bOnEdge, bOnVertex, bOnLineOfEdge);

	return bIsPointInside;
}

bool Polygon::IsConvexType()
{
	Point* pHeadPrevPoint = pHeadPoint->pPrev;
	Point* pHeadNextPoint = pHeadPoint->pNext;

	double vHead1PX = (pHeadPrevPoint->GetX() - pHeadPoint->GetX());
	double vHead1PY = (pHeadPrevPoint->GetY() - pHeadPoint->GetY());
	double vHead2PX = (pHeadNextPoint->GetX() - pHeadPoint->GetX());
	double vHead2PY = (pHeadNextPoint->GetY() - pHeadPoint->GetY());

	double dHeadCrossProd = (vHead1PX * vHead2PY) - (vHead1PY * vHead2PX);

	Point* pCurrPoint = pHeadPoint->pNext;
	while (pCurrPoint != pHeadPoint)
	{
		Point* pPrevPoint = pCurrPoint->pPrev;
		Point* pNextPoint = pCurrPoint->pNext;

		double v1PX = (pPrevPoint->GetX() - pCurrPoint->GetX());
		double v1PY = (pPrevPoint->GetY() - pCurrPoint->GetY());
		double v2PX = (pNextPoint->GetX() - pCurrPoint->GetX());
		double v2PY = (pNextPoint->GetY() - pCurrPoint->GetY());

		double dCrossProd = (v1PX * v2PY) - (v1PY * v2PX);

		if ((dCrossProd > 0 && dHeadCrossProd < 0) || (dCrossProd < 0 && dHeadCrossProd > 0))
			return false;

		pCurrPoint = pCurrPoint->pNext;
	}

	return true;
}