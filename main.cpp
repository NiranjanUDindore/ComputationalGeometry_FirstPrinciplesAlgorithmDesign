#include <iostream>
#include "Vector.h"
#include "Polygon.h"
#include "Circle.h"
#include "TriangularFace.h"
#include "Tetrahedron.h"
#include "Edge.h"
#include "Matrix.h"
#include "Line.h"
#include "Plane.h"
#include "Tria.h"
#include "Node.h"
#include "Edge2.h"
#include <cmath>
#include <algorithm>
#include <GLFW/glfw3.h>

const double PI = 3.14159265358979323846;
const double tolerance = 1E-06;

void drawLine()
{
	/*glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES);
	glVertex2f(-100, -100);
	glVertex2f(200, 200);

	glVertex2f(200, 200);
	glVertex2f(200, 0);

	glVertex2f(200, 0);
	glVertex2f(-200, 0);
	glEnd();
	glFlush();*/
}

void drawPoint()
{
	////glClear(GL_COLOR_BUFFER_BIT);
	//glColor3f(1.0, 0.0, 0.0);
	//glPointSize(10.0);
	//glBegin(GL_POINTS);
	//glVertex2f(100, 100);	
	//glEnd();
	//glFlush();
}

void GetMinYWithIndex(std::vector<Point>& vPoints, double& dMinY, int& iMinYPointIndex)
{
	for (int iPointIndex = 0; iPointIndex < vPoints.size(); iPointIndex++)
	{
		if (dMinY > vPoints.at(iPointIndex).GetY())
		{
			dMinY = vPoints.at(iPointIndex).GetY();
			iMinYPointIndex = iPointIndex;
		}
	}
}

void GetMaxYWithIndex(std::vector<Point>& vPoints, double& dMaxY, int& iMaxYPointIndex)
{
	for (int iPointIndex = 0; iPointIndex < vPoints.size(); iPointIndex++)
	{
		if (dMaxY < vPoints.at(iPointIndex).GetY())
		{
			dMaxY = vPoints.at(iPointIndex).GetY();
			iMaxYPointIndex = iPointIndex;
		}
	}
}

void CreateHalfConvexHullUsingGiftWrap(std::vector<Point>& vPoints, int iXAxisDirection, int iMinYPointIndex, double dMaxY, std::vector<Point>& vHalfConvexHullPoints)
{
	Point pHeadPointRight = vPoints.at(iMinYPointIndex);
	vHalfConvexHullPoints.push_back(pHeadPointRight);

	for (int iConvexHullPointIndex = 0; iConvexHullPointIndex < vHalfConvexHullPoints.size(); iConvexHullPointIndex++)
	{
		double dMinAngle = 180.0;
		int iMinAnglePointIndex = -1;
		for (int iPointIndex = 0; iPointIndex < vPoints.size(); iPointIndex++)
		{
			bool bPointAlreadyInConvexHull = false;
			for (int iConvexHullPointIndex = 0; iConvexHullPointIndex < vHalfConvexHullPoints.size(); iConvexHullPointIndex++)
			{
				if (vPoints.at(iPointIndex) == vHalfConvexHullPoints.at(iConvexHullPointIndex))
				{
					bPointAlreadyInConvexHull = true;
					break;
				}
			}

			if (false == bPointAlreadyInConvexHull)
			{
				double dAngle = 0.0;
				if (vHalfConvexHullPoints.size() < 2)
				{
					// Get cross product with X axis
					double dP2HeadPointX = vPoints.at(iPointIndex).GetX() - vHalfConvexHullPoints.at(iConvexHullPointIndex).GetX();
					double dP2HeadPointY = vPoints.at(iPointIndex).GetY() - vHalfConvexHullPoints.at(iConvexHullPointIndex).GetY();

					double dDotProd = (iXAxisDirection * 1 * dP2HeadPointX) + (0 * dP2HeadPointY) + (0 * 0);
					double dMagXAxis = sqrt(1 * 1 + 0 * 0 + 0 * 0);
					double dMagP2HeadPoint = sqrt(dP2HeadPointX * dP2HeadPointX + dP2HeadPointY * dP2HeadPointY + 0.0);
					dAngle = acos(dDotProd / (dMagXAxis * dMagP2HeadPoint)) * 180 / PI;
				}
				else
				{
					double dP2P1X = vPoints.at(iPointIndex).GetX() - vHalfConvexHullPoints.at(iConvexHullPointIndex).GetX();
					double dP2P1Y = vPoints.at(iPointIndex).GetY() - vHalfConvexHullPoints.at(iConvexHullPointIndex).GetY();

					double dP1P0X = vHalfConvexHullPoints.at(iConvexHullPointIndex).GetX() - vHalfConvexHullPoints.at(iConvexHullPointIndex - 1).GetX();
					double dP1P0Y = vHalfConvexHullPoints.at(iConvexHullPointIndex).GetY() - vHalfConvexHullPoints.at(iConvexHullPointIndex - 1).GetY();

					double dDotProd = (dP2P1X * dP1P0X) + (dP2P1Y * dP1P0Y) + (0 * 0);
					double dMagP1P0 = sqrt(dP1P0X * dP1P0X + dP1P0Y * dP1P0Y + 0.0);
					double dMagP2P1 = sqrt(dP2P1X * dP2P1X + dP2P1Y * dP2P1Y + 0.0);
					dAngle = acos((dDotProd / (dMagP1P0 * dMagP2P1))) * 180 / PI;
				}

				if (dMinAngle > dAngle)
				{
					dMinAngle = dAngle;
					iMinAnglePointIndex = iPointIndex;
				}
			}
		}

		if (-1 != iMinAnglePointIndex)
		{
			vHalfConvexHullPoints.push_back(vPoints.at(iMinAnglePointIndex));
			if (dMaxY == vPoints.at(iMinAnglePointIndex).GetY())
				break;
		}
	}
}

void CombineConvexHullInRightVector(std::vector<Point>& vConvexHullPointsRight, std::vector<Point>& vConvexHullPointsLeft, std::vector<Point>& vConvexHullPoints)
{
	int dSizeConvexHullRight = vConvexHullPointsRight.size();
	int dSizeConvexHullLeft = vConvexHullPointsLeft.size();

	for (int iRightHullPointIndex = 0; iRightHullPointIndex < vConvexHullPointsRight.size(); iRightHullPointIndex++)
	{
		for (int iLeftHullPointIndex = 0; iLeftHullPointIndex < vConvexHullPointsLeft.size(); iLeftHullPointIndex++)
		{
			if (vConvexHullPointsLeft.at(iLeftHullPointIndex) == vConvexHullPointsRight.at(iRightHullPointIndex))
			{
				vConvexHullPointsLeft.erase(vConvexHullPointsLeft.begin() + iLeftHullPointIndex);
				break;
			}
		}
	}
	std::reverse(vConvexHullPointsLeft.begin(), vConvexHullPointsLeft.end());

	for (int iLeftConvexHullPointsIndex = 0; iLeftConvexHullPointsIndex < vConvexHullPointsLeft.size(); iLeftConvexHullPointsIndex++)
	{
		vConvexHullPointsRight.push_back(vConvexHullPointsLeft.at(iLeftConvexHullPointsIndex));
	}
	vConvexHullPoints = vConvexHullPointsRight;
}

void CreateConvexHullUsingGiftWrap(std::vector<Point>& vPoints, std::vector<Point>& vConvexHullPoints)
{
	double dMinY = 0;
	int iMinYPointIndex = 0;
	GetMinYWithIndex(vPoints, dMinY, iMinYPointIndex);

	double dMaxY = 0;
	int iMaxYPointIndex = 0;
	GetMaxYWithIndex(vPoints, dMaxY, iMaxYPointIndex);

	std::vector<Point> vConvexHullPointsRight;
	CreateHalfConvexHullUsingGiftWrap(vPoints, 1, iMinYPointIndex, dMaxY, vConvexHullPointsRight);

	std::vector<Point> vConvexHullPointsLeft;
	CreateHalfConvexHullUsingGiftWrap(vPoints, -1, iMinYPointIndex, dMaxY, vConvexHullPointsLeft);

	// Combine the left and right side of convex hull
	CombineConvexHullInRightVector(vConvexHullPointsRight, vConvexHullPointsLeft, vConvexHullPoints);
}

void SortVectorIncreasinglyWithY(std::vector<Point>& vPointsGrahamScanAlgo)
{
	for (int iLoopIndex = 0; iLoopIndex < vPointsGrahamScanAlgo.size() - 1; iLoopIndex++)
	{
		for (int iPointIndex = 0; iPointIndex < (vPointsGrahamScanAlgo.size() - 1 - iLoopIndex); iPointIndex++)
		{
			if ((vPointsGrahamScanAlgo.at(iPointIndex).GetY() > vPointsGrahamScanAlgo.at(iPointIndex + 1).GetY()) ||
				(vPointsGrahamScanAlgo.at(iPointIndex).GetY() == vPointsGrahamScanAlgo.at(iPointIndex + 1).GetY() 
					&& vPointsGrahamScanAlgo.at(iPointIndex).GetX() > vPointsGrahamScanAlgo.at(iPointIndex + 1).GetX()))
			{
				double dTempY = vPointsGrahamScanAlgo.at(iPointIndex + 1).GetY();
				vPointsGrahamScanAlgo.at(iPointIndex + 1).SetY(vPointsGrahamScanAlgo.at(iPointIndex).GetY());
				vPointsGrahamScanAlgo.at(iPointIndex).SetY(dTempY);

				double dTempX = vPointsGrahamScanAlgo.at(iPointIndex + 1).GetX();
				vPointsGrahamScanAlgo.at(iPointIndex + 1).SetX(vPointsGrahamScanAlgo.at(iPointIndex).GetX());
				vPointsGrahamScanAlgo.at(iPointIndex).SetX(dTempX);

				double dTempZ = vPointsGrahamScanAlgo.at(iPointIndex + 1).GetZ();
				vPointsGrahamScanAlgo.at(iPointIndex + 1).SetZ(vPointsGrahamScanAlgo.at(iPointIndex).GetZ());
				vPointsGrahamScanAlgo.at(iPointIndex).SetZ(dTempZ);
			}
		}
	}
}

bool IsPointAlreadyConsidered(std::vector<Point> vConvexHullPointsUsingGrahamScanConstructed, Point& pPoint)
{
	for (int iPointIndex = 0; iPointIndex < vConvexHullPointsUsingGrahamScanConstructed.size(); iPointIndex++)
	{
		if (pPoint == vConvexHullPointsUsingGrahamScanConstructed.at(iPointIndex))
			return true;
	}

	return false;
}

void CreateHalfConvexHullUsingGrahamScan(std::vector<Point>& vPointsGrahamScanAlgo, std::vector<Point>& vConvexHullPointsUsingGrahamScan, 
	const std::vector<Point>& vConvexHullPointsUsingGrahamScanConstructed = std::vector<Point> ())
{
	std::vector<Point> vConvexHullPointsUsingGrahamScanConstructedTemp;
	if (false == vConvexHullPointsUsingGrahamScanConstructed.empty())
		vConvexHullPointsUsingGrahamScanConstructedTemp = vConvexHullPointsUsingGrahamScanConstructed;

	for (int iPointIndex = 1; iPointIndex < vPointsGrahamScanAlgo.size(); iPointIndex++)
	{
		if (vPointsGrahamScanAlgo.at(vPointsGrahamScanAlgo.size() - 1).GetY() == vPointsGrahamScanAlgo.at(iPointIndex).GetY())
			vConvexHullPointsUsingGrahamScan.push_back(vPointsGrahamScanAlgo.at(iPointIndex));
		else
		{
			if (false == vConvexHullPointsUsingGrahamScanConstructed.empty())
			{
				if (true == IsPointAlreadyConsidered(vConvexHullPointsUsingGrahamScanConstructed, vPointsGrahamScanAlgo.at(iPointIndex)))
					continue;
			}

			double dHeadPointP1X = vConvexHullPointsUsingGrahamScan.at(vConvexHullPointsUsingGrahamScan.size() - 1).GetX() - vPointsGrahamScanAlgo.at(iPointIndex).GetX();
			double dHeadPointP1Y = vConvexHullPointsUsingGrahamScan.at(vConvexHullPointsUsingGrahamScan.size() - 1).GetY() - vPointsGrahamScanAlgo.at(iPointIndex).GetY();
			double dHeadPointP1Z = vConvexHullPointsUsingGrahamScan.at(vConvexHullPointsUsingGrahamScan.size() - 1).GetZ() - vPointsGrahamScanAlgo.at(iPointIndex).GetZ();

			bool bLastEdgeConnectionForLeftHull = false;
			int iNextPointIndex = iPointIndex + 1;
			if (false == vConvexHullPointsUsingGrahamScanConstructed.empty())
			{
				while (true == IsPointAlreadyConsidered(vConvexHullPointsUsingGrahamScanConstructed, vPointsGrahamScanAlgo.at(iNextPointIndex)))
				{
					iNextPointIndex++;
					if ((iNextPointIndex >= vPointsGrahamScanAlgo.size()))
					{
						bLastEdgeConnectionForLeftHull = true;
						break;
					}
				}
			}

			double dP2P1X, dP2P1Y, dP2P1Z;
			if (true == bLastEdgeConnectionForLeftHull)
			{
				dP2P1X = vConvexHullPointsUsingGrahamScanConstructedTemp.at(vConvexHullPointsUsingGrahamScanConstructed.size() - 1).GetX() - vPointsGrahamScanAlgo.at(iPointIndex).GetX();
				dP2P1Y = vConvexHullPointsUsingGrahamScanConstructedTemp.at(vConvexHullPointsUsingGrahamScanConstructed.size() - 1).GetY() - vPointsGrahamScanAlgo.at(iPointIndex).GetY();
				dP2P1Z = vConvexHullPointsUsingGrahamScanConstructedTemp.at(vConvexHullPointsUsingGrahamScanConstructed.size() - 1).GetZ() - vPointsGrahamScanAlgo.at(iPointIndex).GetZ();
			}
			else
			{
				dP2P1X = vPointsGrahamScanAlgo.at(iNextPointIndex).GetX() - vPointsGrahamScanAlgo.at(iPointIndex).GetX();
				dP2P1Y = vPointsGrahamScanAlgo.at(iNextPointIndex).GetY() - vPointsGrahamScanAlgo.at(iPointIndex).GetY();
				dP2P1Z = vPointsGrahamScanAlgo.at(iNextPointIndex).GetZ() - vPointsGrahamScanAlgo.at(iPointIndex).GetZ();
			}

			double dCrossProduct = (dP2P1X * dHeadPointP1Y) - (dP2P1Y * dHeadPointP1X);
			if (true == vConvexHullPointsUsingGrahamScanConstructed.empty())
			{
				if (dCrossProduct > 0)
					vConvexHullPointsUsingGrahamScan.push_back(vPointsGrahamScanAlgo.at(iPointIndex));
			}
			else if (false == vConvexHullPointsUsingGrahamScanConstructed.empty())
			{
				if (dCrossProduct < 0)
					vConvexHullPointsUsingGrahamScan.push_back(vPointsGrahamScanAlgo.at(iPointIndex));
			}

		}
	}
}

void CreateConvexHullUsingGrahamScan(std::vector<Point>& vPointsGrahamScanAlgo, std::vector<Point>& vConvexHullPointsUsingGrahamScan)
{
	// Arrange the Points in the ascsending order
	SortVectorIncreasinglyWithY(vPointsGrahamScanAlgo);

	std::vector<Point> vConvexHullPointsUsingGrahamScanRight;
	vConvexHullPointsUsingGrahamScanRight.push_back(vPointsGrahamScanAlgo.at(0));
	CreateHalfConvexHullUsingGrahamScan(vPointsGrahamScanAlgo, vConvexHullPointsUsingGrahamScanRight);

	std::vector<Point> vConvexHullPointsUsingGrahamScanLeft;
	vConvexHullPointsUsingGrahamScanLeft.push_back(vPointsGrahamScanAlgo.at(0));
	CreateHalfConvexHullUsingGrahamScan(vPointsGrahamScanAlgo, vConvexHullPointsUsingGrahamScanLeft, vConvexHullPointsUsingGrahamScanRight);

	// Combine the left and right side of convex hull
	CombineConvexHullInRightVector(vConvexHullPointsUsingGrahamScanRight, vConvexHullPointsUsingGrahamScanLeft, vConvexHullPointsUsingGrahamScan);
}

bool IsPointInsidePolyhedron(Point& pO, Point& pD, std::vector<TriangularFace>& polyhedronFaces)
{
	int iIntersectionCount = 0;
	for (int iFaceIndex = 0; iFaceIndex < polyhedronFaces.size(); iFaceIndex++)
	{
		TriangularFace tfCurrFace = polyhedronFaces.at(iFaceIndex);

		Point pV0 = tfCurrFace.GetPoint(0);
		Point pV1 = tfCurrFace.GetPoint(1);
		Point pV2 = tfCurrFace.GetPoint(2);

		// (V1 - V0)
		double dV1minusV0X = pV1.GetX() - pV0.GetX();
		double dV1minusV0Y = pV1.GetY() - pV0.GetY();
		double dV1minusV0Z = pV1.GetZ() - pV0.GetZ();

		// (V2 - V0)
		double dV2minusV0X = pV2.GetX() - pV0.GetX();
		double dV2minusV0Y = pV2.GetY() - pV0.GetY();
		double dV2minusV0Z = pV2.GetZ() - pV0.GetZ();

		// (O - V0)
		double dpOminusV2X = pO.GetX() - pV0.GetX();
		double dpOminusV2Y = pO.GetY() - pV0.GetY();
		double dpOminusV2Z = pO.GetZ() - pV0.GetZ();

		Matrix mD(3, 3);
		mD(0, 0) = -pD.GetX();
		mD(0, 1) = dV1minusV0X;
		mD(0, 2) = dV2minusV0X;
		mD(1, 0) = -pD.GetY();
		mD(1, 1) = dV1minusV0Y;
		mD(1, 2) = dV2minusV0Y;
		mD(2, 0) = -pD.GetZ();
		mD(2, 1) = dV1minusV0Z;
		mD(2, 2) = dpOminusV2Z;

		double dDetD = mD.GetDeterminantIf3X3();

		Matrix mDX(3, 3);
		mDX(0, 0) = dpOminusV2X;
		mDX(0, 1) = dV1minusV0X;
		mDX(0, 2) = dV2minusV0X;
		mDX(1, 0) = dpOminusV2Y;
		mDX(1, 1) = dV1minusV0Y;
		mDX(1, 2) = dV2minusV0Y;
		mDX(2, 0) = dpOminusV2Z;
		mDX(2, 1) = dV1minusV0Z;
		mDX(2, 2) = dpOminusV2Z;

		double dDetDX = mDX.GetDeterminantIf3X3();

		Matrix mDY(3, 3);
		mDY(0, 0) = -pD.GetX();
		mDY(0, 1) = dpOminusV2X;
		mDY(0, 2) = dV2minusV0X;
		mDY(1, 0) = -pD.GetY();
		mDY(1, 1) = dpOminusV2Y;
		mDY(1, 2) = dV2minusV0Y;
		mDY(2, 0) = -pD.GetZ();
		mDY(2, 1) = dpOminusV2Z;
		mDY(2, 2) = dpOminusV2Z;

		double dDetDY = mDY.GetDeterminantIf3X3();

		Matrix mDZ(3, 3);
		mDZ(0, 0) = -pD.GetX();
		mDZ(0, 1) = dV1minusV0X;
		mDZ(0, 2) = dpOminusV2X;
		mDZ(1, 0) = -pD.GetY();
		mDZ(1, 1) = dV1minusV0Y;
		mDZ(1, 2) = dpOminusV2Y;
		mDZ(2, 0) = -pD.GetZ();
		mDZ(2, 1) = dV1minusV0Z;
		mDZ(2, 2) = dpOminusV2Z;

		double dDetDZ = mDZ.GetDeterminantIf3X3();

		double t = dDetDX / dDetD;
		double u = dDetDY / dDetD;
		double v = dDetDZ / dDetD;

		if (u > 0 && u < 1)
			if (v > 0 and (u + v) < 1)
				if (t > 0)
					iIntersectionCount++;
	}

	if (iIntersectionCount > 0)
		return true;

	return false;
}

void FindExtremePoints(std::vector<Point>& vPointsClarksonShorAlgo, Point& pMaxX, Point& pMinX, Point& pMaxY, Point& pMinY, Point& pMaxZ, Point& pMinZ)
{
	for (int iPointIndex = 0; iPointIndex < vPointsClarksonShorAlgo.size(); iPointIndex++)
	{
		if (vPointsClarksonShorAlgo.at(iPointIndex).GetX() > pMaxX.GetX())
			pMaxX = vPointsClarksonShorAlgo.at(iPointIndex);

		if (vPointsClarksonShorAlgo.at(iPointIndex).GetX() < pMinX.GetX())
			pMinX = vPointsClarksonShorAlgo.at(iPointIndex);

		if (vPointsClarksonShorAlgo.at(iPointIndex).GetY() > pMaxY.GetY())
			pMaxY = vPointsClarksonShorAlgo.at(iPointIndex);

		if (vPointsClarksonShorAlgo.at(iPointIndex).GetY() < pMinY.GetY())
			pMinY = vPointsClarksonShorAlgo.at(iPointIndex);

		if (vPointsClarksonShorAlgo.at(iPointIndex).GetZ() > pMaxZ.GetZ())
			pMaxZ = vPointsClarksonShorAlgo.at(iPointIndex);

		if (vPointsClarksonShorAlgo.at(iPointIndex).GetZ() < pMinZ.GetZ())
			pMinZ = vPointsClarksonShorAlgo.at(iPointIndex);
	}
}

void GetFirstSecondPoints(std::vector<Point>& vExtremePoints, Point& pFirstPoint, Point& pSecondPoint)
{
	double dMaxDist = 0;

	for (int iPoint1Index = 0; iPoint1Index < vExtremePoints.size(); iPoint1Index++)
	{
		for (int iPoint2Index = 0; iPoint2Index < vExtremePoints.size(); iPoint2Index++)
		{
			double dDist = sqrt((vExtremePoints.at(iPoint2Index).GetX() - vExtremePoints.at(iPoint1Index).GetX()) * (vExtremePoints.at(iPoint2Index).GetX() - vExtremePoints.at(iPoint1Index).GetX())
				+ (vExtremePoints.at(iPoint2Index).GetY() - vExtremePoints.at(iPoint1Index).GetY()) * (vExtremePoints.at(iPoint2Index).GetY() - vExtremePoints.at(iPoint1Index).GetY())
				+ (vExtremePoints.at(iPoint2Index).GetZ() - vExtremePoints.at(iPoint1Index).GetZ()) * (vExtremePoints.at(iPoint2Index).GetZ() - vExtremePoints.at(iPoint1Index).GetZ()));

			if (dMaxDist < dDist)
			{
				dMaxDist = dDist;
				pFirstPoint = vExtremePoints.at(iPoint1Index);
				pSecondPoint = vExtremePoints.at(iPoint2Index);
			}
		}
	}
}

void GetThirdPoint(std::vector<Point>& vExtremePoints, Point& pFirstPoint, Point& pSecondPoint, Point& pThirdPoint)
{
	double dMaxDistFromLine = -std::numeric_limits<double>::infinity();

	for (int iPointIndex = 0; iPointIndex < vExtremePoints.size(); iPointIndex++)
	{
		Point pCurrPoint = vExtremePoints.at(iPointIndex);

		Vector vFirstToCurrVector;
		vFirstToCurrVector.Create3DVector(pFirstPoint, pCurrPoint);

		Vector vSecondToCurrVector;
		vSecondToCurrVector.Create3DVector(pSecondPoint, pCurrPoint);

		Vector vFirstToSecondVector;
		vFirstToSecondVector.Create3DVector(pFirstPoint, pSecondPoint);
		double dMagFirstToSecondVector = vFirstToCurrVector.GetMagnitude();

		Vector vCrossedFirstSecondCurrVectors = vFirstToCurrVector.CrossProduct3D(vSecondToCurrVector);
		double dMagCrossedFirstSecondCurrVectors = vCrossedFirstSecondCurrVectors.GetMagnitude();

		double dDistFromLine = dMagCrossedFirstSecondCurrVectors / dMagFirstToSecondVector;

		if (dDistFromLine > dMaxDistFromLine)
		{
			pThirdPoint = pCurrPoint;
			dMaxDistFromLine = dDistFromLine;
		}
	}
}

void GetApexPoint(std::vector<Point>& vPointsClarksonShorAlgo, Point& pFirstPoint, Point& pSecondPoint, Point& pThirdPoint, Point& pApexPoint)
{
	double dMaxDistFromPlane = -std::numeric_limits<double>::infinity();

	for (int iPointIndex = 0; iPointIndex < vPointsClarksonShorAlgo.size(); iPointIndex++)
	{
		Point pCurrPoint = vPointsClarksonShorAlgo.at(iPointIndex);

		Vector vFirstToSecondVector;
		vFirstToSecondVector.Create3DVector(pFirstPoint, pSecondPoint);

		Vector vFirstToThirdVector;
		vFirstToThirdVector.Create3DVector(pFirstPoint, pThirdPoint);

		Vector vNormalVector = vFirstToSecondVector.CrossProduct3D(vFirstToThirdVector);

		Vector vFirstToCurrVector;
		vFirstToCurrVector.Create3DVector(pFirstPoint, pCurrPoint);

		double dDotFirstCurrNormalVectors = vFirstToCurrVector.DotProduct3D(vNormalVector);
		double dDistFromPlane = abs(dDotFirstCurrNormalVectors) / vNormalVector.GetMagnitude();
		if (dDistFromPlane > dMaxDistFromPlane)
		{
			pApexPoint = pCurrPoint;
			dMaxDistFromPlane = dDistFromPlane;
		}
	}
}

void CreateSetOfPointsFacewise(std::vector<Point>& vPointsClarksonShorAlgo, std::vector<TriangularFace>& tfTetrahedronFaces, Point& pTetraCentroid, std::vector<Point>& vPointsNOTToConsider)
{
	for (int iTriaFaceIndex = 0; iTriaFaceIndex < tfTetrahedronFaces.size(); iTriaFaceIndex++)
	{
		TriangularFace tfCurrFace = tfTetrahedronFaces.at(iTriaFaceIndex);
		if (true == tfTetrahedronFaces.at(iTriaFaceIndex).HasSetOfPoints())
			continue;

		for (int iPointIndex = 0; iPointIndex < vPointsClarksonShorAlgo.size(); iPointIndex++)
		{
			Point pConsiderCurrPoint = vPointsClarksonShorAlgo.at(iPointIndex);
			double dDotProdForSameSide = -1;

			// Is point in vPointsNOTToConsider
			bool bPointNOTToConsider = false;
			for (int iPointNOTToConsiderIndex = 0; iPointNOTToConsiderIndex < vPointsNOTToConsider.size(); iPointNOTToConsiderIndex++)
			{
				if (pConsiderCurrPoint == vPointsNOTToConsider.at(iPointNOTToConsiderIndex))
				{
					bPointNOTToConsider = true;
					break;
				}
			}
			if (true == bPointNOTToConsider)
				continue;

			Point pConsiderFirst = tfCurrFace.GetPoint(0);
			Point pConsiderSecond = tfCurrFace.GetPoint(1);
			Point pConsiderThird = tfCurrFace.GetPoint(2);

			Vector vConsiderFirstToSecondVector;
			vConsiderFirstToSecondVector.Create3DVector(pConsiderFirst, pConsiderSecond);

			Vector vConsiderFirstToThirdVector;
			vConsiderFirstToThirdVector.Create3DVector(pConsiderFirst, pConsiderThird);

			// Use the centroid of the tetrahedron from args
			// Create a vector from centroid to a point of the face
			// Dot product vector (centroid to face pt) with normal vector
			// if negative dot product, reverse the direction of the normal vector
			// reversing means adding minus (-) sign to each coordinate
			Vector vConsiderNormalVector = vConsiderFirstToSecondVector.CrossProduct3D(vConsiderFirstToThirdVector);
			Vector vCentroidToFirstVector;
			vCentroidToFirstVector.Create3DVector(pTetraCentroid, pConsiderFirst);
			double dNormalCentroidSide = vCentroidToFirstVector.DotProduct3D(vConsiderNormalVector);
			if (dNormalCentroidSide < 0)
				vConsiderNormalVector.Reverse();

			Vector vConsiderFirstToCurrVector;
			vConsiderFirstToCurrVector.Create3DVector(pConsiderFirst, pConsiderCurrPoint);

			dDotProdForSameSide = vConsiderFirstToCurrVector.DotProduct3D(vConsiderNormalVector);

			if (dDotProdForSameSide > 0)
			{
				tfTetrahedronFaces.at(iTriaFaceIndex).AddPointInSet(pConsiderCurrPoint);
				vPointsNOTToConsider.push_back(pConsiderCurrPoint);
			}
		}

		tfTetrahedronFaces.at(iTriaFaceIndex).SetHasSetOfPoints(true);
	}
}

Point GetTetraCentroid(Tetrahedron tTetrahedron)
{
	std::vector<TriangularFace> vTetraFaces = tTetrahedron.GetFaces();

	double dTotalX = 0;
	double dTotalY = 0;
	double dTotalZ = 0;
	for (int iTetraFaceIndex = 0; iTetraFaceIndex < vTetraFaces.size(); iTetraFaceIndex++)
	{
		for (int iFaceVertex = 0; iFaceVertex < 3; iFaceVertex++)
		{
			dTotalX += vTetraFaces.at(iTetraFaceIndex).GetPoint(iFaceVertex).GetX();
			dTotalY += vTetraFaces.at(iTetraFaceIndex).GetPoint(iFaceVertex).GetY();
			dTotalZ += vTetraFaces.at(iTetraFaceIndex).GetPoint(iFaceVertex).GetZ();
		}
	}

	double dAvgX = dTotalX / (vTetraFaces.size() * 3);
	double dAvgY = dTotalY / (vTetraFaces.size() * 3);
	double dAvgZ = dTotalZ / (vTetraFaces.size() * 3);

	return Point(dAvgX, dAvgY, dAvgZ);
}

void Create3DConvexHullUsingQuickHull(std::vector<Point>& vPointsClarksonShorAlgo, std::vector<TriangularFace>& tfConvexHullFaces)
{
	double dPositiveInf = std::numeric_limits<double>::infinity();
	double dNegativeInf = -std::numeric_limits<double>::infinity();

	Point pMaxX(dNegativeInf, dNegativeInf, dNegativeInf);
	Point pMinX(dPositiveInf, dPositiveInf, dPositiveInf);
	Point pMaxY(dNegativeInf, dNegativeInf, dNegativeInf);
	Point pMinY(dPositiveInf, dPositiveInf, dPositiveInf);
	Point pMaxZ(dNegativeInf, dNegativeInf, dNegativeInf);
	Point pMinZ(dPositiveInf, dPositiveInf, dPositiveInf);
	FindExtremePoints(vPointsClarksonShorAlgo, pMaxX, pMinX, pMaxY, pMinY, pMaxZ, pMinZ);

	std::vector<Point> vExtremePoints;
	vExtremePoints.push_back(pMaxX);
	vExtremePoints.push_back(pMinX);
	vExtremePoints.push_back(pMaxY);
	vExtremePoints.push_back(pMinY);
	vExtremePoints.push_back(pMaxZ);
	vExtremePoints.push_back(pMinZ);

	Point pFirstPoint;
	Point pSecondPoint;
	GetFirstSecondPoints(vExtremePoints, pFirstPoint, pSecondPoint);

	// Points not to be considered
	std::vector<Point> vPointsNOTToConsider;
	vPointsNOTToConsider.push_back(pFirstPoint);
	vPointsNOTToConsider.push_back(pSecondPoint);

	// Get Third Point
	Point pThirdPoint;
	GetThirdPoint(vExtremePoints, pFirstPoint, pSecondPoint, pThirdPoint);

	vPointsNOTToConsider.push_back(pThirdPoint);

	// Get Apex of Tetrahedron
	Point pApexPoint;
	GetApexPoint(vPointsClarksonShorAlgo, pFirstPoint, pSecondPoint, pThirdPoint, pApexPoint);

	vPointsNOTToConsider.push_back(pApexPoint);

	// Create 4 faces
	TriangularFace tf123(Point(pFirstPoint.GetX(), pFirstPoint.GetY(), pFirstPoint.GetZ()),
		Point(pSecondPoint.GetX(), pSecondPoint.GetY(), pSecondPoint.GetZ()),
		Point(pThirdPoint.GetX(), pThirdPoint.GetY(), pThirdPoint.GetZ()));

	TriangularFace tf12Apex(Point(pFirstPoint.GetX(), pFirstPoint.GetY(), pFirstPoint.GetZ()),
		Point(pSecondPoint.GetX(), pSecondPoint.GetY(), pSecondPoint.GetZ()),
		Point(pApexPoint.GetX(), pApexPoint.GetY(), pApexPoint.GetZ()));

	TriangularFace tf1Apex3(Point(pFirstPoint.GetX(), pFirstPoint.GetY(), pFirstPoint.GetZ()),
		Point(pApexPoint.GetX(), pApexPoint.GetY(), pApexPoint.GetZ()),
		Point(pThirdPoint.GetX(), pThirdPoint.GetY(), pThirdPoint.GetZ()));

	TriangularFace tfApex23(Point(pApexPoint.GetX(), pApexPoint.GetY(), pApexPoint.GetZ()),
		Point(pSecondPoint.GetX(), pSecondPoint.GetY(), pSecondPoint.GetZ()),
		Point(pThirdPoint.GetX(), pThirdPoint.GetY(), pThirdPoint.GetZ()));

	tfConvexHullFaces.push_back(tf123);
	tfConvexHullFaces.push_back(tf12Apex);
	tfConvexHullFaces.push_back(tf1Apex3);
	tfConvexHullFaces.push_back(tfApex23);

	// Find the centroid of the tetrahedron
	Tetrahedron tMainTEtrahedron(tf123, tf12Apex, tf1Apex3, tfApex23);
	Point pMainTetraCentroid = GetTetraCentroid(tMainTEtrahedron);

	// Create facewise set of points
	CreateSetOfPointsFacewise(vPointsClarksonShorAlgo, tfConvexHullFaces, pMainTetraCentroid, vPointsNOTToConsider);

	// Consider each face with its points' set and create convex hull faces vector.
	// Find next apex if still there are some points outside the new tetrahedrons continue finding apex points until the nummber is zero.
	// Create new vectors for the points not to be considered for each face separately.
	for (int iTetrahedronFaceIndex = 0; iTetrahedronFaceIndex < tfConvexHullFaces.size(); iTetrahedronFaceIndex++)
	{
		TriangularFace tfCurrFace = tfConvexHullFaces.at(iTetrahedronFaceIndex);
		std::vector<Point> vCurrFaceSetOfPoints = tfCurrFace.GetSetOfPoints();
		if (0 == vCurrFaceSetOfPoints.size())
			continue;

		std::vector<Point> vPointsNOTToConsiderAgain;

		double dMaxDistFromPlane = -std::numeric_limits<double>::infinity();
		
		Point pFirstPoint = tfCurrFace.GetPoint(0);
		Point pSecondPoint = tfCurrFace.GetPoint(1);
		Point pThirdPoint = tfCurrFace.GetPoint(2);

		Point p2ndaryApexPoint;
		/*for (int iPointIndex = 0; iPointIndex < vCurrFaceSetOfPoints.size(); iPointIndex++)
		{
			Point pCurrPoint = vCurrFaceSetOfPoints.at(iPointIndex);

			double dDotProdForSameSide = -1;


			
			Vector vFirstToSecondVector;
			vFirstToSecondVector.Create3DVector(pFirstPoint, pSecondPoint);

			Vector vFirstToThirdVector;
			vFirstToThirdVector.Create3DVector(pFirstPoint, pThirdPoint);

			Vector vNormalVector = vFirstToSecondVector.CrossProduct3D(vFirstToThirdVector);

			Vector vFirstToCurrVector;
			vFirstToCurrVector.Create3DVector(pFirstPoint, pCurrPoint);

			double dDotFirstCurrNormalVectors = vFirstToCurrVector.DotProduct3D(vNormalVector);
			double dDistFromPlane = abs(dDotFirstCurrNormalVectors) / vNormalVector.GetMagnitude();
			if (dDistFromPlane > dMaxDistFromPlane)
			{
				p2ndaryApexPoint = pCurrPoint;
				dMaxDistFromPlane = dDistFromPlane;
			}
		}*/
		GetApexPoint(vCurrFaceSetOfPoints, pFirstPoint, pSecondPoint, pThirdPoint, p2ndaryApexPoint);

		TriangularFace tf122ndaryApex(Point(pFirstPoint.GetX(), pFirstPoint.GetY(), pFirstPoint.GetZ()),
			Point(pSecondPoint.GetX(), pSecondPoint.GetY(), pSecondPoint.GetZ()),
			Point(p2ndaryApexPoint.GetX(), p2ndaryApexPoint.GetY(), p2ndaryApexPoint.GetZ()));

		TriangularFace tf12ndaryApex3(Point(pFirstPoint.GetX(), pFirstPoint.GetY(), pFirstPoint.GetZ()),
			Point(p2ndaryApexPoint.GetX(), p2ndaryApexPoint.GetY(), p2ndaryApexPoint.GetZ()),
			Point(pThirdPoint.GetX(), pThirdPoint.GetY(), pThirdPoint.GetZ()));

		TriangularFace tf2ndaryApex23(Point(p2ndaryApexPoint.GetX(), p2ndaryApexPoint.GetY(), p2ndaryApexPoint.GetZ()),
			Point(pSecondPoint.GetX(), pSecondPoint.GetY(), pSecondPoint.GetZ()),
			Point(pThirdPoint.GetX(), pThirdPoint.GetY(), pThirdPoint.GetZ()));

		tfConvexHullFaces.push_back(tf122ndaryApex);
		tfConvexHullFaces.push_back(tf12ndaryApex3);
		tfConvexHullFaces.push_back(tf2ndaryApex23);

		// Find the centroid of the tetrahedron
		Tetrahedron tRestTetrahedron(tf122ndaryApex, tf12ndaryApex3, tf2ndaryApex23, tfConvexHullFaces.at(iTetrahedronFaceIndex));
		Point pRestTetraCentroid = GetTetraCentroid(tRestTetrahedron);

		// Set Points' set for new faces
		/*for (int iTriaFaceIndex = 0; iTriaFaceIndex < tfConvexHullFaces.size(); iTriaFaceIndex++)
		{
			TriangularFace tfCurrFace = tfConvexHullFaces.at(iTriaFaceIndex);
			if (true == tfCurrFace.HasSetOfPoints())
				continue;

			for (int iPointIndex = 0; iPointIndex < vCurrFaceSetOfPoints.size(); iPointIndex++)
			{
				Point pConsiderCurrPoint = vCurrFaceSetOfPoints.at(iPointIndex);
				double dDotProdForSameSide = -1;

				bool bPointNOTToConsider = false;
				for (int iPointNOTToConsiderIndex = 0; iPointNOTToConsiderIndex < vPointsNOTToConsiderAgain.size(); iPointNOTToConsiderIndex++)
				{
					if (pConsiderCurrPoint == vPointsNOTToConsiderAgain.at(iPointNOTToConsiderIndex))
					{
						bPointNOTToConsider = true;
						break;
					}
				}
				if (true == bPointNOTToConsider)
					continue;

				Point pConsiderFirst = tfCurrFace.GetPoint(0);
				Point pConsiderSecond = tfCurrFace.GetPoint(1);
				Point pConsiderThird = tfCurrFace.GetPoint(2);

				Vector vConsiderFirstToSecondVector;
				vConsiderFirstToSecondVector.Create3DVector(pConsiderFirst, pConsiderSecond);

				Vector vConsiderFirstToThirdVector;
				vConsiderFirstToThirdVector.Create3DVector(pConsiderFirst, pConsiderThird);

				Vector vConsiderNormalVector = vConsiderFirstToSecondVector.CrossProduct3D(vConsiderFirstToThirdVector);

				Vector vConsiderFirstToCurrVector;
				vConsiderFirstToCurrVector.Create3DVector(pConsiderFirst, pConsiderCurrPoint);

				dDotProdForSameSide = vConsiderFirstToCurrVector.DotProduct3D(vConsiderNormalVector);

				if (dDotProdForSameSide > 0)
				{
					tfCurrFace.AddPointInSet(pConsiderCurrPoint);
					vPointsNOTToConsiderAgain.push_back(pConsiderCurrPoint);
				}
			}

			tfCurrFace.SetHasSetOfPoints(true);
		}*/
		CreateSetOfPointsFacewise(vCurrFaceSetOfPoints, tfConvexHullFaces, pMainTetraCentroid, vPointsNOTToConsiderAgain);
	}

	int i = 0;
}

bool IsPointWithinSegmentRange(Point& pIntersectionPoint, Line& line, bool bIs3D)
{
	if (true == bIs3D)
	{
		if (pIntersectionPoint.GetX() > std::min(line.GetStart().GetX(), line.GetEnd().GetX()) &&
			pIntersectionPoint.GetY() > std::min(line.GetStart().GetY(), line.GetEnd().GetY()) &&
			pIntersectionPoint.GetZ() > std::min(line.GetStart().GetZ(), line.GetEnd().GetZ()))
		{
			return true;
		}
	}
	else
	{
		if (pIntersectionPoint.GetX() > std::min(line.GetStart().GetX(), line.GetEnd().GetX()) &&
			pIntersectionPoint.GetY() > std::min(line.GetStart().GetY(), line.GetEnd().GetY()))
		{
			return true;
		}
	}
	
	return false;
}

bool GetLinesIntersectionPoint(Line line1, Line line2, bool bIs3D, Point& pIntersectionPoint)
{
	// Get Syemmteric form of 2 lines
	// (x - x1) / (x2 - x1) = (y - y1) / (y2 - y1) = (z - z1) / (z2 - z1) = t
	// x = ((x2 - x1) * t) + x1; y = ((y2 - y1) * t) + y1; z = ((z2 - z1) * t) + z1 ______ (equn)
	// x = ((x4 - x3) * s) + x3;
	// t * (x2 - x1) + x1 = s * (x4 - x3) + x3;
	// t * (y2 - y1) + y1 = s * (y4 - y3) + y3;

	double x1 = line1.GetStart().GetX();
	double y1 = line1.GetStart().GetY();
	double z1 = line1.GetStart().GetZ();

	double x2 = line1.GetEnd().GetX();
	double y2 = line1.GetEnd().GetY();
	double z2 = line1.GetEnd().GetZ();

	double x3 = line2.GetStart().GetX();
	double y3 = line2.GetStart().GetY();
	double z3 = line2.GetStart().GetZ();

	double x4 = line2.GetEnd().GetX();
	double y4 = line2.GetEnd().GetY();
	double z4 = line2.GetEnd().GetZ();

	double t = ((x3 - x1) * (y4 - y3) + (y1 - y3) * (x4 - x3)) /
		((x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1));
	if (std::numeric_limits<double>::infinity() == t)
		return false;
	// Put t in (equn)

	double s = (t * (y2 - y1) + y1 - y3) /
		(y4 - y3);
	if (std::numeric_limits<double>::infinity() == s)
		return false;

	if (true == bIs3D)
	{
		// Check validity. s and t should satisfy all the equns involving s and t i.e. 3 equns (2 aare used for finding s and t).
		double dValidzLHS = (static_cast<int>((t * (z2 - z1) + z1) * 100)) / 100;
		double dValidzRHS = (static_cast<int>((s * (z4 - z3) + z3) * 100)) / 100;

		if (dValidzLHS != dValidzRHS)
			return false;
	}
		
	pIntersectionPoint.SetX((t * (x2 - x1) + x1));
	pIntersectionPoint.SetY((t * (y2 - y1) + y1));
	pIntersectionPoint.SetZ((t * (z2 - z1) + z1));

	bool bPointWithinline1Range = false;
	if (true == line1.IsSegment())
		bPointWithinline1Range = IsPointWithinSegmentRange(pIntersectionPoint, line1, bIs3D);

	bool bPointWithinline2Range = false;
	if (true == line2.IsSegment())
		bPointWithinline2Range = IsPointWithinSegmentRange(pIntersectionPoint, line2, bIs3D);

	if (false == line1.IsSegment() && false == line2.IsSegment())
		return true;
	else if (true == line1.IsSegment() && false == line2.IsSegment())
	{
		if (true == bPointWithinline1Range)
			return true;
		else
			return false;
	}
	else if (false == line1.IsSegment() && true == line2.IsSegment())
	{
		if (true == bPointWithinline2Range)
			return true;
		else
			return false;
	}
	else if (true == line1.IsSegment() && true == line2.IsSegment())
	{
		if (true == bPointWithinline1Range && bPointWithinline2Range)
			return true;
		else
			return false;
	}

	return false;
}

Point GetPlanesIntersectionCoordinates(Plane& plane1, Plane& plane2, double zCoord)
{
	double dPlane1A = plane1.GetA();
	double dPlane1B = plane1.GetB();
	double dPlane1C = plane1.GetC();
	double dPlane1D = plane1.GetD();

	double dPlane2A = plane2.GetA();
	double dPlane2B = plane2.GetB();
	double dPlane2C = plane2.GetC();
	double dPlane2D = plane2.GetD();

	double y = ((dPlane2D - dPlane2D) - (zCoord) * (dPlane2C - dPlane1C)) / (dPlane2B - dPlane1B);
	double x = (-(dPlane1B * y) - (dPlane1C * zCoord) + dPlane1D) / dPlane1A;

	return Point(x, y, zCoord);
}

bool ArePlanesParallel(Plane plane1, Plane plane2)
{
	std::vector<double> vPlane1VectorCompo;
	vPlane1VectorCompo.push_back(plane1.GetA());
	vPlane1VectorCompo.push_back(plane1.GetB());
	vPlane1VectorCompo.push_back(plane1.GetC());
	Vector vPLane1NormalVector(vPlane1VectorCompo);

	std::vector<double> vPlane2VectorCompo;
	vPlane2VectorCompo.push_back(plane2.GetA());
	vPlane2VectorCompo.push_back(plane2.GetB());
	vPlane2VectorCompo.push_back(plane2.GetC());
	Vector vPLane2NormalVector(vPlane2VectorCompo);

	Vector vCrossedPlanes = vPLane1NormalVector.CrossProduct3D(vPLane2NormalVector);
	if (0 == vCrossedPlanes.GetMagnitude())
		return true;

	return false;
}

bool GetPlanesIntersetionLine(Plane& plane1, Plane& plane2, Line& lIntersectionLine)
{
	double dPlane1A = plane1.GetA();
	double dPlane1B = plane1.GetB();
	double dPlane1C = plane1.GetC();

	double dPlane2A = plane2.GetA();
	double dPlane2B = plane2.GetB();
	double dPlane2C = plane2.GetC();

	
	bool bArePlanesParallel = ArePlanesParallel(plane1, plane2);
	if (true == bArePlanesParallel)
		return false;

	// Eliminate x
	plane2.SetA(dPlane2A * dPlane1A / dPlane2A);
	plane2.SetB(dPlane2B * dPlane1A / dPlane2A);
	plane2.SetC(dPlane2C * dPlane1A / dPlane2A);

	Point pStartIntersectionLine = GetPlanesIntersectionCoordinates(plane1, plane2, 0);
	Point pEndIntersectionLine = GetPlanesIntersectionCoordinates(plane1, plane2, 1);

	lIntersectionLine = Line(pStartIntersectionLine, pEndIntersectionLine);
	
	return true;
}

void FindExtremePointsForSuperTriangle(std::vector<Node>& vPointsDelaunayAlgo, Node& pMaxX, Node& pMinX, Node& pMaxY, Node& pMinY)
{
	for (int iPointIndex = 0; iPointIndex < vPointsDelaunayAlgo.size(); iPointIndex++)
	{
		if (vPointsDelaunayAlgo.at(iPointIndex).GetX() > pMaxX.GetX())
			pMaxX = vPointsDelaunayAlgo.at(iPointIndex);

		if (vPointsDelaunayAlgo.at(iPointIndex).GetX() < pMinX.GetX())
			pMinX = vPointsDelaunayAlgo.at(iPointIndex);

		if (vPointsDelaunayAlgo.at(iPointIndex).GetY() > pMaxY.GetY())
			pMaxY = vPointsDelaunayAlgo.at(iPointIndex);

		if (vPointsDelaunayAlgo.at(iPointIndex).GetY() < pMinY.GetY())
			pMinY = vPointsDelaunayAlgo.at(iPointIndex);
	}
}

Tria CreateSuperTriangle(std::vector<Node>& vPointsDelaunayAlgo, double& dAllPtsInsideSurityConst)
{
	double dPositiveInf = std::numeric_limits<double>::infinity();
	double dNegativeInf = -std::numeric_limits<double>::infinity();

	Node pMaxX(dNegativeInf, dNegativeInf, dNegativeInf);
	Node pMinX(dPositiveInf, dPositiveInf, dPositiveInf);
	Node pMaxY(dNegativeInf, dNegativeInf, dNegativeInf);
	Node pMinY(dPositiveInf, dPositiveInf, dPositiveInf);

	FindExtremePointsForSuperTriangle(vPointsDelaunayAlgo, pMaxX, pMinX, pMaxY, pMinY);
	if (0 == pMinX.GetX())
		pMinX.SetX(1);

	if (0 == pMinX.GetY())
		pMinX.SetY(1);

	if (0 == pMaxX.GetX())
		pMaxX.SetX(-1);

	if (0 == pMaxX.GetY())
		pMaxX.SetY(-1);

	if (0 == pMinY.GetX())
		pMinY.SetX(1);

	if (0 == pMinY.GetY())
		pMinY.SetY(1);

	if (0 == pMaxY.GetX())
		pMaxY.SetX(-1);

	if (0 == pMaxY.GetY())
		pMaxY.SetY(-1);

	// Create vertex of Super Triangle
	Node n1(pMinX.GetX() - (std::abs(pMinX.GetX()) * dAllPtsInsideSurityConst), pMinY.GetY() - (std::abs(pMinY.GetY()) * dAllPtsInsideSurityConst), 0.0);
	Node n2(pMaxX.GetX() + (std::abs(pMaxX.GetX()) * dAllPtsInsideSurityConst), pMinY.GetY() - (std::abs(pMinY.GetY()) * dAllPtsInsideSurityConst), 0.0);
	Node n3(((pMinX.GetX() + pMaxX.GetX()) / 2), pMaxY.GetY()  + (std::abs(pMaxY.GetY()) * dAllPtsInsideSurityConst), 0.0);

	Tria tSuperTria(Edge2(n1, n2), Edge2(n2, n3), Edge2(n3, n1));

	return tSuperTria;
}

bool AreNodesInside(std::vector<Node>& vPointsDelaunayAlgo, Tria& tSuperTria)
{
	Node P1 = tSuperTria.P1();
	Node P2 = tSuperTria.P2();
	Node P3 = tSuperTria.P3();
	std::vector<Node> vNodeVec;
	vNodeVec.push_back(P1);
	vNodeVec.push_back(P2);
	vNodeVec.push_back(P3);

	bool bIsNodeInside = true;
	for (int iPointIndex = 0; iPointIndex < vPointsDelaunayAlgo.size(); iPointIndex++)
	{
		Node pCurrNode = vPointsDelaunayAlgo.at(iPointIndex);

		double dPreviousCrossProd = 0;
		for (int iCount = 0; iCount < vNodeVec.size(); iCount++)
		{
			Node pCurrTriaNode = vNodeVec.at(iCount);
			Node pCurrTriaNextNode;
			if ((vNodeVec.size() - 1) == iCount)
				pCurrTriaNextNode = vNodeVec.at(0);
			else
				pCurrTriaNextNode = vNodeVec.at(iCount+1);

			Vector vCurrPStrt;
			vCurrPStrt.Create3DVector(pCurrNode, pCurrTriaNode);

			Vector vPEndPStrt;
			vPEndPStrt.Create3DVector(pCurrTriaNextNode, pCurrTriaNode);

			double dInsideOrOutside = vPEndPStrt.CrossProduct2D(vCurrPStrt);
			if (0 == iCount)
			{
				dPreviousCrossProd = dInsideOrOutside;
			}
			else if ((dInsideOrOutside > 0 && dPreviousCrossProd > 0) || (dInsideOrOutside < 0 && dPreviousCrossProd < 0))
			{
				//dPreviousCrossProd = dInsideOrOutside;
			}
			else
			{
				bIsNodeInside = false;
				break;
			}
		}
	}

	return bIsNodeInside;
}

bool IsNodeInsideElement(Tria& pCurrElement, Node& pCurrNode)
{
	bool bIsPointInside = true;

	Node P1 = pCurrElement.P1();
	Node P2 = pCurrElement.P2();
	Node P3 = pCurrElement.P3();
	std::vector<Node> vNodeVec;
	vNodeVec.push_back(P1);
	vNodeVec.push_back(P2);
	vNodeVec.push_back(P3);

	double dPreviousCrossProd = NULL;
	for (int iCount = 0; iCount < vNodeVec.size(); iCount++)
	{
		Node pCurrTriaNode = vNodeVec.at(iCount);
		Node pCurrTriaNextNode;
		if ((vNodeVec.size() - 1) == iCount)
			pCurrTriaNextNode = vNodeVec.at(0);
		else
			pCurrTriaNextNode = vNodeVec.at(iCount + 1);

		Vector vCurrPStrt;
		vCurrPStrt.Create3DVector(pCurrNode, pCurrTriaNode);
		vCurrPStrt.Normalize();

		Vector vPEndPStrt;
		vPEndPStrt.Create3DVector(pCurrTriaNextNode, pCurrTriaNode);
		vPEndPStrt.Normalize();

		double dInsideOrOutside = vPEndPStrt.CrossProduct2D(vCurrPStrt);
		if (0 == iCount)
		{
			dPreviousCrossProd = dInsideOrOutside;
		}
		else if ((dInsideOrOutside > 0 && dPreviousCrossProd > 0) || (dInsideOrOutside < 0 && dPreviousCrossProd < 0))
		{
			//dPreviousCrossProd = dInsideOrOutside;
		}
		else
		{
			bIsPointInside = false;
			break;
		}
	}

	return bIsPointInside;
}

Point GetMidPtOfNodes(Node& n1, Node& n2)
{
	double xMidCoord = (n1.GetX() + n2.GetX()) / 2;
	double yMidCoord = (n1.GetY() + n2.GetY()) / 2;
	double zMidCoord = (n1.GetZ() + n2.GetZ()) / 2;

	return Point(xMidCoord, yMidCoord, zMidCoord);
}

Line GetPerpendicularBisectorOfElemEdge(Node& n1, Node& n2, Point& pMidPoint)
{
	double n1X = n1.GetX();
	double n1Y = n1.GetY();

	double n2X = n2.GetX();
	double n2Y = n2.GetY();

	double dSlopeOfEdge = (n2Y - n1Y) / (n2X - n1X);
	double dSlopeOfBisector = -(1 / dSlopeOfEdge);

	Point pAlongBisector(pMidPoint.GetX() + 50, pMidPoint.GetY() + (50 * dSlopeOfBisector), 0);

	return Line(pMidPoint, pAlongBisector);
}

double GetDistBetnPoints(Point& pCCCenter, Node& node)
{
	double dDist = sqrt((node.GetX() - pCCCenter.GetX()) * (node.GetX() - pCCCenter.GetX())
		+ (node.GetY() - pCCCenter.GetY()) * (node.GetY() - pCCCenter.GetY())
		+ (node.GetZ() - pCCCenter.GetZ()) * (node.GetZ() - pCCCenter.GetZ()));

	return dDist;
}

bool IsNodeInsideTriaCC(Tria& pNeighbourElem, Node& pCurrNode)
{
	Node P1 = pNeighbourElem.P1();
	Node P2 = pNeighbourElem.P2();
	Node P3 = pNeighbourElem.P3();

	if (pCurrNode == P1 || pCurrNode == P2 || pCurrNode == P3)
		return false;

	Point pMidP1P2 = GetMidPtOfNodes(P1, P2);
	Point pMidP2P3 = GetMidPtOfNodes(P2, P3);

	Line lP1P2Bisector = GetPerpendicularBisectorOfElemEdge(P1, P2, pMidP1P2);
	Line lP2P3Bisector = GetPerpendicularBisectorOfElemEdge(P2, P3, pMidP2P3);

	Point pCCCenter;
	GetLinesIntersectionPoint(lP1P2Bisector, lP2P3Bisector, true, pCCCenter);

	double dCCRadius = GetDistBetnPoints(pCCCenter, P1);

	Circle cElemCircumCircle(pCCCenter, dCCRadius);
	Point pCurrPoint(pCurrNode.GetX(), pCurrNode.GetY(), pCurrNode.GetZ());
	bool bOnCCEdge = false;
	bool bIsNodeInsideTriaCC = cElemCircumCircle.IsPointInside(pCurrPoint, bOnCCEdge);

	if (true == bOnCCEdge && true == bIsNodeInsideTriaCC)
		return false;

	return bIsNodeInsideTriaCC;
}

int GetNewTriaHavingCommonEdge(std::vector<Tria>& vElements, Node& pCommonNode1, Node& pCommonNode2)
{
	int iTriaHavingCommonEdge = -1;

	for (int iTriaIndex = 0; iTriaIndex < vElements.size(); iTriaIndex++)
	{
		if (true == vElements.at(iTriaIndex).GetToBeDeleted())
			continue;

		Node P1 = vElements.at(iTriaIndex).P1();
		Node P2 = vElements.at(iTriaIndex).P2();
		Node P3 = vElements.at(iTriaIndex).P3();

		bool bCommonNode1Found = false;
		if (pCommonNode1 == P1 || pCommonNode1 == P2 || pCommonNode1 == P3)
			bCommonNode1Found = true;

		bool bCommonNode2Found = false;
		if (pCommonNode2 == P1 || pCommonNode2 == P2 || pCommonNode2 == P3)
			bCommonNode2Found = true;

		if (true == bCommonNode1Found && true == bCommonNode2Found)
		{
			iTriaHavingCommonEdge = iTriaIndex;
			break;
		}
	}

	return iTriaHavingCommonEdge;
}

void DeleteSuperTriaConnectedElements(std::vector<Tria>& vElements, std::vector<Node> vSuperTriaNodes)
{
	for (int iElemIndex = 0; iElemIndex < vElements.size(); iElemIndex++)
	{
		Node P1 = vElements.at(iElemIndex).P1();
		Node P2 = vElements.at(iElemIndex).P2();
		Node P3 = vElements.at(iElemIndex).P3();

		for (int iSuperTriaNodeIndex = 0; iSuperTriaNodeIndex < vSuperTriaNodes.size(); iSuperTriaNodeIndex++)
		{
			Node nCurrNode = vSuperTriaNodes.at(iSuperTriaNodeIndex);

			if (nCurrNode == P1 || nCurrNode == P2 || nCurrNode == P3)
				vElements.at(iElemIndex).SetToBeDeleted(true);
		}
	}
}

double CalculateAreaUsingTriaElems(std::vector<Tria> vElements)
{
	double dTotalArea = 0;

	for (int iElemIndex = 0; iElemIndex < vElements.size(); iElemIndex++)
	{
		Node P1 = vElements.at(iElemIndex).P1();
		Node P2 = vElements.at(iElemIndex).P2();
		Node P3 = vElements.at(iElemIndex).P3();

		Vector vP1P2;
		vP1P2.Create3DVector(P2, P1);

		Vector vP3P2;
		vP3P2.Create3DVector(P2, P3);

		double dCrossProduct = vP1P2.CrossProduct2D(vP3P2);
		double dTriaArea = std::abs(dCrossProduct) / 2;
		dTotalArea += dTriaArea;
	}

	return dTotalArea;
}

bool MeshUsingDelaunayAlgo(std::vector<Node>& vPointsDelaunayAlgo, std::vector<Tria>& vElements, double& dPolygonArea)
{
	if (vPointsDelaunayAlgo.size() < 4)
		return false;

	double dAllPtsInsideSurityConst = 10;
	Tria tSuperTria = CreateSuperTriangle(vPointsDelaunayAlgo, dAllPtsInsideSurityConst);
	
	// Check all points are inside the Super triangle or not
	while (false == AreNodesInside(vPointsDelaunayAlgo, tSuperTria))
	{
		dAllPtsInsideSurityConst = dAllPtsInsideSurityConst * 10;
		tSuperTria = CreateSuperTriangle(vPointsDelaunayAlgo, dAllPtsInsideSurityConst);
	}

	vElements.push_back(tSuperTria);

	std::vector<Node> vSuperTriaNodes;
	vSuperTriaNodes.push_back(tSuperTria.P1());
	vSuperTriaNodes.push_back(tSuperTria.P2());
	vSuperTriaNodes.push_back(tSuperTria.P3());

	for (int iPointIndex = 0; iPointIndex < vPointsDelaunayAlgo.size(); iPointIndex++)
	{
		Node pCurrNode = vPointsDelaunayAlgo.at(iPointIndex);

		Tria pCurrElem;
		bool bIsNodeInsideElem = false;
		int iElementIndexForDeletionMark = -1;
		for (int iElementIndex = 0; iElementIndex < vElements.size(); iElementIndex++)
		{
			pCurrElem = vElements.at(iElementIndex);
			if (true == pCurrElem.GetToBeDeleted())
				continue;

			bIsNodeInsideElem = IsNodeInsideElement(pCurrElem, pCurrNode);
			if (true == bIsNodeInsideElem)
			{
				iElementIndexForDeletionMark = iElementIndex;
				break;
			}
		}

		if (false == bIsNodeInsideElem)
			continue;

		// Create 3 elems in the current Tria using current Node
		if (-1 != iElementIndexForDeletionMark)
			vElements.at(iElementIndexForDeletionMark).SetToBeDeleted(true);

		Node P1 = pCurrElem.P1();
		Node P2 = pCurrElem.P2();
		Node P3 = pCurrElem.P3();

		Tria tTria1(Edge2(pCurrNode, P1), Edge2(P1, P2), Edge2(P2, pCurrNode));
		vElements.push_back(tTria1);

		Tria tTria2(Edge2(pCurrNode, P2), Edge2(P2, P3), Edge2(P3, pCurrNode));
		vElements.push_back(tTria2);

		Tria tTria3(Edge2(pCurrNode, P3), Edge2(P3, P1), Edge2(P1, pCurrNode));
		vElements.push_back(tTria3);

		for (int iElementIndex = 0; iElementIndex < vElements.size(); iElementIndex++)
		{
			Tria pCurrNeighbourElem = vElements.at(iElementIndex);
			if (true == pCurrNeighbourElem.GetToBeDeleted())
				continue;

			bool bIsNodeInsideTriaCC = IsNodeInsideTriaCC(pCurrNeighbourElem, pCurrNode);
			if (false == bIsNodeInsideTriaCC)
				continue;

			Node pCommonNode1;
			Node pCommonNode2;
			bool bGetCommonNodes = pCurrElem.GetCommonNodes(pCurrNeighbourElem, pCommonNode1, pCommonNode2);
			if (false == bGetCommonNodes)
				continue;

			Node pUncommonNode;
			Node nNeighbourElemNodeP1 = pCurrNeighbourElem.P1();
			Node nNeighbourElemNodeP2 = pCurrNeighbourElem.P2();
			Node nNeighbourElemNodeP3 = pCurrNeighbourElem.P3();
			std::vector<Node> vNodeVecNeighbourElem;
			vNodeVecNeighbourElem.push_back(nNeighbourElemNodeP1);
			vNodeVecNeighbourElem.push_back(nNeighbourElemNodeP2);
			vNodeVecNeighbourElem.push_back(nNeighbourElemNodeP3);

			for (int iCount = 0; iCount < vNodeVecNeighbourElem.size(); iCount++)
			{
				Node nCurrNeighbourNode = vNodeVecNeighbourElem.at(iCount);

				if (nCurrNeighbourNode == pCommonNode1 || nCurrNeighbourNode == pCommonNode2)
					continue;

				pUncommonNode = nCurrNeighbourNode;
			}

			vElements.at(iElementIndex).SetToBeDeleted(true);
			// Get New Tria index for deletion mark addition
			int iNewTriaIndexForDeletionMark = GetNewTriaHavingCommonEdge(vElements, pCommonNode1, pCommonNode2);

			if (-1 != iNewTriaIndexForDeletionMark)
				vElements.at(iNewTriaIndexForDeletionMark).SetToBeDeleted(true);
			else
			{
				vElements.at(iElementIndex).SetToBeDeleted(false);
				continue;
			}

			Tria tTriaNew1(Edge2(pCurrNode, pUncommonNode), Edge2(pUncommonNode, pCommonNode1), Edge2(pCommonNode1, pCurrNode));
			vElements.push_back(tTriaNew1);

			Tria tTriaNew2(Edge2(pUncommonNode, pCurrNode), Edge2(pCurrNode, pCommonNode2), Edge2(pCommonNode2, pUncommonNode));
			vElements.push_back(tTriaNew2);
		}
	}

	DeleteSuperTriaConnectedElements(vElements, vSuperTriaNodes);

	for (int iElementIndex = vElements.size() - 1; iElementIndex >= 0; iElementIndex--)
	{
		Tria pElem = vElements.at(iElementIndex);

		if (true == pElem.GetToBeDeleted())
			vElements.erase(vElements.begin() + iElementIndex);
	}

	dPolygonArea = CalculateAreaUsingTriaElems(vElements);

	return true;
}

std::vector<Line> GetEdgeLines(Polygon& pPolygon)
{
	std::vector<Line> vEdgeList;

	Point pCurrPoint = pPolygon.GetHeadPoint();
	for (int iPointIndex = 0; iPointIndex < pPolygon.GetEdgesCount(); iPointIndex++)
	{
		Point pNextPoint = *pCurrPoint.pNext;
		Line pEdgeLine(pCurrPoint, pNextPoint);
		pEdgeLine.SetSegment(true);
		vEdgeList.push_back(pEdgeLine);

		pCurrPoint = *pCurrPoint.pNext;
	}

	return vEdgeList;
}

void AddPointInPolygonLinkedList(Line& pLine, Point* pIntersectionPoint)
{
	/*Point* pStartPoint = &pLine.GetStart();
	pStartPoint->pNext = pIntersectionPoint;
	pIntersectionPoint->pPrev = pStartPoint;

	Point* pEndPoint = &pLine.GetEnd();
	pEndPoint->pPrev = pIntersectionPoint;
	pIntersectionPoint->pNext = pEndPoint;*/
}

std::vector<Point> GetIntersections(std::vector<Polygon>& vPolygons)
{
	// Polygons should be convex and not concave
	//

	std::vector<Point> vIntersections;

	if (vPolygons.size() < 2 || vPolygons.size() > 3)
		return vIntersections;

	Polygon pLeft = vPolygons.at(0);
	Polygon pRight = vPolygons.at(1);

	std::vector<Line> vEdgesLeftPolygon = GetEdgeLines(pLeft);
	std::vector<Line> vEdgesRightPolygon = GetEdgeLines(pRight);

	for (int iLeftLineIndex = 0; iLeftLineIndex < vEdgesLeftPolygon.size(); iLeftLineIndex++)
	{
		Line pLeftPolygonLine = vEdgesLeftPolygon.at(iLeftLineIndex);

		for (int iRightLineIndex = 0; iRightLineIndex < vEdgesRightPolygon.size(); iRightLineIndex++)
		{
			Line pRightPolygonLine = vEdgesRightPolygon.at(iRightLineIndex);

			Point pIntersection;
			bool bIsIntersection = GetLinesIntersectionPoint(pLeftPolygonLine, pRightPolygonLine, false, pIntersection);
			if (true == bIsIntersection)
			{
				vIntersections.push_back(pIntersection);

				// Left polygon Point addition
				Point pIntersectionCopyLeftTemp = pIntersection.Copy();
				Point* pIntersectionCopyLeft = &pIntersectionCopyLeftTemp;

				Point pStartPointLeft = pLeftPolygonLine.GetStart();
				for (int iLeftPointIndex = 0; iLeftPointIndex < pLeft.GetEdgesCount(); iLeftPointIndex++)
				{
					Point pCurrPointLeftTemp = pLeft.GetHeadPoint();
					Point* pCurrPointLeft = &pCurrPointLeftTemp;
					if (pStartPointLeft == *pCurrPointLeft)
					{
						Point* pCurrPointNextTemp = pCurrPointLeft->pNext;

						pCurrPointLeft->pNext = pIntersectionCopyLeft;
						pIntersectionCopyLeft->pPrev = pCurrPointNextTemp;

						pCurrPointNextTemp->pPrev = pIntersectionCopyLeft;
						pIntersectionCopyLeft->pNext = pCurrPointNextTemp;
					}
				}

				// Right polygon Point addition
				Point pIntersectionCopyRightTemp = pIntersection.Copy();
				Point* pIntersectionCopyRight = &pIntersectionCopyRightTemp;

				Point pStartPointRight = pRightPolygonLine.GetStart();
				for (int iLeftPointIndex = 0; iLeftPointIndex < pLeft.GetEdgesCount(); iLeftPointIndex++)
				{
					Point pCurrPointRIghtTemp = pLeft.GetHeadPoint();
					Point* pCurrPointRight = &pCurrPointRIghtTemp;
					if (pStartPointRight == *pCurrPointRight)
					{
						Point* pCurrPointNextTemp = pCurrPointRight->pNext;

						pCurrPointRight->pNext = pIntersectionCopyLeft;
						pIntersectionCopyLeft->pPrev = pCurrPointRight;

						pCurrPointNextTemp->pPrev = pIntersectionCopyLeft;
						pIntersectionCopyLeft->pNext = pCurrPointNextTemp;
					}
				}
			}
		}
	}

	// Create 4 chains
	std::vector<std::vector<Point>> vvFourChains;
	std::vector<Point> vChain1;
	std::vector<Point> vChain2;
	std::vector<Point> vChain3;
	std::vector<Point> vChain4;
	vvFourChains.push_back(vChain1);
	vvFourChains.push_back(vChain2);
	vvFourChains.push_back(vChain3);
	vvFourChains.push_back(vChain4);

	int iChainIndex = 0;
	bool bStartStopChain1Creation = false;
	bool bStartStopChain2Creation = false;
	for (int iPolygonIndex = 0; iPolygonIndex < vPolygons.size(); iPolygonIndex++)
	{
		Point pCurrPointChain1Temp = vPolygons.at(iPolygonIndex).GetHeadPoint();
		Point* pCurrPointChain1 = &pCurrPointChain1Temp;
		Point pCuurPointChain2Temp = vPolygons.at(iPolygonIndex).GetHeadPoint();
		Point* pCurrPointChain2 = &pCuurPointChain2Temp;
		for (int iPointIndex = 0; iPointIndex < vPolygons.at(iPolygonIndex).GetEdgesCount(); iPointIndex++)
		{
			if (*pCurrPointChain1 == vIntersections.at(0))
				bStartStopChain1Creation = true;

			if (*pCurrPointChain2 == vIntersections.at(0))
				bStartStopChain2Creation = true;

			if (true == bStartStopChain1Creation)
				vvFourChains.at(iChainIndex).push_back(*pCurrPointChain1);

			if (true == bStartStopChain2Creation)
				vvFourChains.at(iChainIndex+1).push_back(*pCurrPointChain2);

			if (*pCurrPointChain1 == vIntersections.at(1))
				bStartStopChain1Creation = false;

			if (*pCurrPointChain2 == vIntersections.at(1))
				bStartStopChain2Creation = false;
			
			pCurrPointChain1 = pCurrPointChain1->pNext;
			pCurrPointChain2 = pCurrPointChain2->pPrev;
		}

		iChainIndex = iChainIndex + 2;
	}


}

double GetCommonArea(std::vector<Polygon>& vPolygons)
{
	if (vPolygons.size() < 2 || vPolygons.size() > 3)
		return NULL;

	Polygon pLeft = vPolygons.at(0);
	Polygon pRight = vPolygons.at(1);

	// Get intersections
	std::vector<Point> dIntersectionsCount = GetIntersections(vPolygons);

}

int main()
{
	//Vector Algebra
	/*Vector vVector;
	double vVecSize = vVector.GetSize();
	double vVecValue0 = 0.0;
	vVector.GetValue(0, vVecValue0);
	double vVecValue1 = 0.0;
	vVector.GetValue(1, vVecValue1);

	Vector vVec1(3, 1);
	Vector vVec2(3, 1);
	if (vVec1 == vVec2)
		std::cout << "Vectors are same." << std::endl;
	else
		std::cout << "Vectors are not same." << std::endl;*/

	// IsPointInside a Polygon
	/*Point pIsPointInside(4, 4, 0);

	Point pHeadPoint(8, 0.0, 0.0);
	Polygon pPolygon(pHeadPoint);
	Point pVerifyHeadPoint = pPolygon.GetHeadPoint();
	pPolygon.AddPoint(Point(12, 4, 0));
	pPolygon.AddPoint(Point(12, 8, 0));
	pPolygon.AddPoint(Point(8, 12, 0));
	pPolygon.AddPoint(Point(6, 6, 0));
	pPolygon.AddPoint(Point(4, 12, 0));
	pPolygon.AddPoint(Point(0, 0, 0), true);

	bool bOnEdge = false;
	bool bOnVertex = false;
	bool bOnLineOfEdge = false;
	Polygon::MethodType eMethodType = Polygon::MethodType::VECTOR_ALGEBRA;//Polygon::MethodType::RAY_CASTING
	bool bIsPointInside = pPolygon.IsPointInside(pIsPointInside, eMethodType, bOnEdge, bOnVertex, bOnLineOfEdge);*/

	// IsPointInside a Circle
	/*Point pIsPointInside(1, 0, 0);
	
	Circle cCircle(Point(0.0, 0.0, 0.0), 5);

	bool bOnCircum = false;
	bool bIsPointInside = cCircle.IsPointInside(pIsPointInside, bOnCircum);*/

	// Finding if a Polygon has a concave angle (Then add the condition in the isPointInside Methods of Vector Algebra in Polygon and for Area)
	/*Point pHeadPoint(0.0, 0.0, 0.0);
	Polygon pPolygon(pHeadPoint);
	Point pVerifyHeadPoint = pPolygon.GetHeadPoint();
	pPolygon.AddPoint(Point(10, 0, 0));
	pPolygon.AddPoint(Point(11, 10, 0));
	pPolygon.AddPoint(Point(5, 4, 0));
	pPolygon.AddPoint(Point(0, 10, 0), true);

	bool bIsPolygonConvex = pPolygon.IsConvexType();*/

	// Generation of 2D convex hull by Gift Wrapping algorithm
	/*std::vector<Point> vPointsGiftWrapAlgo = {{8, 0, 0},
									{12, 4, 0},
									{13, 8, 0},
									{8, 12, 0},
									{4, 13, 0},
									{0, 8, 0},
									{2, 2, 0},
									{8, 4, 0},
									{6, 6, 0},
									{4, 8, 0} };
	
	std::vector<Point> vConvexHullPointsUsingGiftWrapping;
	CreateConvexHullUsingGiftWrap(vPointsGiftWrapAlgo, vConvexHullPointsUsingGiftWrapping);

	// Generation of 2D convex hull by Graham Scan algorithm
	std::vector<Point> vPointsGrahamScanAlgo = {	{8, 0, 0},
									{12, 4, 0},
									{13, 8, 0},
									{8, 12, 0},
									{4, 13, 0},
									{0, 8, 0},
									{2, 2, 0},
									{8, 4, 0},
									{6, 6, 0},
									{4, 8, 0} };

	std::vector<Point> vConvexHullPointsUsingGrahamScan;
	CreateConvexHullUsingGrahamScan(vPointsGrahamScanAlgo, vConvexHullPointsUsingGrahamScan);*/
	
	// IsPointInside a 3D Polygon/Polyhedron
	/*|\
	  | \
	  | I\
______|_. \______> Is I inside the triangular face
	  |____\

	//[ -DX, (V1 - V0)X, (V2 - V0)X][  t  ]   [ (O - V0)X ]
	//[ -DY, (V1 - V0)Y, (V2 - V0)Y][  u  ] = [ (O - V0)Y ]
	//[ -DZ, (V1 - V0)Z, (V2 - V0)Z][  v  ]   [ (O - V0)Z ]

	// Conditions to be satidfied. (u, v, w are barycentric coordinates and u + v + w = 1)
	// u > 0 and u < 1
	// v > 0 and (u + v) < 1
	// t > 0 becasue we need intersection in the direction of ray starting from ray's Origin O

	// Moeller-Trumbore Algorithm and Cramer's rule

	std::vector<TriangularFace> polyhedronFaces;
	polyhedronFaces.push_back(TriangularFace(Point(0, 0, 0), Point(3, 0, 0), Point(3, 5, 0)));
	polyhedronFaces.push_back(TriangularFace(Point(0, 0, 0), Point(3, 0, 0), Point(2, 3, 5)));
	polyhedronFaces.push_back(TriangularFace(Point(2, 3, 5), Point(3, 0, 0), Point(3, 5, 0)));
	polyhedronFaces.push_back(TriangularFace(Point(0, 0, 0), Point(2, 3, 5), Point(3, 5, 0)));

	Point pO(2, 2.5, 3); // Inside
	//Point pO(5, 2.5, 3); // Outside

	Point pD(1, 0, 0);

	bool bIsPintInsidePolyhedron = IsPointInsidePolyhedron(pO, pD, polyhedronFaces);*/

	// Generation of 3D convex hull
	/*std::vector<Point> vPointsClarksonShorAlgo = {  {0, 0, 0},
													{1, 0, 1},
													{0.5, 4, 3},
													{2, 2, 5},
													{3, 3.5, 2},
													{-1, 2, 3},
													{-2, -1, 3},
													{3, -3, 3},
													{-2, 3, 0},
													{-3, 1, 4},
													{3, 1, 3} };*/

	/*std::vector<Point> vPointsClarksonShorAlgo = {  {0, 3, 0},
													{0, 2, 1},
													{4, 2, 1},
													{5, 2, 0},
													{6, 0, 0},
													{5, 0, 1},
													{0, 0, 1.5},
													{0, 0, 0},
													{3, 1, 0.5},
													{3.5, 2, 1} };

	std::vector<TriangularFace> vConvexHullFaces;
	Create3DConvexHullUsingQuickHull(vPointsClarksonShorAlgo, vConvexHullFaces);*/
	
	// Area of an irregular polygon
	// https://www.pre-scient.com/knowledge-center/product-development-by-reverse-engineering/meshing-algorithms/
	/* Shoelace formula :
	// Take one point. 
	// Create a triangle with each set of 2 consecutive points.
	// Take cross prodict and divide by 2 each time.
	// Add all to find the area of polygon
	// It is not usitable for polygon having concave angle.*/
	// Tria Meshing using Delaunay Triangulation Algorithm
	/*std::vector<Node> vPointsDelaunayAlgo = {{8, 0, 0},
									{12, 4, 0},
									{13, 8, 0},
									{8, 12, 0},
									{4, 13, 0},
									{0, 8, 0},
									{2, 2, 0},
									{8, 4, 0},
									{6, 6, 0},
									{4, 8, 0} };*/

	/*std::vector<Node> vPointsDelaunayAlgo = {	{3, 1, 0},
												{5, 2, 0},
												{1, 3, 0},
												{3, 4, 0},
												{7, 3, 0},
												{7, -3, 0} };*/

	/*std::vector<Node> vPointsDelaunayAlgo = {	{3, 1, 0},
												{3, 2, 0},
												{4, 1, 0} };
	
	std::vector<Tria> vElements;
	double dPolygonArea = 0;
	bool bMeshUsingDelaunayAlgo = MeshUsingDelaunayAlgo(vPointsDelaunayAlgo, vElements, dPolygonArea);*/


	// Intersection line-line, plane-plane, line-plane
	// Intersection of lines and segments in 3D space
	/*Line line1(Point(2, 0, 0), Point(-3.01, 2.01, 4.24));
	line1.SetSegment(true);
	//Line line2(Point(-2, 0, 0), Point(1.74, 2.16, 4.56)); 2.16 intersects in Y and 5 does not.
	//line2.SetSegment(false);

	// For checking Parallelism case
	Line line2(Point(4, 0, 0), Point(-1.01, 2.01, 4.24));
	line1.SetSegment(true);

	double dInitialiseInf = std::numeric_limits<double>::infinity();
	Point pIntersectionPoint(dInitialiseInf, dInitialiseInf, dInitialiseInf);
	bool bIntersectionPoint = GetLinesIntersectionPoint(line1, line2, true, pIntersectionPoint);*/

	// intersection of planes in 3D space
	/*Plane plane1(1, 2, 3, 10);
	//Plane plane2(2, 3, 1, 15); // Non-parallel plane
	Plane plane2(2, 4, 6, 20);// Parallel plane

	Line lIntersectionLine(Point(0, 0, 0), Point(0, 0, 0));

	bool bIntersectionLine = GetPlanesIntersetionLine(plane1, plane2, lIntersectionLine);*/

	// Intersection of 2D shapes and finding common area
	Point pLeftHeadPoint(8, 2, 0);
	Polygon pLeft;
	pLeft.AddPoint(pLeftHeadPoint);
	pLeft.AddPoint(Point(4, 6, 0));
	pLeft.AddPoint(Point(8, 10, 0));
	pLeft.AddPoint(Point(12, 6, 0), true);

	Point pRightHeadPoint(16, -2, 0);
	Polygon pRight;
	pRight.AddPoint(pRightHeadPoint);
	pRight.AddPoint(Point(8, 6, 0));
	pRight.AddPoint(Point(16, 14, 0), true);

	std::vector<Polygon> vPolygons;
	vPolygons.push_back(pLeft);
	vPolygons.push_back(pRight);

	double dIntersectionArea = GetCommonArea(vPolygons);

	// Intersection of 3D shapes and finding common volume


	// OpenGL
	/*if (!glfwInit())
		return 0;

	GLFWwindow* window = glfwCreateWindow(640, 480, "OpenGL Graphics Window", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return 0;
	}

	glfwMakeContextCurrent(window);

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glOrtho(-width/2, width / 2, -height/2, height / 2, -1, 1);

	while (!glfwWindowShouldClose(window))
	{
		drawLine();
		drawPoint();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();*/

	return 1;
}

