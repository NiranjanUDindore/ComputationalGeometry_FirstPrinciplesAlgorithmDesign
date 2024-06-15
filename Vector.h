#pragma once

#include <vector>
#include "Point.h"
#include "Node.h"

class Vector
{
private:
	std::vector<double> vData = { 0 };

public:
	Vector() {};

	Vector(int _iDimensions, double _dValue) : vData(_iDimensions, _dValue) {};

	Vector(std::vector<double> _vData) : vData(_vData) {}; 

	void Create3DVector(Point& pStartPoint, Point& pEndPoint);
	void Create3DVector(Node& pStartPoint, Node& pEndPoint);

	void SetValue(int iIndex, double dValue);
	bool GetValue(int iIndex, double& dValue);
	double GetSize();
	double GetMagnitude();
	void Normalize();

	double GetX();
	double GetY();
	double GetZ();

	bool operator==(const Vector& other) const;

	double CrossProduct2D(Vector& vOtherVector);
	Vector CrossProduct3D(Vector& vOtherVector);
	double DotProduct3D(Vector& vOtherVector);
	void Reverse();
};

