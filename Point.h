#pragma once
class Point
{
private:
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;

public:
	Point() : pNext(nullptr), pPrev(nullptr) {};//

	Point(double _dX, double _dY, double _dZ) : dX(_dX), dY(_dY), dZ(_dZ) {};

	Point* pNext = nullptr;
	Point* pPrev = nullptr;//

	void SetX(double dSetX) { dX = dSetX; }
	double GetX() { return dX; }

	void SetY(double dSetY) { dY = dSetY; }
	double GetY() { return dY; }

	void SetZ(double dSetZ) { dZ = dSetZ; }
	double GetZ() { return dZ; }

	bool operator!=(const Point& pOtherPoint) const;
	bool operator==(const Point& pOtherPoint) const;

	Point Copy();
};

