#pragma once
class Plane
{
private:
	double a = 0;
	double b = 0;
	double c = 0;
	double d = 0;

public:
	Plane() {};

	Plane(double _a, double _b, double _c, double _d) : a(_a), b(_b), c(_c), d(_d) {};

	double GetA() { return a; }
	double GetB() { return b; }
	double GetC() { return c; }
	double GetD() { return d; }

	void SetA(double dSetA) { a = dSetA; }
	void SetB(double dSetB) { a = dSetB; }
	void SetC(double dSetC) { a = dSetC; }
	void SetD(double dSetD) { a = dSetD; }
};

