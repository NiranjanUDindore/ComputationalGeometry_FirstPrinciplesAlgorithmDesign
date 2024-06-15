#pragma once
#include "Node.h"
class Edge2
{
private:
	Node nNode1;
	Node nNode2;

public:
	Edge2() : nNode1(0.0, 0.0, 0.0), nNode2(0.0, 0.0, 0.0) {};

	Edge2(Node _nNode1, Node _nNode2) : nNode1(_nNode1), nNode2(_nNode2) {};

	void SetX1(double dSetX1) { nNode1.SetX(dSetX1); }
	double GetX1() { return nNode1.GetX(); }

	void SetY1(double dSetY1) { nNode1.SetX(dSetY1); }
	double GetY1() { return nNode1.GetY(); }

	void SetZ1(double dSetZ1) { nNode1.SetX(dSetZ1); }
	double GetZ1() { return nNode1.GetZ(); }

	void SetX2(double dSetX2) { nNode2.SetX(dSetX2); }
	double GetX2() { return nNode2.GetX(); }

	void SetY2(double dSetY2) { nNode2.SetX(dSetY2); }
	double GetY2() { return nNode2.GetY(); }

	void SetZ2(double dSetZ2) { nNode2.SetX(dSetZ2); }
	double GetZ2() { return nNode2.GetZ(); }
};

