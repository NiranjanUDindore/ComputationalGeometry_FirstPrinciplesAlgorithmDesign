#pragma once
#include <iostream>
#include <vector>
#include "Node.h"
#include "Point.h"
#include "Edge.h"
#include "Edge2.h"
class Tria
{
private:
	Node p1;
	Node p2;
	Node p3;

	Edge2 e1;
	Edge2 e2;
	Edge2 e3;

	bool bToBeDeleted = false;

public:
	Tria()/* : p1(0.0, 0.0, 0.0), p2(0.0, 0.0, 0.0), p3(0.0, 0.0, 0.0) */{};

	Tria(Edge2 _e1, Edge2 _e2, Edge2 _e3) : e1(_e1), e2(_e2), e3(_e3)
	{
		Node _p1(e1.GetX1(), e1.GetY1(), e1.GetZ1());
		p1 = _p1;

		Node _p2(e1.GetX2(), e1.GetY2(), e1.GetZ2());
		p2 = _p2;

		Node _p3(e2.GetX2(), e2.GetY2(), e2.GetZ2());
		p3 = _p3;
	};

	Node P1() { return p1; }
	Node P2() { return p2; }
	Node P3() { return p3; }

	Edge2 GetEdgeP1P2() { return Edge2(p1, p2); }
	Edge2 GetEdgeP2P3() { return Edge2(p2, p3); }

	void SetToBeDeleted(bool bSetToBeDeleted) { bToBeDeleted = bSetToBeDeleted; }
	bool GetToBeDeleted() { return bToBeDeleted; }

	bool GetCommonNodes(Tria& pOtherElem, Node& node1, Node& node2);
};

