#pragma once
//#include "Tria.h"
class Node
{
private:
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;

	//std::vector<Tria> vConnectedElems;

public:
	Node() /*: pNext(nullptr), pPrev(nullptr)*/ {};//

	Node(double _dX, double _dY, double _dZ) : dX(_dX), dY(_dY), dZ(_dZ) {};

	//Node* pNext = nullptr;
	//Node* pPrev = nullptr;//

	void SetX(double dSetX) { dX = dSetX; }
	double GetX() { return dX; }

	void SetY(double dSetY) { dY = dSetY; }
	double GetY() { return dY; }

	void SetZ(double dSetZ) { dZ = dSetZ; }
	double GetZ() { return dZ; }

	//void AddElementConnection(Tria& pElement) { vConnectedElems.push_back(pElement); }

	bool operator==(const Node& pOtherNode) const;
	bool operator==(const Node* pOtherNode) const;

	Node Copy();
};

