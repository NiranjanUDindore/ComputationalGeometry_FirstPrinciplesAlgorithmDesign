#pragma once
#include <iostream>
#include <vector>
#include "TriangularFace.h"

class Tetrahedron
{
private:
	std::vector<TriangularFace> vTriaFaces;

public:
	Tetrahedron() {};

	Tetrahedron(TriangularFace triaFace1, TriangularFace triaFace2, TriangularFace triaFace3, TriangularFace triaFace4)
	{
		vTriaFaces.push_back(triaFace1);
		vTriaFaces.push_back(triaFace2);
		vTriaFaces.push_back(triaFace3);
		vTriaFaces.push_back(triaFace4);
	}

	std::vector<TriangularFace> GetFaces() { return vTriaFaces; }
};

