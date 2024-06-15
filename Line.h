#pragma once
#include <iostream>
#include "Point.h"
class Line
{
private:
	Point pStart;
	Point pEnd;

	bool bIsSegment = false;

public:
	Line();

	Line(Point _pStart, Point _pEnd) : pStart(_pStart), pEnd(_pEnd) {};

	Point GetStart() { return pStart; }
	Point GetEnd() { return pEnd; }
	void SetStart(Point pSetStart) { pStart = pSetStart; }
	void SetEnd(Point pSetEnd) { pEnd = pSetEnd; }
	void SetSegment(bool bSetSeg) { bIsSegment = bSetSeg; }
	bool IsSegment() { return bIsSegment; }
};

