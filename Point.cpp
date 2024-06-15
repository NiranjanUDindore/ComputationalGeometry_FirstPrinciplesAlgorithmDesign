#include "Point.h"

bool Point::operator!=(const Point& pOtherPoint) const
{
    if ((this->dX == pOtherPoint.dX) && (this->dY == pOtherPoint.dY) && (this->dZ == pOtherPoint.dZ))
        return false;

    return true;
}

bool Point::operator==(const Point& pOtherPoint) const
{
    if ((this->dX == pOtherPoint.dX) && (this->dY == pOtherPoint.dY) && (this->dZ == pOtherPoint.dZ))
        return true;

    return false;
}

Point Point::Copy()
{
    Point pCopyPoint(this->GetX(), this->GetY(), this->GetZ());

    return pCopyPoint;
}
