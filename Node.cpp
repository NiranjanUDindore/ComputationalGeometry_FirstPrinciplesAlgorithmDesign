#include "Node.h"

bool Node::operator==(const Node& pOtherPoint) const
{
    if ((this->dX == pOtherPoint.dX) && (this->dY == pOtherPoint.dY) && (this->dZ == pOtherPoint.dZ))
        return true;

    return false;
}

bool Node::operator==(const Node* pOtherPoint) const
{
    if ((this->dX == pOtherPoint->dX) && (this->dY == pOtherPoint->dY) && (this->dZ == pOtherPoint->dZ))
        return true;

    return false;
}

Node Node::Copy()
{
    Node nCopyNode(this->GetX(), this->GetY(), this->GetZ());

    return nCopyNode;
}
