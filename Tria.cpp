#include "Tria.h"

bool Tria::GetCommonNodes(Tria& pOtherElem, Node& node1, Node& node2)
{
	bool bIsNode1Assigned = false;
	bool bIsNode2Assigned = false;

	std::vector<Node> vNodeVecThis;
	vNodeVecThis.push_back(p1);
	vNodeVecThis.push_back(p2);
	vNodeVecThis.push_back(p3);

	Node P1 = pOtherElem.P1();
	Node P2 = pOtherElem.P2();
	Node P3 = pOtherElem.P3();
	std::vector<Node> vNodeVecOther;
	vNodeVecOther.push_back(P1);
	vNodeVecOther.push_back(P2);
	vNodeVecOther.push_back(P3);

	for (int iThisCount = 0; iThisCount < vNodeVecThis.size(); iThisCount++)
	{
		Node nThisP = vNodeVecThis.at(iThisCount);

		for (int iOtherCount = 0; iOtherCount < vNodeVecOther.size(); iOtherCount++)
		{
			Node nOtherP = vNodeVecOther.at(iOtherCount);

			if (nThisP == nOtherP)
			{
				if (false == bIsNode1Assigned)
				{
					node1 = nThisP;
					bIsNode1Assigned = true;
				}
				else
				{
					node2 = nThisP;
					bIsNode2Assigned = true;
					break;
				}
			}
		}
		if (true == bIsNode2Assigned)
			break;
	}

	if (false == bIsNode1Assigned || false == bIsNode2Assigned)
		return false;
	else
		return true;
}
