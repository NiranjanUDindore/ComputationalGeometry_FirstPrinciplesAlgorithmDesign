#include "Vector.h"
#include "Point.h"

void Vector::Create3DVector(Point& pStartPoint, Point& pEndPoint)
{
	vData.resize(3);
	vData[0] = pEndPoint.GetX() - pStartPoint.GetX();
	vData[1] = pEndPoint.GetY() - pStartPoint.GetY();
	vData[2] = pEndPoint.GetZ() - pStartPoint.GetZ();
}

void Vector::Create3DVector(Node& pStartPoint, Node& pEndPoint)
{
	vData.resize(3);
	vData[0] = pEndPoint.GetX() - pStartPoint.GetX();
	vData[1] = pEndPoint.GetY() - pStartPoint.GetY();
	vData[2] = pEndPoint.GetZ() - pStartPoint.GetZ();
}

void Vector::SetValue(int iIndex, double dValue)
{
	vData.at(iIndex) = dValue;
}

bool Vector::GetValue(int iIndex, double& dValue)
{
	if (iIndex < vData.size())
	{
		dValue = vData.at(iIndex);
		return true;
	}
	
	dValue = -1.0;

	return false;
}

double Vector::GetSize()
{
	return vData.size();
}

double Vector::GetMagnitude()
{
	double dMag = 0;
	for (int iIndex = 0; iIndex < vData.size(); iIndex++)
	{
		dMag += vData.at(iIndex) * vData.at(iIndex);
	}

	return sqrt(dMag);
}

void Vector::Normalize()
{
	double dMag = GetMagnitude();
	for (int iIndex = 0; iIndex < vData.size(); iIndex++)
	{
		vData.at(iIndex) = vData.at(iIndex) / dMag;
	}
}

double Vector::GetX()
{
	return vData.at(0);
}

double Vector::GetY()
{
	return vData.at(1);
}

double Vector::GetZ()
{
	return vData.at(2);
}

bool Vector::operator==(const Vector& other) const
{
	if (vData.size() != other.vData.size())
		return false;

	for (int iIndex = 0; iIndex < vData.size(); iIndex++)
	{
		if (vData.at(iIndex) != other.vData.at(iIndex))
			return false;
	}

	return true;
}

double Vector::CrossProduct2D(Vector& vOtherVector)
{
	double dCrossProd = (this->GetX() * vOtherVector.GetY()) - (this->GetY() * vOtherVector.GetX());

	return dCrossProd;
}

Vector Vector::CrossProduct3D(Vector& vOtherVector)
{
	std::vector<double> vCrossedVectorCoefficients;

	double dCoeffi = (this->GetY() * vOtherVector.GetZ()) - (this->GetZ() * vOtherVector.GetY());
	double dCoeffj = (this->GetZ() * vOtherVector.GetX()) - (this->GetX() * vOtherVector.GetZ());
	double dCoeffk = (this->GetX() * vOtherVector.GetY()) - (this->GetY() * vOtherVector.GetX());

	vCrossedVectorCoefficients.push_back(dCoeffi);
	vCrossedVectorCoefficients.push_back(dCoeffj);
	vCrossedVectorCoefficients.push_back(dCoeffk);

	return Vector(vCrossedVectorCoefficients);
}

double Vector::DotProduct3D(Vector& vOtherVector)
{
	double dDotProduct = this->GetX() * vOtherVector.GetX()
		+ this->GetY() * vOtherVector.GetY()
		+ this->GetZ() * vOtherVector.GetZ();

	return dDotProduct;
}

void Vector::Reverse()
{
	this->SetValue(0, -this->GetX());
	this->SetValue(1, -this->GetY());
	this->SetValue(2, -this->GetZ());
}
