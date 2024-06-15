#pragma once
#include <iostream>
#include <vector>

class Matrix
{
private:
	std::vector<std::vector<double>> vData;
	int iRows;
	int iCols;

public:
	Matrix();

	Matrix(int _iRows, int _iCols) : iRows(_iRows), iCols(_iCols)
	{
		vData.resize(iRows, std::vector<double>(iCols, 0.0));
	}

	double& operator()(int iRows, int iCols)
	{
		return vData[iRows][iCols];
	}

	double operator()(int iRows, int iCols) const
	{
		return vData[iRows][iCols];
	}

	double GetDeterminantIf3X3()
	{
		if (3 != iRows || 3 != iCols)
			return NULL;

		double dDet = vData[0][0] * (vData[1][1] * vData[2][2] - vData[1][2] * vData[2][1])
			- vData[0][1] * (vData[1][0] * vData[2][2] - vData[1][2] * vData[2][0])
			+ vData[0][2] * (vData[1][0] * vData[2][1] - vData[1][1] * vData[2][0]);
		
		return dDet;
	}

};

