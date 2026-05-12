//////////////////////////////////////////////////////////////////////////
/// Copyright(c) 2010-2015 CJBUAA
/// All right reserved.
///
/// @file: Math.h
///
/// @brief: Definition of some useful math functions under namespace sc_math.
///
/// @author: ChengJian
/// @date: March 30, 2012
/// @version: 1.0
/// @date: Noverber 13, 2013
/// @version: 1.1
//////////////////////////////////////////////////////////////////////////
#ifndef __MATH_H_
#define __MATH_H_

#include "Array.h"
#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

#ifndef NDEBUG
#include <cassert>
#include <cstdlib>
#endif

namespace sc_math
{
	static const double E = 2.7182818284590452353602874713527; /**<自然对数的底{exp(1)}*/
	static const double PI = 4.0 * atan(1.0);				   /**< 圆周率{2*arcsin(1)}*/
	static const double SQRT_2 = sqrt(2.0);					   /**< sqrt(2)*/
	static const double SQRT_3 = sqrt(3.0);					   /**< sqrt(3)*/

	/** \brief get machine zero */
	inline double getMachineZero(void)
	{
		return (std::numeric_limits<double>::epsilon());
	}

	/** \brief sign function */
	inline int Sign(int var)
	{
		return (var >= 0 ? 1 : -1);
	}
	inline int Sign(double var)
	{
		return (var >= 0 ? 1 : -1);
	}
	/** \brief min & max function */
	inline double Min(double a, double b)
	{
		return (a < b ? a : b);
	}
	inline double Min(double a, double b, double c)
	{
		return (a < Min(b, c) ? a : Min(b, c));
	}
	inline double Max(double a, double b)
	{
		return (a > b ? a : b);
	}
	inline double Max(double a, double b, double c)
	{
		return (a > Max(b, c) ? a : Max(b, c));
	}
	inline void Max(unsigned int n, Array1D<double> &v, double &max, unsigned int &index)
	{
		max = v(0);
		index = 0;

		for (unsigned int i = 1; i != n; ++i)
		{
			if (v(i) > max)
			{
				max = v(i);
				index = i;
			}
		}
	}
	inline void Min(unsigned int n, Array1D<double> &v, double &min, unsigned int &index)
	{
		min = v(0);
		index = 0;

		for (unsigned int i = 1; i != n; i++)
		{
			if (v(i) < min)
			{
				min = v(i);
				index = i;
			}
		}
	}
	/** \brief minmod function */
	inline double Minmod(double a, double b)
	{
		double temp = fabs(a);

		if ((a < 0 && b < 0) || (a > 0 && b > 0))
		{
			return Sign(a) * Min(temp, fabs(b));
		}
		else
		{
			return 0;
		}
	}
	inline double Minmod(double a, double b, double c)
	{
		if ((a < 0 && b < 0 && c < 0) || (a > 0 && b > 0 && c > 0))
		{
			double temp = fabs(a);
			temp = Min(temp, fabs(b));
			temp = Min(temp, fabs(c));

			return Sign(a) * temp;
		}
		else
		{
			return 0;
		}
	}

	inline double Minmod(double a, double b, double c, double M, double h)
	{
		if (fabs(a) < M * h * h)
		{
			return a;
		}
		if (fabs(a) < 1.0e-6 && fabs(b) < 1.0e-6 && fabs(c) < 1.0e-6)
		{
			return a;
		}
		else
		{
			if ((a < 0 && b < 0 && c < 0) || (a > 0 && b > 0 && c > 0))
			{
				double temp = fabs(a);
				temp = Min(temp, fabs(b));
				temp = Min(temp, fabs(c));

				return Sign(a) * temp;
			}
			else
			{
				return 0;
			}
		}
	}

	inline double Minmod(double a, double b, double c, double d)
	{
		double temp = fabs(a);
		if ((a < 0 && b < 0 && c < 0 && d < 0) || (a > 0 && b > 0 && c > 0 && d > 0))
		{
			temp = Min(temp, fabs(b));
			temp = Min(temp, fabs(c));
			temp = Min(temp, fabs(d));
			return Sign(a) * temp;
		}
		else
		{
			return 0;
		}
	}
	/** \brief convert degree to arc and vice versa */
	inline double Deg2arc(double x)
	{
		return x * PI / 180.0;
	}
	inline double Arc2deg(double x)
	{
		return x / PI * 180.0;
	}
	/** \brief square & cubic*/
	inline double Square(double x)
	{
		return x * x;
	}
	inline double Cubic(double x)
	{
		return x * x * x;
	}
	/** \brief log function*/
	inline double Log10(double x)
	{
		return (log(x) / log(10.0));
	}
	inline double Log2(double x)
	{
		return (log(x) / log(2.0));
	}
	/**
	 * Get vector norm-2
	 */
	inline double norm_2(int dim, const Array1D<double> &v)
	{
		double s(0);

		for (int i = 0; i != dim; ++i)
		{
			s += v[i] * v[i];
		}

		return (sqrt(s));
	}
	inline double norm_2(int dim, double *v)
	{
		double s(0);

		for (int i = 0; i != dim; ++i)
		{
			s += v[i] * v[i];
		}

		return (sqrt(s));
	}

	/**
	 * Calculate detetminant of a 3*3 matrix
	 * \param[in] A : matrix of 3*3
	 * \return: det(A)
	 */
	inline double MatrixDet3(const Array2D<double> &A)
	{
		double s(0);

		s = A[0][0] * A[1][1] * A[2][2] + A[1][0] * A[2][1] * A[0][2] + A[2][0] * A[0][1] * A[1][2] - A[0][0] * A[2][1] * A[1][2] - A[1][0] * A[0][1] * A[2][2] - A[2][0] * A[1][1] * A[0][2];

		return s;
	}

	/**
	 * Calculate detetminant of a 4*4 matrix
	 * \param[in] A : matrix of 4*4
	 * \return: det(A)
	 */
	inline double MatrixDet4(const Array2D<double> &A)
	{
		double s(0);

		s = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * (A[2][2] * A[3][3] - A[3][2] * A[2][3]) - (A[0][0] * A[1][2] - A[1][0] * A[0][2]) * (A[2][1] * A[3][3] - A[3][1] * A[2][3]) + (A[0][0] * A[1][3] - A[1][0] * A[0][3]) * (A[2][1] * A[3][2] - A[3][1] * A[2][2]) + (A[0][1] * A[1][2] - A[1][1] * A[0][2]) * (A[2][0] * A[3][3] - A[3][0] * A[2][3]) - (A[0][1] * A[1][3] - A[1][1] * A[0][3]) * (A[2][0] * A[3][2] - A[3][0] * A[2][2]) + (A[0][2] * A[1][3] - A[1][2] * A[0][3]) * (A[2][0] * A[3][1] - A[3][0] * A[2][1]);

		return s;
	}

	/**
	 * \brief  Computes the QR decomposition using Givens Rotations
	 *         A = mxn matrix
	 *         Q = mxm matrix orthogonal matrix
	 *         R = mxn matrix upper triangular matrix
	 *
	 * \param[in]: m,n the number of row and column of A
	 * \param[in]: A
	 * \param[out]: Qt the orthogonal matrix
	 * \param[out]: R the upper triangular matrix
	 */
	inline void QR_decomposition(int m, int n, Array2D<double> &A, Array2D<double> &Qt, Array2D<double> &R)
	{
		double eps(1.0e-6);
		double theta(0), c(0), s(0);
		double temp(0);

		// copy A to R
		for (int i = 0; i != m; ++i)
		{
			for (int j = 0; j != n; ++j)
			{
				R[i][j] = A[i][j];
			}
		}

		// set zeros
		Qt.setZero();
		for (int i = 0; i != m; ++i)
		{
			Qt[i][i] = 1.0;
		}

		// QR decomposition
		for (int j = 0; j != n; ++j)
		{
			for (int i = m - 2; i != j - 1; --i)
			{
				theta = sqrt(R[i + 1][j] * R[i + 1][j] + R[i][j] * R[i][j]);

				if (theta > eps)
				{
					// Calculate the rotation matrix
					theta = 1.0 / theta;
					c = -R[i][j] * theta;
					s = R[i + 1][j] * theta;

					// Update R
					for (int k = j; k != n; ++k)
					{
						temp = R[i][k] * c + R[i + 1][k] * (-s);
						R[i + 1][k] = R[i][k] * s + R[i + 1][k] * c;
						R[i][k] = temp;
					}

					// Update Q^t
					for (int k = 0; k != m; ++k)
					{
						temp = Qt[i][k] * c + Qt[i + 1][k] * (-s);
						Qt[i + 1][k] = Qt[i][k] * s + Qt[i + 1][k] * c;
						Qt[i][k] = temp;
					}
				}
			}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature for line
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights :  weight of each gauss point
	 */
	inline void GaussLine(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 0.5;
			gauss_points[0][1] = 0.5;

			gauss_weights[0] = 1.0;

			break;
		}
		case 2: //(P3)
		{
			gauss_points[0][0] = 0.5 * (1.0 + sqrt(1.0 / 3.0));
			gauss_points[0][1] = 0.5 * (1.0 - sqrt(1.0 / 3.0));

			gauss_points[1][0] = 0.5 * (1.0 - sqrt(1.0 / 3.0));
			gauss_points[1][1] = 0.5 * (1.0 + sqrt(1.0 / 3.0));

			gauss_weights[0] = 0.5;
			gauss_weights[1] = 0.5;

			break;
		}
		case 3: //(P5)
		{
			gauss_points[0][0] = 0.5 * (1.0 + sqrt(3.0 / 5.0));
			gauss_points[0][1] = 0.5 * (1.0 - sqrt(3.0 / 5.0));

			gauss_points[1][0] = 0.5;
			gauss_points[1][1] = 0.5;

			gauss_points[2][0] = 0.5 * (1.0 - sqrt(3.0 / 5.0));
			gauss_points[2][1] = 0.5 * (1.0 + sqrt(3.0 / 5.0));

			gauss_weights[0] = 5.0 / 18.0;
			gauss_weights[1] = 8.0 / 18.0;
			gauss_weights[2] = 5.0 / 18.0;

			break;
		}
		case 4: //(P7)
		{
			double gp1 = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
			double gp2 = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));

			double gw1 = (18.0 + sqrt(30.0)) / 36.0;
			double gw2 = (18.0 - sqrt(30.0)) / 36.0;

			gauss_points[0][0] = 0.5 * (1.0 + gp2);
			gauss_points[0][1] = 0.5 * (1.0 - gp2);

			gauss_points[1][0] = 0.5 * (1.0 + gp1);
			gauss_points[1][1] = 0.5 * (1.0 - gp1);

			gauss_points[2][0] = 0.5 * (1.0 - gp1);
			gauss_points[2][1] = 0.5 * (1.0 + gp1);

			gauss_points[3][0] = 0.5 * (1.0 - gp2);
			gauss_points[3][1] = 0.5 * (1.0 + gp2);

			gauss_weights[0] = 0.5 * gw2;
			gauss_weights[1] = 0.5 * gw1;
			gauss_weights[2] = 0.5 * gw1;
			gauss_weights[3] = 0.5 * gw2;

			break;
		}
		case 5: // P(9)
		{
			double gp1 = 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
			double gp2 = 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));

			double gw1 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
			double gw2 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;

			gauss_points[0][0] = 0.5;
			gauss_points[0][1] = 0.5;

			gauss_points[1][0] = 0.5 * (1.0 + gp2);
			gauss_points[1][1] = 0.5 * (1.0 - gp2);

			gauss_points[2][0] = 0.5 * (1.0 + gp1);
			gauss_points[2][1] = 0.5 * (1.0 - gp1);

			gauss_points[3][0] = 0.5 * (1.0 - gp1);
			gauss_points[3][1] = 0.5 * (1.0 + gp1);

			gauss_points[4][0] = 0.5 * (1.0 - gp2);
			gauss_points[4][1] = 0.5 * (1.0 + gp2);

			gauss_weights[0] = 0.5 * 128.0 / 225.0;
			gauss_weights[1] = 0.5 * gw2;
			gauss_weights[2] = 0.5 * gw1;
			gauss_weights[3] = 0.5 * gw1;
			gauss_weights[4] = 0.5 * gw2;

			break;
		}
		default:
		{
			// Can not be here!
			exit(-1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature for quadrangle
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights :  weight of each gauss point
	 *
	 * \note: https://en.wikipedia.org/wiki/Gaussian_quadrature
	 */
	inline void GaussQuad(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 0.0;
			gauss_points[0][1] = 0.0;

			gauss_weights[0] = 4.0;

			break;
		}
		case 4: //(P3)
		{
			double gp = 1.0 / sqrt(3.0);

			gauss_points[0][0] = gp;
			gauss_points[0][1] = gp;
			gauss_points[1][0] = gp;
			gauss_points[1][1] = -gp;
			gauss_points[2][0] = -gp;
			gauss_points[2][1] = gp;
			gauss_points[3][0] = -gp;
			gauss_points[3][1] = -gp;

			gauss_weights[0] = 1.0;
			gauss_weights[1] = 1.0;
			gauss_weights[2] = 1.0;
			gauss_weights[3] = 1.0;

			break;
		}
		case 9: //(P5)
		{
			double gp0 = sqrt(3.0 / 5.0);
			double gp1 = 0.0;
			double gp2 = -sqrt(3.0 / 5.0);

			double w0 = 5.0 / 9.0;
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;

			gauss_points[0][0] = gp0;
			gauss_points[0][1] = gp0;
			gauss_points[1][0] = gp0;
			gauss_points[1][1] = gp1;
			gauss_points[2][0] = gp0;
			gauss_points[2][1] = gp2;

			gauss_points[3][0] = gp1;
			gauss_points[3][1] = gp0;
			gauss_points[4][0] = gp1;
			gauss_points[4][1] = gp1;
			gauss_points[5][0] = gp1;
			gauss_points[5][1] = gp2;

			gauss_points[6][0] = gp2;
			gauss_points[6][1] = gp0;
			gauss_points[7][0] = gp2;
			gauss_points[7][1] = gp1;
			gauss_points[8][0] = gp2;
			gauss_points[8][1] = gp2;

			gauss_weights[0] = w0 * w0;
			gauss_weights[1] = w0 * w1;
			gauss_weights[2] = w0 * w2;

			gauss_weights[3] = w1 * w0;
			gauss_weights[4] = w1 * w1;
			gauss_weights[5] = w1 * w2;

			gauss_weights[6] = w2 * w0;
			gauss_weights[7] = w2 * w1;
			gauss_weights[8] = w2 * w2;

			break;
		}
		case 16: //(P7)
		{
			double gp0 = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
			double gp1 = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
			double gp2 = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
			double gp3 = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));

			double w0 = (18.0 - sqrt(30.0)) / 36.0;
			double w1 = (18.0 + sqrt(30.0)) / 36.0;
			double w2 = (18.0 + sqrt(30.0)) / 36.0;
			double w3 = (18.0 - sqrt(30.0)) / 36.0;

			gauss_points[0][0] = gp0;
			gauss_points[0][1] = gp0;
			gauss_points[1][0] = gp0;
			gauss_points[1][1] = gp1;
			gauss_points[2][0] = gp0;
			gauss_points[2][1] = gp2;
			gauss_points[3][0] = gp0;
			gauss_points[3][1] = gp3;

			gauss_points[4][0] = gp1;
			gauss_points[4][1] = gp0;
			gauss_points[5][0] = gp1;
			gauss_points[5][1] = gp1;
			gauss_points[6][0] = gp1;
			gauss_points[6][1] = gp2;
			gauss_points[7][0] = gp1;
			gauss_points[7][1] = gp3;

			gauss_points[8][0] = gp2;
			gauss_points[8][1] = gp0;
			gauss_points[9][0] = gp2;
			gauss_points[9][1] = gp1;
			gauss_points[10][0] = gp2;
			gauss_points[10][1] = gp2;
			gauss_points[11][0] = gp2;
			gauss_points[11][1] = gp3;

			gauss_points[12][0] = gp3;
			gauss_points[12][1] = gp0;
			gauss_points[13][0] = gp3;
			gauss_points[13][1] = gp1;
			gauss_points[14][0] = gp3;
			gauss_points[14][1] = gp2;
			gauss_points[15][0] = gp3;
			gauss_points[15][1] = gp3;

			gauss_weights[0] = w0 * w0;
			gauss_weights[1] = w0 * w1;
			gauss_weights[2] = w0 * w2;
			gauss_weights[3] = w0 * w3;

			gauss_weights[4] = w1 * w0;
			gauss_weights[5] = w1 * w1;
			gauss_weights[6] = w1 * w2;
			gauss_weights[7] = w1 * w3;

			gauss_weights[8] = w2 * w0;
			gauss_weights[9] = w2 * w1;
			gauss_weights[10] = w2 * w2;
			gauss_weights[11] = w2 * w3;

			gauss_weights[12] = w3 * w0;
			gauss_weights[13] = w3 * w1;
			gauss_weights[14] = w3 * w2;
			gauss_weights[15] = w3 * w3;

			break;
		}
		case 25: // P(9)
		{
			double gp0 = -1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
			double gp1 = -1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
			double gp2 = 0.0;
			double gp3 = 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
			double gp4 = 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));

			double w0 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
			double w1 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
			double w2 = 128.0 / 225.0;
			double w3 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
			double w4 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;

			gauss_points[0][0] = gp0;
			gauss_points[0][1] = gp0;
			gauss_points[1][0] = gp0;
			gauss_points[1][1] = gp1;
			gauss_points[2][0] = gp0;
			gauss_points[2][1] = gp2;
			gauss_points[3][0] = gp0;
			gauss_points[3][1] = gp3;
			gauss_points[4][0] = gp0;
			gauss_points[4][1] = gp4;

			gauss_points[5][0] = gp1;
			gauss_points[5][1] = gp0;
			gauss_points[6][0] = gp1;
			gauss_points[6][1] = gp1;
			gauss_points[7][0] = gp1;
			gauss_points[7][1] = gp2;
			gauss_points[8][0] = gp1;
			gauss_points[8][1] = gp3;
			gauss_points[9][0] = gp1;
			gauss_points[9][1] = gp4;

			gauss_points[10][0] = gp2;
			gauss_points[10][1] = gp0;
			gauss_points[11][0] = gp2;
			gauss_points[11][1] = gp1;
			gauss_points[12][0] = gp2;
			gauss_points[12][1] = gp2;
			gauss_points[13][0] = gp2;
			gauss_points[13][1] = gp3;
			gauss_points[14][0] = gp2;
			gauss_points[14][1] = gp4;

			gauss_points[15][0] = gp3;
			gauss_points[15][1] = gp0;
			gauss_points[16][0] = gp3;
			gauss_points[16][1] = gp1;
			gauss_points[17][0] = gp3;
			gauss_points[17][1] = gp2;
			gauss_points[18][0] = gp3;
			gauss_points[18][1] = gp3;
			gauss_points[19][0] = gp3;
			gauss_points[19][1] = gp4;

			gauss_points[20][0] = gp4;
			gauss_points[20][1] = gp0;
			gauss_points[21][0] = gp4;
			gauss_points[21][1] = gp1;
			gauss_points[22][0] = gp4;
			gauss_points[22][1] = gp2;
			gauss_points[23][0] = gp4;
			gauss_points[23][1] = gp3;
			gauss_points[24][0] = gp4;
			gauss_points[24][1] = gp4;

			gauss_weights[0] = w0 * w0;
			gauss_weights[1] = w0 * w1;
			gauss_weights[2] = w0 * w2;
			gauss_weights[3] = w0 * w3;
			gauss_weights[4] = w0 * w4;

			gauss_weights[5] = w1 * w0;
			gauss_weights[6] = w1 * w1;
			gauss_weights[7] = w1 * w2;
			gauss_weights[8] = w1 * w3;
			gauss_weights[9] = w1 * w4;

			gauss_weights[10] = w2 * w0;
			gauss_weights[11] = w2 * w1;
			gauss_weights[12] = w2 * w2;
			gauss_weights[13] = w2 * w3;
			gauss_weights[14] = w2 * w4;

			gauss_weights[15] = w3 * w0;
			gauss_weights[16] = w3 * w1;
			gauss_weights[17] = w3 * w2;
			gauss_weights[18] = w3 * w3;
			gauss_weights[19] = w3 * w4;

			gauss_weights[20] = w4 * w0;
			gauss_weights[21] = w4 * w1;
			gauss_weights[22] = w4 * w2;
			gauss_weights[23] = w4 * w3;
			gauss_weights[24] = w4 * w4;

			break;
		}
		default:
		{
			// Can not be here!
			exit(-1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature for triangle
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note:
	 */
	inline void GaussTriangle(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 1.0 / 3.0;
			gauss_points[0][1] = 1.0 / 3.0;
			gauss_points[0][2] = 1.0 / 3.0;

			gauss_weights[0] = 1.0;

			break;
		}
		case 3: //(P2)
		{
			// gauss_points[0][0] = 0.0; gauss_points[0][1] = 0.5; gauss_points[0][2] = 0.5;
			// gauss_points[1][0] = 0.5; gauss_points[1][1] = 0.0; gauss_points[1][2] = 0.5;
			// gauss_points[2][0] = 0.5; gauss_points[2][1] = 0.5; gauss_points[2][2] = 0.0;

			gauss_points[0][0] = 2.0 / 3.0;
			gauss_points[0][1] = 1.0 / 6.0;
			gauss_points[0][2] = 1.0 / 6.0;
			gauss_points[1][0] = 1.0 / 6.0;
			gauss_points[1][1] = 2.0 / 3.0;
			gauss_points[1][2] = 1.0 / 6.0;
			gauss_points[2][0] = 1.0 / 6.0;
			gauss_points[2][1] = 1.0 / 6.0;
			gauss_points[2][2] = 2.0 / 3.0;

			gauss_weights[0] = 1.0 / 3.0;
			gauss_weights[1] = 1.0 / 3.0;
			gauss_weights[2] = 1.0 / 3.0;

			break;
		}
		case 4: //(P3)
		{
			gauss_points[0][0] = 1.0 / 3.0;
			gauss_points[0][1] = 1.0 / 3.0;
			gauss_points[0][2] = 1.0 / 3.0;

			gauss_points[1][0] = 0.6;
			gauss_points[1][1] = 0.2;
			gauss_points[1][2] = 0.2;
			gauss_points[2][0] = 0.2;
			gauss_points[2][1] = 0.6;
			gauss_points[2][2] = 0.2;
			gauss_points[3][0] = 0.2;
			gauss_points[3][1] = 0.2;
			gauss_points[3][2] = 0.6;

			gauss_weights[0] = -27.0 / 48.0;

			gauss_weights[1] = 25.0 / 48.0;
			gauss_weights[2] = 25.0 / 48.0;
			gauss_weights[3] = 25.0 / 48.0;

			break;
		}
		case 7: //(P5)
		{
			double alpha1 = 0.059715871789770;
			double alpha2 = 0.797426985353087;
			double beta1 = 0.470142064105115;
			double beta2 = 0.101286507323456;

			gauss_points[0][0] = 1.0 / 3.0;
			gauss_points[0][1] = 1.0 / 3.0;
			gauss_points[0][2] = 1.0 / 3.0;

			gauss_points[1][0] = alpha1;
			gauss_points[1][1] = beta1;
			gauss_points[1][2] = beta1;
			gauss_points[2][0] = beta1;
			gauss_points[2][1] = alpha1;
			gauss_points[2][2] = beta1;
			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = alpha1;

			gauss_points[4][0] = alpha2;
			gauss_points[4][1] = beta2;
			gauss_points[4][2] = beta2;
			gauss_points[5][0] = beta2;
			gauss_points[5][1] = alpha2;
			gauss_points[5][2] = beta2;
			gauss_points[6][0] = beta2;
			gauss_points[6][1] = beta2;
			gauss_points[6][2] = alpha2;

			gauss_weights[0] = 0.225000000000000;

			gauss_weights[1] = 0.132394152788506;
			gauss_weights[2] = 0.132394152788506;
			gauss_weights[3] = 0.132394152788506;

			gauss_weights[4] = 0.125939180544827;
			gauss_weights[5] = 0.125939180544827;
			gauss_weights[6] = 0.125939180544827;

			break;
		}
		case 13: //(P7)
		{
			double beta1 = 0.260345966079038;
			double alpha1 = 0.479308067841923;

			double beta2 = 0.065130102902216;
			double alpha2 = 0.869739794195568;

			double alpha3 = 0.638444188569809;
			double beta3 = 0.312865496004875;
			double gamma3 = 0.048690315425316;

			gauss_points[0][0] = 1.0 / 3.0;
			gauss_points[0][1] = 1.0 / 3.0;
			gauss_points[0][2] = 1.0 / 3.0;

			gauss_points[1][0] = alpha1;
			gauss_points[1][1] = beta1;
			gauss_points[1][2] = beta1;
			gauss_points[2][0] = beta1;
			gauss_points[2][1] = alpha1;
			gauss_points[2][2] = beta1;
			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = alpha1;

			gauss_points[4][0] = alpha2;
			gauss_points[4][1] = beta2;
			gauss_points[4][2] = beta2;
			gauss_points[5][0] = beta2;
			gauss_points[5][1] = alpha2;
			gauss_points[5][2] = beta2;
			gauss_points[6][0] = beta2;
			gauss_points[6][1] = beta2;
			gauss_points[6][2] = alpha2;

			gauss_points[7][0] = alpha3;
			gauss_points[7][1] = beta3;
			gauss_points[7][2] = gamma3;
			gauss_points[8][0] = alpha3;
			gauss_points[8][1] = gamma3;
			gauss_points[8][2] = beta3;
			gauss_points[9][0] = beta3;
			gauss_points[9][1] = alpha3;
			gauss_points[9][2] = gamma3;
			gauss_points[10][0] = beta3;
			gauss_points[10][1] = gamma3;
			gauss_points[10][2] = alpha3;
			gauss_points[11][0] = gamma3;
			gauss_points[11][1] = alpha3;
			gauss_points[11][2] = beta3;
			gauss_points[12][0] = gamma3;
			gauss_points[12][1] = beta3;
			gauss_points[12][2] = alpha3;

			gauss_weights[0] = -0.149570044467670;

			gauss_weights[1] = 0.175615257433204;
			gauss_weights[2] = 0.175615257433204;
			gauss_weights[3] = 0.175615257433204;

			gauss_weights[4] = 0.053347235608839;
			gauss_weights[5] = 0.053347235608839;
			gauss_weights[6] = 0.053347235608839;

			gauss_weights[7] = 0.077113760890257;
			gauss_weights[8] = 0.077113760890257;
			gauss_weights[9] = 0.077113760890257;
			gauss_weights[10] = 0.077113760890257;
			gauss_weights[11] = 0.077113760890257;
			gauss_weights[12] = 0.077113760890257;
			break;
		}
		default:
		{
			// Can not be here!
			exit(-1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature for standard tetrahedron element
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: P Keast, Moderate degree tetrahedral quadrature formulas, CMAME 55: 339-348 (1986)
	 */
	inline void GaussTetrahedron(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P0)
		{
			gauss_points[0][0] = 1.0 / 4.0;
			gauss_points[0][1] = 1.0 / 4.0;
			gauss_points[0][2] = 1.0 / 4.0;
			gauss_points[0][3] = 1.0 / 4.0;

			gauss_weights[0] = 1.0;

			break;
		}
		case 4: //(P1)
		{
			double alpha1 = 0.5854101966249685;
			double beta1 = 0.1381966011250105;

			gauss_points[0][0] = alpha1;
			gauss_points[0][1] = beta1;
			gauss_points[0][2] = beta1;
			gauss_points[0][3] = beta1;

			gauss_points[1][0] = beta1;
			gauss_points[1][1] = alpha1;
			gauss_points[1][2] = beta1;
			gauss_points[1][3] = beta1;

			gauss_points[2][0] = beta1;
			gauss_points[2][1] = beta1;
			gauss_points[2][2] = alpha1;
			gauss_points[2][3] = beta1;

			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = beta1;
			gauss_points[3][3] = alpha1;

			gauss_weights[0] = 1.0 / 4.0;
			gauss_weights[1] = 1.0 / 4.0;
			gauss_weights[2] = 1.0 / 4.0;
			gauss_weights[3] = 1.0 / 4.0;

			break;
		}
		case 5: //(P2)
		{
			double alpha1 = 1.0 / 2.0;
			double beta1 = 1.0 / 6.0;

			gauss_points[0][0] = alpha1;
			gauss_points[0][1] = beta1;
			gauss_points[0][2] = beta1;
			gauss_points[0][3] = beta1;

			gauss_points[1][0] = beta1;
			gauss_points[1][1] = alpha1;
			gauss_points[1][2] = beta1;
			gauss_points[1][3] = beta1;

			gauss_points[2][0] = beta1;
			gauss_points[2][1] = beta1;
			gauss_points[2][2] = alpha1;
			gauss_points[2][3] = beta1;

			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = beta1;
			gauss_points[3][3] = alpha1;

			gauss_points[4][0] = 1.0 / 4.0;
			gauss_points[4][1] = 1.0 / 4.0;
			gauss_points[4][2] = 1.0 / 4.0;
			gauss_points[4][3] = 1.0 / 4.0;

			gauss_weights[0] = 9.0 / 20.0;
			gauss_weights[1] = 9.0 / 20.0;
			gauss_weights[2] = 9.0 / 20.0;
			gauss_weights[3] = 9.0 / 20.0;
			gauss_weights[4] = -4.0 / 5.0;

			break;
		}
		case 11: //(P4)
		{
			double alpha1 = 0.7857142857142857;
			double beta1 = 0.0714285714285714;
			double alpha2 = 0.3994035761667992;
			double beta2 = 0.1005964238332008;

			// 4
			gauss_points[0][0] = alpha1;
			gauss_points[0][1] = beta1;
			gauss_points[0][2] = beta1;
			gauss_points[0][3] = beta1;

			gauss_points[1][0] = beta1;
			gauss_points[1][1] = alpha1;
			gauss_points[1][2] = beta1;
			gauss_points[1][3] = beta1;

			gauss_points[2][0] = beta1;
			gauss_points[2][1] = beta1;
			gauss_points[2][2] = alpha1;
			gauss_points[2][3] = beta1;

			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = beta1;
			gauss_points[3][3] = alpha1;

			// 6
			gauss_points[4][0] = alpha2;
			gauss_points[4][1] = alpha2;
			gauss_points[4][2] = beta2;
			gauss_points[4][3] = beta2;

			gauss_points[5][0] = alpha2;
			gauss_points[5][1] = beta2;
			gauss_points[5][2] = alpha2;
			gauss_points[5][3] = beta2;

			gauss_points[6][0] = alpha2;
			gauss_points[6][1] = beta2;
			gauss_points[6][2] = beta2;
			gauss_points[6][3] = alpha2;

			gauss_points[7][0] = beta2;
			gauss_points[7][1] = alpha2;
			gauss_points[7][2] = alpha2;
			gauss_points[7][3] = beta2;

			gauss_points[8][0] = beta2;
			gauss_points[8][1] = alpha2;
			gauss_points[8][2] = beta2;
			gauss_points[8][3] = alpha2;

			gauss_points[9][0] = beta2;
			gauss_points[9][1] = beta2;
			gauss_points[9][2] = alpha2;
			gauss_points[9][3] = alpha2;

			// 1
			gauss_points[10][0] = 1.0 / 4.0;
			gauss_points[10][1] = 1.0 / 4.0;
			gauss_points[10][2] = 1.0 / 4.0;
			gauss_points[10][3] = 1.0 / 4.0;

			gauss_weights[0] = 0.0457333333333333;
			gauss_weights[1] = 0.0457333333333333;
			gauss_weights[2] = 0.0457333333333333;
			gauss_weights[3] = 0.0457333333333333;

			gauss_weights[4] = 0.1493333333333333;
			gauss_weights[5] = 0.1493333333333333;
			gauss_weights[6] = 0.1493333333333333;
			gauss_weights[7] = 0.1493333333333333;
			gauss_weights[8] = 0.1493333333333333;
			gauss_weights[9] = 0.1493333333333333;

			gauss_weights[10] = -0.0789333333333333;

			break;
		}
		case 15: // P(6)
		{
			double alpha1 = 0.0;
			double beta1 = 1.0 / 3.0;
			double alpha2 = 0.7272727272727272;
			double beta2 = 0.0909090909090909;
			double alpha3 = 0.0665501535736643;
			double beta3 = 0.4334498464263357;

			// 4
			gauss_points[0][0] = alpha1;
			gauss_points[0][1] = beta1;
			gauss_points[0][2] = beta1;
			gauss_points[0][3] = beta1;

			gauss_points[1][0] = beta1;
			gauss_points[1][1] = alpha1;
			gauss_points[1][2] = beta1;
			gauss_points[1][3] = beta1;

			gauss_points[2][0] = beta1;
			gauss_points[2][1] = beta1;
			gauss_points[2][2] = alpha1;
			gauss_points[2][3] = beta1;

			gauss_points[3][0] = beta1;
			gauss_points[3][1] = beta1;
			gauss_points[3][2] = beta1;
			gauss_points[3][3] = alpha1;

			// 4
			gauss_points[4][0] = alpha2;
			gauss_points[4][1] = beta2;
			gauss_points[4][2] = beta2;
			gauss_points[4][3] = beta2;

			gauss_points[5][0] = beta2;
			gauss_points[5][1] = alpha2;
			gauss_points[5][2] = beta2;
			gauss_points[5][3] = beta2;

			gauss_points[6][0] = beta2;
			gauss_points[6][1] = beta2;
			gauss_points[6][2] = alpha2;
			gauss_points[6][3] = beta2;

			gauss_points[7][0] = beta2;
			gauss_points[7][1] = beta2;
			gauss_points[7][2] = beta2;
			gauss_points[7][3] = alpha2;

			// 6
			gauss_points[8][0] = alpha3;
			gauss_points[8][1] = alpha3;
			gauss_points[8][2] = beta3;
			gauss_points[8][3] = beta3;

			gauss_points[9][0] = alpha3;
			gauss_points[9][1] = beta3;
			gauss_points[9][2] = alpha3;
			gauss_points[9][3] = beta3;

			gauss_points[10][0] = alpha3;
			gauss_points[10][1] = beta3;
			gauss_points[10][2] = beta3;
			gauss_points[10][3] = alpha3;

			gauss_points[11][0] = beta3;
			gauss_points[11][1] = alpha3;
			gauss_points[11][2] = alpha3;
			gauss_points[11][3] = beta3;

			gauss_points[12][0] = beta3;
			gauss_points[12][1] = alpha3;
			gauss_points[12][2] = beta3;
			gauss_points[12][3] = alpha3;

			gauss_points[13][0] = beta3;
			gauss_points[13][1] = beta3;
			gauss_points[13][2] = alpha3;
			gauss_points[13][3] = alpha3;

			// 1
			gauss_points[14][0] = 1.0 / 4.0;
			gauss_points[14][1] = 1.0 / 4.0;
			gauss_points[14][2] = 1.0 / 4.0;
			gauss_points[14][3] = 1.0 / 4.0;

			gauss_weights[0] = 0.0361607142857143;
			gauss_weights[1] = 0.0361607142857143;
			gauss_weights[2] = 0.0361607142857143;
			gauss_weights[3] = 0.0361607142857143;

			gauss_weights[4] = 0.0698714945161738;
			gauss_weights[5] = 0.0698714945161738;
			gauss_weights[6] = 0.0698714945161738;
			gauss_weights[7] = 0.0698714945161738;

			gauss_weights[8] = 0.0656948493683187;
			gauss_weights[9] = 0.0656948493683187;
			gauss_weights[10] = 0.0656948493683187;
			gauss_weights[11] = 0.0656948493683187;
			gauss_weights[12] = 0.0656948493683187;
			gauss_weights[13] = 0.0656948493683187;

			gauss_weights[14] = 0.1817020685825351;

			break;
		}
		default:
		{
			cout << "Error: Gauss quadrature rule is not supported for Tetrahedron!" << endl;
			cout << "File: " << __FILE__ << endl;
			cout << "Line: " << __LINE__ << endl;
			exit(1);
		}
		}
	}

	/**
	 * Gaussian elimination for solving systems of linear equations (in C style)
	 * \param[in] N : dimension of linear equations
	 * \param[in] A : 1D left-hand matrix of Ax=B (note: A is stored in one dimension)
	 * \param[in|out] B: 1D right-hand vector/solution of Ax=B
	 */
	inline void GaussElimination(int N, double *A, double *B)
	{
		// 分配内存空间
		double **a;
		int i, j, k, n;
		a = new double *[N];
		for (i = 0; i != N; i++)
		{
			a[i] = new double[N];
		}
		double *b;
		b = new double[N];
		double *x;
		x = new double[N];

		// 赋值
		for (i = 0; i != N; i++)
		{
			for (j = 0; j != N; j++)
			{
				a[i][j] = A[i * N + j];
			}
		}

		for (j = 0; j != N; j++)
		{
			b[j] = B[j];
		}

		// 选列主元Gauss消去法求解线性方程组
		int num;
		double major, temp, sum;
		double *m;
		m = new double[N];

		for (k = 0; k != N - 1; ++k)
		{
			// step-1:选主元
			num = k;
			major = a[k][k];
			for (i = k + 1; i != N; i++)
			{
				if (fabs(a[i][k]) > fabs(major))
				{
					num = i;
					major = a[i][k];
				}
			}

			// step-2:行交换
			for (j = k; j != N; j++)
			{
				temp = a[k][j];
				a[k][j] = a[num][j];
				a[num][j] = temp;
			}
			temp = b[k];
			b[k] = b[num];
			b[num] = temp;

			// step-3:Gauss消元法
			for (i = k + 1; i != N; i++)
			{
				m[i] = a[i][k] / a[k][k];

				for (j = k + 1; j != N; j++)
				{
					a[i][j] = a[i][j] - m[i] * a[k][j];
				}

				b[i] = b[i] - m[i] * b[k];
			}
		}

		// step-4: 回代过程
		x[N - 1] = b[N - 1] / a[N - 1][N - 1];
		for (n = N - 2; n != -1; n--)
		{
			sum = 0.0;
			for (j = n + 1; j != N; j++)
			{
				sum = sum + a[n][j] * x[j];
			}

			x[n] = (b[n] - sum) / a[n][n];
		}

		// 返回解向量
		for (i = 0; i != N; i++)
		{
			B[i] = x[i];
		}

		delete[] m;
		delete[] b;
		delete[] x;
		for (i = 0; i != N; i++)
		{
			delete[] a[i];
		}
		delete[] a;
	}

	/**
	 * Gaussian elimination for solving systems of linear equations (in C++ style)
	 * \param[in] N : dimension of linear equations
	 * \param[in] A : 2D left-hand matrix of Ax=B
	 * \param[in|out] B: 1D right-hand vector/solution of Ax=B
	 */
	inline void GaussElimination(int N, const Array2D<double> &A, Array1D<double> &B)
	{
		// 分配内存空间
		int i(0), j(0), k(0), n(0);
		Array2D<double> a(N, N);
		Array1D<double> b(N);
		Array1D<double> x(N);

		// 赋值
		for (i = 0; i != N; i++)
		{
			for (j = 0; j != N; j++)
			{
				a[i][j] = A[i][j];
			}
		}

		for (j = 0; j != N; j++)
		{
			b[j] = B[j];
		}

		// 选列主元Gauss消去法求解线性方程组
		int num(0);
		double temp(0.0), sum(0.0);
		Array1D<double> m(N);

		for (k = 0; k != N - 1; ++k)
		{
			// step-1:选主元
			num = k;
			double major = a[k][k];
			for (i = k + 1; i != N; i++)
			{
				if (fabs(a[i][k]) > fabs(major))
				{
					num = i;
					major = a[i][k];
				}
			}

			// step-2:行交换
			for (j = k; j != N; j++)
			{
				temp = a[k][j];
				a[k][j] = a[num][j];
				a[num][j] = temp;
			}
			temp = b[k];
			b[k] = b[num];
			b[num] = temp;

			// step-3:Gauss消元法
			for (i = k + 1; i != N; i++)
			{
				m[i] = a[i][k] / a[k][k];

				for (j = k + 1; j != N; j++)
				{
					a[i][j] = a[i][j] - m[i] * a[k][j];
				}

				b[i] = b[i] - m[i] * b[k];
			}
		}

		// step-4: 回代过程
		x[N - 1] = b[N - 1] / a[N - 1][N - 1];
		for (n = N - 2; n != -1; n--)
		{
			sum = 0.0;
			for (j = n + 1; j != N; j++)
			{
				sum = sum + a[n][j] * x[j];
			}

			x[n] = (b[n] - sum) / a[n][n];
		}

		// 返回解向量
		for (i = 0; i != N; i++)
		{
			B[i] = x[i];
		}
	}

	/**
	 * Get inverse matrix of array A by Gaussian elimination
	 * \param[in] N : dimension of array
	 * \param[in] A : array A
	 * \param[in|out] B: inverse matrix of A
	 */
	inline void getInverseMatrix_GE(int N, const Array2D<double> &A, Array2D<double> &B)
	{
		Array1D<double> temp(N);

		for (int k = 0; k != N; ++k)
		{
			// initialize right hand vector
			for (int i = 0; i != N; ++i)
			{
				temp[i] = 0.0;
			}
			temp[k] = 1.0;

			// solve the k_th column of inverse matrix
			GaussElimination(N, A, temp);

			// set the k_th column of the inverse matrix
			for (int i = 0; i != N; ++i)
			{
				B[i][k] = temp[i];
			}
		}
	}

	/**
	 * Partial pivot LU decomposition for suqare matrix using Doolittle algorithm
	 *
	 * \param[in] N : dim of A
	 * \param[in] A :
	 * \param[out] L,U,Q: QA = LU
	 */
	inline void LUDecomposition(int N, Array2D<double> &A, Array2D<double> &L, Array2D<double> &U, Array2D<double> &Q)
	{
		// Set Q
		L.setZero();
		U.setZero();
		Q.setZero();
		for (int i = 0; i != N; ++i)
		{
			Q[i][i] = 1;
		}

		// temporal variables
		Array1D<double> s(N);
		double temp(0);
		int ik(0);

		for (int k = 0; k < N; ++k)
		{
			for (int i = k; i < N; ++i)
			{
				s[i] = A[i][k];
				for (int t = 0; t <= k - 1; ++t)
				{
					s[i] = s[i] - L[i][t] * U[t][k];
				}
			}

			temp = fabs(s[k]);
			ik = k;
			for (int i = k + 1; i < N; ++i)
			{
				if (fabs(s[i]) > temp)
				{
					temp = fabs(s[i]);
					ik = i;
				}
			}

			if (ik != k) // exchange rows
			{
				temp = s[k];
				s[k] = s[ik];
				s[ik] = temp;

				for (int t = 0; t <= k - 1; ++t)
				{
					temp = L[k][t];
					L[k][t] = L[ik][t];
					L[ik][t] = temp;
				}
				for (int t = k; t < N; ++t)
				{
					temp = A[k][t];
					A[k][t] = A[ik][t];
					A[ik][t] = temp;
				}

				for (int t = 0; t < N; ++t)
				{
					temp = Q[k][t];
					Q[k][t] = Q[ik][t];
					Q[ik][t] = temp;
				}
			}

			U[k][k] = s[k];

			for (int j = k + 1; j < N; ++j)
			{
				U[k][j] = A[k][j];
				for (int t = 0; t <= k - 1; ++t)
				{
					U[k][j] = U[k][j] - L[k][t] * U[t][j];
				}
			}

			L[k][k] = 1.0;
			for (int i = k + 1; i < N; ++i)
			{
				L[i][k] = s[i] / U[k][k];
			}
		}
	}

	/**
	 * Get inverse matrix of array A by LU decomposition
	 * \param[in] N : dimension of array
	 * \param[in] A : array A
	 * \param[in|out] B: inverse matrix of A
	 */
	inline void getInverseMatrix_LU(int N, const Array2D<double> &A, Array2D<double> &B)
	{
		Array2D<double> L(N, N, 0.0);
		Array2D<double> U(N, N, 0.0);
		Array2D<double> Q(N, N, 0.0);
		Array2D<double> T(N, N, 0.0);
		Array1D<double> y(N);

		for (int i = 0; i != N; ++i)
		{
			for (int j = 0; j != N; ++j)
			{
				T[i][j] = A[i][j];
			}
		}

		// Step-1: QA=LU
		LUDecomposition(N, T, L, U, Q);

		// Step-2: B=inv(A)
		B.setZero();
		for (int k = 0; k != N; ++k)
		{
			// forward iteration
			y[0] = Q[0][k];
			for (int i = 1; i != N; ++i)
			{
				y[i] = Q[i][k];
				for (int t = 0; t <= i - 1; ++t)
				{
					y[i] = y[i] - L[i][t] * y[t];
				}
			}

			// backward iteration
			B[N - 1][k] = y[N - 1] / U[N - 1][N - 1];
			for (int i = N - 2; i != -1; --i)
			{
				B[i][k] = y[i];
				for (int t = i + 1; t <= N - 1; ++t)
				{
					B[i][k] = B[i][k] - U[i][t] * B[t][k];
				}
				B[i][k] = B[i][k] / U[i][i];
			}
		}
	}

	/**
	 * Get the inverse of square matrix by recursive Schur Complement
	 * \param[in] N : dimension of array
	 * \param[in] M : matrix M
	 * \param[in|out] IM: inverse matrix of M
	 * Author: J. Cheng, F. Zhang(Who provides the algorithm in matlab)
	 * Reference: https://en.wikipedia.org/wiki/Schur_complement
	 * Date: Aug. 17, 2016
	 */
	inline void getInverseMatrix_SC(int N, const Array2D<double> &M, Array2D<double> &IM)
	{
		IM.setZero();

		if (N == 1)
		{
			IM[0][0] = 1.0 / M[0][0];

			if (fabs(M[0][0]) < 1.0e-16)
			{
				std::cout << "Warning: Divided by machine zero, Schur complement may fail..." << std::endl;
			}
		}
		else if (N == 2)
		{
			double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];

			if (fabs(det) < 1.0e-16)
			{
				std::cout << "Warning: Divided by machine zero, Schur complement may fail..." << std::endl;
			}

			IM[0][0] = M[1][1] / det;
			IM[1][1] = M[0][0] / det;
			IM[0][1] = -M[0][1] / det;
			IM[1][0] = -M[1][0] / det;
		}
		else // N > 2
		{
			int p = N / 2;
			int q = N - p;

			Array2D<double> A(p, p), B(p, q), C(q, p), D(q, q);

			// A[p,p]
			for (int i = 0; i != p; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					A[i][j] = M[i][j];
				}
			}

			// B[p,q]
			for (int i = 0; i != p; ++i)
			{
				for (int j = p; j != N; ++j)
				{
					B[i][j - p] = M[i][j];
				}
			}

			// C[q,p]
			for (int i = p; i != N; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					C[i - p][j] = M[i][j];
				}
			}

			// D[q,q]
			for (int i = p; i != N; ++i)
			{
				for (int j = p; j != N; ++j)
				{
					D[i - p][j - p] = M[i][j];
				}
			}

			// IA[p,p]
			Array2D<double> IA(p, p);
			getInverseMatrix_SC(p, A, IA);

			// IAB[p,q]
			Array2D<double> IAB(p, q, 0);
			for (int i = 0; i != p; ++i)
			{
				for (int j = 0; j != q; ++j)
				{
					// IAB[i][j]
					for (int k = 0; k != p; ++k)
					{
						IAB[i][j] += IA[i][k] * B[k][j];
					}
				}
			}

			// CIA[q,p]
			Array2D<double> CIA(q, p, 0);
			for (int i = 0; i != q; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					// CIA[i][j]
					for (int k = 0; k != p; ++k)
					{
						CIA[i][j] += C[i][k] * IA[k][j];
					}
				}
			}

			// D_CIAB[q][q]
			Array2D<double> D_CIAB(q, q, 0);
			for (int i = 0; i != q; ++i)
			{
				for (int j = 0; j != q; ++j)
				{
					D_CIAB[i][j] = D[i][j];

					// D_CIAB[i][j]
					for (int k = 0; k != p; ++k)
					{
						D_CIAB[i][j] -= CIA[i][k] * B[k][j];
					}
				}
			}

			Array2D<double> T(q, q);
			getInverseMatrix_SC(q, D_CIAB, T);

			// TCIA[q,p]
			Array2D<double> TCIA(q, p, 0);
			for (int i = 0; i != q; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					// TCIA[i][j]
					for (int k = 0; k != q; ++k)
					{
						TCIA[i][j] += T[i][k] * CIA[k][j];
					}
				}
			}

			// IM
			//(IA + IAB * TCIA)[p,p]
			for (int i = 0; i != p; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					IM[i][j] = IA[i][j];

					for (int k = 0; k != q; ++k)
					{
						IM[i][j] += IAB[i][k] * TCIA[k][j];
					}
				}
			}

			//(-IAB*T)[p,q]
			for (int i = 0; i != p; ++i)
			{
				for (int j = 0; j != q; ++j)
				{
					for (int k = 0; k != q; ++k)
					{
						IM[i][j + p] += -IAB[i][k] * T[k][j];
					}
				}
			}

			//-TCIA[q,p]
			for (int i = p; i != N; ++i)
			{
				for (int j = 0; j != p; ++j)
				{
					IM[i][j] = -TCIA[i - p][j];
				}
			}

			// T[q,q]
			for (int i = p; i != N; ++i)
			{
				for (int j = p; j != N; ++j)
				{
					IM[i][j] = T[i - p][j - p];
				}
			}
		}
	}
}
#endif
//////////////////////////////////////////////////////////////////////////
// End of Math.h