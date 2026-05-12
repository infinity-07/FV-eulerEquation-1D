//////////////////////////////////////////////////////////////////////////
/// Copyright(c) 2010-2015 CJBUAA
/// All right reserved.
///
/// @file: GaussQuadratureInRef.hpp
///
/// @brief: Gauss quadrature rule on 2D/3D reference element.
///
/// @author: Dr. Jian Cheng (chengjian@buaa.edu.cn)
/// @date: Dec 29, 2015
/// @modified: Aug 15, 2017, Jan, 01, 2018
//////////////////////////////////////////////////////////////////////////
#ifndef __GAUSSQUADRATUREINREF_HPP__
#define __GAUSSQUADRATUREINREF_HPP__
#include <iostream>
#include <fstream>
#include <string>

namespace sc_math
{
	///////////////////////////////////////////////////////////////////////////////////
	// 1D rules
	///////////////////////////////////////////////////////////////////////////////////
	/**
	 * Get points and weights of Gauss quadrature on reference 1D line [-1,1]
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Gauss-Legendre quadrature: https://github.com/vincentlab/PyFR/tree/develop/pyfr/quadrules/line
	 */
	inline void GaussLegendre_ref(int num_point, Array1D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		gauss_points.Resize(num_point);
		gauss_weights.Resize(num_point);
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0] = 0.0;
			gauss_weights[0] = 2.0;

			break;
		}
		case 2: //(P3)
		{
			gauss_points[0] = -0.577350269189625764509148780502;
			gauss_weights[0] = 1.0;
			gauss_points[1] = 0.577350269189625764509148780502;
			gauss_weights[1] = 1.0;

			break;
		}
		case 3: //(P5)
		{
			gauss_points[0] = -0.774596669241483377035853079956;
			gauss_weights[0] = 0.555555555555555555555555555556;
			gauss_points[1] = 0.0;
			gauss_weights[1] = 0.888888888888888888888888888889;
			gauss_points[2] = 0.774596669241483377035853079956;
			gauss_weights[2] = 0.555555555555555555555555555556;

			break;
		}
		case 4: //(P7)
		{
			gauss_points[0] = -0.861136311594052575223946488893;
			gauss_weights[0] = 0.347854845137453857373063949222;
			gauss_points[1] = -0.339981043584856264802665759103;
			gauss_weights[1] = 0.652145154862546142626936050778;
			gauss_points[2] = 0.339981043584856264802665759103;
			gauss_weights[2] = 0.652145154862546142626936050778;
			gauss_points[3] = 0.861136311594052575223946488893;
			gauss_weights[3] = 0.347854845137453857373063949222;

			break;
		}
		case 5: //(P9)
		{
			gauss_points[0] = -0.906179845938663992797626878299;
			gauss_weights[0] = 0.23692688505618908751426404072;
			gauss_points[1] = -0.5384693101056830910363144207;
			gauss_weights[1] = 0.478628670499366468041291514836;
			gauss_points[2] = 0.0;
			gauss_weights[2] = 0.568888888888888888888888888889;
			gauss_points[3] = 0.5384693101056830910363144207;
			gauss_weights[3] = 0.478628670499366468041291514836;
			gauss_points[4] = 0.906179845938663992797626878299;
			gauss_weights[4] = 0.23692688505618908751426404072;

			break;
		}
		case 6: //(P11)
		{
			gauss_points[0] = -0.932469514203152027812301554494;
			gauss_weights[0] = 0.171324492379170345040296142173;
			gauss_points[1] = -0.66120938646626451366139959502;
			gauss_weights[1] = 0.360761573048138607569833513838;
			gauss_points[2] = -0.238619186083196908630501721681;
			gauss_weights[2] = 0.46791393457269104738987034399;
			gauss_points[3] = 0.238619186083196908630501721681;
			gauss_weights[3] = 0.46791393457269104738987034399;
			gauss_points[4] = 0.66120938646626451366139959502;
			gauss_weights[4] = 0.360761573048138607569833513838;
			gauss_points[5] = 0.932469514203152027812301554494;
			gauss_weights[5] = 0.171324492379170345040296142173;

			break;
		}
		case 7: //(P13)
		{
			gauss_points[0] = -0.949107912342758524526189684048;
			gauss_weights[0] = 0.129484966168869693270611432679;
			gauss_points[1] = -0.741531185599394439863864773281;
			gauss_weights[1] = 0.279705391489276667901467771424;
			gauss_points[2] = -0.405845151377397166906606412077;
			gauss_weights[2] = 0.381830050505118944950369775489;
			gauss_points[3] = 0.0;
			gauss_weights[3] = 0.417959183673469387755102040816;
			gauss_points[4] = 0.405845151377397166906606412077;
			gauss_weights[4] = 0.381830050505118944950369775489;
			gauss_points[5] = 0.741531185599394439863864773281;
			gauss_weights[5] = 0.279705391489276667901467771424;
			gauss_points[6] = 0.949107912342758524526189684048;
			gauss_weights[6] = 0.129484966168869693270611432679;

			break;
		}
		case 8: //(P15)
		{
			gauss_points[0] = -0.960289856497536231683560868569;
			gauss_weights[0] = 0.10122853629037625915253135431;
			gauss_points[1] = -0.796666477413626739591553936476;
			gauss_weights[1] = 0.222381034453374470544355994426;
			gauss_points[2] = -0.525532409916328985817739049189;
			gauss_weights[2] = 0.313706645877887287337962201987;
			gauss_points[3] = -0.18343464249564980493947614236;
			gauss_weights[3] = 0.362683783378361982965150449277;
			gauss_points[4] = 0.18343464249564980493947614236;
			gauss_weights[4] = 0.362683783378361982965150449277;
			gauss_points[5] = 0.525532409916328985817739049189;
			gauss_weights[5] = 0.313706645877887287337962201987;
			gauss_points[6] = 0.796666477413626739591553936476;
			gauss_weights[6] = 0.222381034453374470544355994426;
			gauss_points[7] = 0.960289856497536231683560868569;
			gauss_weights[7] = 0.10122853629037625915253135431;

			break;
		}
		case 9: //(P17)
		{
			gauss_points[0] = -0.968160239507626089835576202904;
			gauss_weights[0] = 0.0812743883615744119718921581105;
			gauss_points[1] = -0.83603110732663579429942978807;
			gauss_weights[1] = 0.180648160694857404058472031243;
			gauss_points[2] = -0.613371432700590397308702039341;
			gauss_weights[2] = 0.260610696402935462318742869419;
			gauss_points[3] = -0.324253423403808929038538014643;
			gauss_weights[3] = 0.312347077040002840068630406584;
			gauss_points[4] = 0.0;
			gauss_weights[4] = 0.330239355001259763164525069287;
			gauss_points[5] = 0.324253423403808929038538014643;
			gauss_weights[5] = 0.312347077040002840068630406584;
			gauss_points[6] = 0.613371432700590397308702039341;
			gauss_weights[6] = 0.260610696402935462318742869419;
			gauss_points[7] = 0.83603110732663579429942978807;
			gauss_weights[7] = 0.180648160694857404058472031243;
			gauss_points[8] = 0.968160239507626089835576202904;
			gauss_weights[8] = 0.0812743883615744119718921581105;

			break;
		}
		case 10: //(P19)
		{
			gauss_points[0] = -0.973906528517171720077964012084;
			gauss_weights[0] = 0.0666713443086881375935688098933;
			gauss_points[1] = -0.865063366688984510732096688423;
			gauss_weights[1] = 0.149451349150580593145776339658;
			gauss_points[2] = -0.679409568299024406234327365115;
			gauss_weights[2] = 0.219086362515982043995534934228;
			gauss_points[3] = -0.433395394129247190799265943166;
			gauss_weights[3] = 0.269266719309996355091226921569;
			gauss_points[4] = -0.14887433898163121088482600113;
			gauss_weights[4] = 0.295524224714752870173892994651;
			gauss_points[5] = 0.14887433898163121088482600113;
			gauss_weights[5] = 0.295524224714752870173892994651;
			gauss_points[6] = 0.433395394129247190799265943166;
			gauss_weights[6] = 0.269266719309996355091226921569;
			gauss_points[7] = 0.679409568299024406234327365115;
			gauss_weights[7] = 0.219086362515982043995534934228;
			gauss_points[8] = 0.865063366688984510732096688423;
			gauss_weights[8] = 0.149451349150580593145776339658;
			gauss_points[9] = 0.973906528517171720077964012084;
			gauss_weights[9] = 0.0666713443086881375935688098933;

			break;
		}
		default:
		{
			std::cout << "Error: Gauss Legendre rule is not supported for 1D line!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}

	inline double GaussLocalPoint(double a, double b, double gpoint_ref)
	{
		double gpoint_local(0);
		gpoint_local = 0.5 * (a + b) + 0.5 * (b - a) * gpoint_ref;
		return gpoint_local;
	}

	inline double GaussLocalWeight(double a, double b, double gweight_ref)
	{
		double gweight_local(0);
		gweight_local = 0.5 * (a - b) * gweight_ref;
		return gweight_local;
	}

	/**
	 * Get points and weights of Gauss quadrature on reference 1D line [-1,1]
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Gauss-Lobatto quadrature: https://github.com/vincentlab/PyFR/tree/develop/pyfr/quadrules/line
	 */
	inline void GaussLobatto_ref(int num_point, Array1D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 2: //(P1)
		{
			gauss_points[0] = -1;
			gauss_weights[0] = 1.0;
			gauss_points[1] = 1;
			gauss_weights[1] = 1.0;

			break;
		}
		case 3: //(P3)
		{
			gauss_points[0] = -1;
			gauss_weights[0] = 0.333333333333333333333333333333;
			gauss_points[1] = 0;
			gauss_weights[1] = 1.333333333333333333333333333333;
			gauss_points[2] = 1;
			gauss_weights[2] = 0.333333333333333333333333333333;

			break;
		}
		case 4: //(P5)
		{
			gauss_points[0] = -1;
			gauss_weights[0] = 0.166666666666666666666666666667;
			gauss_points[1] = -0.447213595499957939281834733746;
			gauss_weights[1] = 0.833333333333333333333333333333;
			gauss_points[2] = 0.447213595499957939281834733746;
			gauss_weights[2] = 0.833333333333333333333333333333;
			gauss_points[3] = 1;
			gauss_weights[3] = 0.166666666666666666666666666667;

			break;
		}
		case 5: //(P7)
		{
			gauss_points[0] = -1;
			gauss_weights[0] = 0.1;
			gauss_points[1] = -0.654653670707977143798292456247;
			gauss_weights[1] = 0.544444444444444444444444444444;
			gauss_points[2] = 0.0;
			gauss_weights[2] = 0.711111111111111111111111111111;
			gauss_points[3] = 0.654653670707977143798292456247;
			gauss_weights[3] = 0.544444444444444444444444444444;
			gauss_points[4] = 1;
			gauss_weights[4] = 0.1;

			break;
		}
		case 6: //(P9)
		{
			gauss_points[0] = -1;
			gauss_weights[0] = 0.0666666666666666666666666666667;
			gauss_points[1] = -0.765055323929464692851002973959;
			gauss_weights[1] = 0.378474956297846980316612808212;
			gauss_points[2] = -0.285231516480645096314150994041;
			gauss_weights[2] = 0.554858377035486353016720525121;
			gauss_points[3] = 0.285231516480645096314150994041;
			gauss_weights[3] = 0.554858377035486353016720525121;
			gauss_points[4] = 0.765055323929464692851002973959;
			gauss_weights[4] = 0.378474956297846980316612808212;
			gauss_points[5] = 1;
			gauss_weights[5] = 0.0666666666666666666666666666667;

			break;
		}
		default:
		{
			std::cout << "Error: Gauss Lobatto rule is not supported for 1D line!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			std::cin.get();
			exit(1);
		}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////
	// 2D rules
	///////////////////////////////////////////////////////////////////////////////////
	/**
	 * Get points and weights of Gauss quadrature on reference 2D quadrilateral element [-1,1]x[-1,1]
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: generate based on 1D Gauss-Legendre rule
	 */
	inline void GaussQuadrilateral_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1:	  //(P1)
		case 4:	  //(P3)
		case 9:	  //(P5)
		case 16:  //(P7)
		case 25:  //(P9)
		case 36:  //(P11)
		case 49:  //(P13)
		case 64:  //(P15)
		case 81:  //(P17)
		case 100: //(P19)
		{
			int num_point_1d = (int)std::sqrt(num_point * 1.0);
			Array1D<double> gps(num_point_1d), gpw(num_point_1d);

			GaussLegendre_ref(num_point_1d, gps, gpw);

			for (int i = 0; i != num_point_1d; ++i)
			{
				for (int j = 0; j != num_point_1d; ++j)
				{
					int count = i * num_point_1d + j;

					gauss_points[count][0] = gps[i];
					gauss_points[count][1] = gps[j];

					gauss_weights[count] = gpw[i] * gpw[j];
				}
			}

			break;
		}
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 2D quadrilateral!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature on reference 2D triangle element [-1,-1]<->[1,-1]<->[-1,1]
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: https://github.com/vincentlab/PyFR/tree/develop/pyfr/quadrules/tri
	 */
	inline void GaussTriangle_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		std::string prefix_name = "witherden-vincent-";
		std::string suffix_name("");
		std::string file_name = "./scTools/tri/";
		fstream filein;

		switch (num_point)
		{
		case 1: //(P1)
			suffix_name = "n1-d1-sp.txt";
			break;
		case 3: //(P2)
			suffix_name = "n3-d2-sp.txt";
			break;
		case 6: //(P4)
			suffix_name = "n6-d4-sp.txt";
			break;
		case 7: //(P5)
			suffix_name = "n7-d5-sp.txt";
			break;
		case 12: //(P6)
			suffix_name = "n12-d6-sp.txt";
			break;
		case 15: //(P7)
			suffix_name = "n15-d7-sp.txt";
			break;
		case 16: //(P8)
			suffix_name = "n16-d8-sp.txt";
			break;
		case 19: //(P9)
			suffix_name = "n19-d9-sp.txt";
			break;
		case 25: //(P10)
			suffix_name = "n25-d10-sp.txt";
			break;
		case 28: //(P11)
			suffix_name = "n28-d11-sp.txt";
			break;
		case 33: //(P12)
			suffix_name = "n33-d12-sp.txt";
			break;
		case 37: //(P13)
			suffix_name = "n37-d13-sp.txt";
			break;
		case 42: //(P14)
			suffix_name = "n42-d14-sp.txt";
			break;
		case 49: //(P15)
			suffix_name = "n49-d15-sp.txt";
			break;
		case 55: //(P16)
			suffix_name = "n55-d16-sp.txt";
			break;
		case 60: //(P17)
			suffix_name = "n60-d17-sp.txt";
			break;
		case 67: //(P18)
			suffix_name = "n67-d18-sp.txt";
			break;
		case 73: //(P19)
			suffix_name = "n73-d19-sp.txt";
			break;
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 2D triangle!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}

		// set gauss points and weights via reading target data file
		file_name += (prefix_name + suffix_name);
		filein.open(file_name.c_str(), ios_base::in);
		if (!filein.is_open())
		{
			std::cout << "Error: Gauss Quadrature 2D triangle rules, can not load file: " << file_name << " !" << std::endl;
			exit(-1);
		}
		for (int n = 0; n != num_point; ++n)
		{
			filein >> gauss_points[n][0] >> gauss_points[n][1] >> gauss_weights[n];
		}
		filein.close();
	}

	///////////////////////////////////////////////////////////////////////////////////
	// 3D rules (under construction...)
	///////////////////////////////////////////////////////////////////////////////////
	/**
	 * Get points and weights of Gauss quadrature on reference 3D tetrahedron element
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Shape functions and points of integration of the finite element (v11370) (Code_Aster)
	 *        Libmesh: http://libmesh.github.io/
	 */
	inline void GaussTetrahedron_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P0-P1)
		{
			gauss_points[0][0] = 0.25;
			gauss_points[0][1] = 0.25;
			gauss_points[0][2] = 0.25;

			gauss_weights[0] = 1.0 / 6.0;

			break;
		}
		case 4: //(P2)
		{
			double a = (5.0 - sqrt(5.0)) / 20.0;
			double b = (5.0 + 3 * sqrt(5.0)) / 20.0;

			gauss_points[0][0] = a;
			gauss_points[0][1] = a;
			gauss_points[0][2] = a;
			gauss_points[1][0] = a;
			gauss_points[1][1] = a;
			gauss_points[1][2] = b;
			gauss_points[2][0] = a;
			gauss_points[2][1] = b;
			gauss_points[2][2] = a;
			gauss_points[3][0] = b;
			gauss_points[3][1] = a;
			gauss_points[3][2] = a;

			gauss_weights[0] = 1.0 / 24.0;
			gauss_weights[1] = 1.0 / 24.0;
			gauss_weights[2] = 1.0 / 24.0;
			gauss_weights[3] = 1.0 / 24.0;

			break;
		}
		case 5: //(P3)-Note: this rule has a negative weight and may be unsuitable for some problems.
		{
			double a = 0.25;
			double b = 1.0 / 6.0;
			double c = 0.5;

			gauss_points[0][0] = a;
			gauss_points[0][1] = a;
			gauss_points[0][2] = a;
			gauss_points[1][0] = b;
			gauss_points[1][1] = b;
			gauss_points[1][2] = b;
			gauss_points[2][0] = b;
			gauss_points[2][1] = b;
			gauss_points[2][2] = c;
			gauss_points[3][0] = b;
			gauss_points[3][1] = c;
			gauss_points[3][2] = b;
			gauss_points[4][0] = c;
			gauss_points[4][1] = b;
			gauss_points[4][2] = b;

			gauss_weights[0] = -2.0 / 15.0;
			gauss_weights[1] = 3.0 / 40.0;
			gauss_weights[2] = 3.0 / 40.0;
			gauss_weights[3] = 3.0 / 40.0;
			gauss_weights[4] = 3.0 / 40.0;

			break;
		}
		case 14: //(P5) Note: The 14-points fifth-order rule comes from libmesh (quadrature_gauss_3D.C of libmesh)
		{

			double a[3] = {0.31088591926330060980, 0.092735250310891226402, 0.045503704125649649492};
			double w[3] = {0.01878132095300264180, 0.012248840519393658257, 0.007091003462846911073};
			double b(0);

			//(1)points 1-4
			b = 1.0 - 3.0 * a[0];

			gauss_points[0][0] = a[0];
			gauss_points[0][1] = a[0];
			gauss_points[0][2] = a[0];
			gauss_points[1][0] = a[0];
			gauss_points[1][1] = b;
			gauss_points[1][2] = a[0];
			gauss_points[2][0] = b;
			gauss_points[2][1] = a[0];
			gauss_points[2][2] = a[0];
			gauss_points[3][0] = a[0];
			gauss_points[3][1] = a[0];
			gauss_points[3][2] = b;

			gauss_weights[0] = w[0];
			gauss_weights[1] = w[0];
			gauss_weights[2] = w[0];
			gauss_weights[3] = w[0];

			//(2)points 4-7
			b = 1.0 - 3.0 * a[1];

			gauss_points[4][0] = a[1];
			gauss_points[4][1] = a[1];
			gauss_points[4][2] = a[1];
			gauss_points[5][0] = a[1];
			gauss_points[5][1] = b;
			gauss_points[5][2] = a[1];
			gauss_points[6][0] = b;
			gauss_points[6][1] = a[1];
			gauss_points[6][2] = a[1];
			gauss_points[7][0] = a[1];
			gauss_points[7][1] = a[1];
			gauss_points[7][2] = b;

			gauss_weights[4] = w[1];
			gauss_weights[5] = w[1];
			gauss_weights[6] = w[1];
			gauss_weights[7] = w[1];

			//(3)points 8-14
			b = 0.5 * (1.0 - 2.0 * a[2]);

			gauss_points[8][0] = b;
			gauss_points[8][1] = b;
			gauss_points[8][2] = a[2];
			gauss_points[9][0] = b;
			gauss_points[9][1] = a[2];
			gauss_points[9][2] = a[2];
			gauss_points[10][0] = a[2];
			gauss_points[10][1] = a[2];
			gauss_points[10][2] = b;
			gauss_points[11][0] = a[2];
			gauss_points[11][1] = b;
			gauss_points[11][2] = a[2];
			gauss_points[12][0] = b;
			gauss_points[12][1] = a[2];
			gauss_points[12][2] = b;
			gauss_points[13][0] = a[2];
			gauss_points[13][1] = b;
			gauss_points[13][2] = b;

			gauss_weights[8] = w[2];
			gauss_weights[9] = w[2];
			gauss_weights[10] = w[2];
			gauss_weights[11] = w[2];
			gauss_weights[12] = w[2];
			gauss_weights[13] = w[2];

			break;
		}
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 3D tetrahedron!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature on reference 3D hexahedron element
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Shape functions and points of integration of the finite element (v11370) (Code_Aster)
	 */
	inline void GaussHexahedron_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 0.0;
			gauss_points[0][1] = 0.0;
			gauss_points[0][2] = 0.0;

			gauss_weights[0] = 8.0;

			break;
		}
		case 8: //(P3)
		{
			double alpha = 1.0 / sqrt(3.0);

			gauss_points[0][0] = -alpha;
			gauss_points[0][1] = -alpha;
			gauss_points[0][2] = -alpha;
			gauss_points[1][0] = -alpha;
			gauss_points[1][1] = alpha;
			gauss_points[1][2] = -alpha;
			gauss_points[2][0] = alpha;
			gauss_points[2][1] = -alpha;
			gauss_points[2][2] = -alpha;
			gauss_points[3][0] = alpha;
			gauss_points[3][1] = alpha;
			gauss_points[3][2] = -alpha;

			gauss_points[4][0] = -alpha;
			gauss_points[4][1] = -alpha;
			gauss_points[4][2] = alpha;
			gauss_points[5][0] = -alpha;
			gauss_points[5][1] = alpha;
			gauss_points[5][2] = alpha;
			gauss_points[6][0] = alpha;
			gauss_points[6][1] = -alpha;
			gauss_points[6][2] = alpha;
			gauss_points[7][0] = alpha;
			gauss_points[7][1] = alpha;
			gauss_points[7][2] = alpha;

			gauss_weights[0] = 1.0;
			gauss_weights[1] = 1.0;
			gauss_weights[2] = 1.0;
			gauss_weights[3] = 1.0;
			gauss_weights[4] = 1.0;
			gauss_weights[5] = 1.0;
			gauss_weights[6] = 1.0;
			gauss_weights[7] = 1.0;

			break;
		}
		case 27: //(P5)
		{
			double alpha = sqrt(3.0 / 5.0);
			double c12c2(200.0 / 729.0), c1c22(320.0 / 729.0);
			double c13(125.0 / 729.0), c23(512.0 / 729.0);

			gauss_points[0][0] = -alpha;
			gauss_points[0][1] = -alpha;
			gauss_points[0][2] = -alpha;
			gauss_points[1][0] = -alpha;
			gauss_points[1][1] = -alpha;
			gauss_points[1][2] = 0;
			gauss_points[2][0] = -alpha;
			gauss_points[2][1] = -alpha;
			gauss_points[2][2] = alpha;
			gauss_points[3][0] = -alpha;
			gauss_points[3][1] = 0;
			gauss_points[3][2] = -alpha;
			gauss_points[4][0] = -alpha;
			gauss_points[4][1] = 0;
			gauss_points[4][2] = 0;
			gauss_points[5][0] = -alpha;
			gauss_points[5][1] = 0;
			gauss_points[5][2] = alpha;
			gauss_points[6][0] = -alpha;
			gauss_points[6][1] = alpha;
			gauss_points[6][2] = -alpha;
			gauss_points[7][0] = -alpha;
			gauss_points[7][1] = alpha;
			gauss_points[7][2] = 0;
			gauss_points[8][0] = -alpha;
			gauss_points[8][1] = alpha;
			gauss_points[8][2] = alpha;

			gauss_points[9][0] = 0;
			gauss_points[9][1] = -alpha;
			gauss_points[9][2] = -alpha;
			gauss_points[10][0] = 0;
			gauss_points[10][1] = -alpha;
			gauss_points[10][2] = 0;
			gauss_points[11][0] = 0;
			gauss_points[11][1] = -alpha;
			gauss_points[11][2] = alpha;
			gauss_points[12][0] = 0;
			gauss_points[12][1] = 0;
			gauss_points[12][2] = -alpha;
			gauss_points[13][0] = 0;
			gauss_points[13][1] = 0;
			gauss_points[13][2] = 0;
			gauss_points[14][0] = 0;
			gauss_points[14][1] = 0;
			gauss_points[14][2] = alpha;
			gauss_points[15][0] = 0;
			gauss_points[15][1] = alpha;
			gauss_points[15][2] = -alpha;
			gauss_points[16][0] = 0;
			gauss_points[16][1] = alpha;
			gauss_points[16][2] = 0;
			gauss_points[17][0] = 0;
			gauss_points[17][1] = alpha;
			gauss_points[17][2] = alpha;

			gauss_points[18][0] = alpha;
			gauss_points[18][1] = -alpha;
			gauss_points[18][2] = -alpha;
			gauss_points[19][0] = alpha;
			gauss_points[19][1] = -alpha;
			gauss_points[19][2] = 0;
			gauss_points[20][0] = alpha;
			gauss_points[20][1] = -alpha;
			gauss_points[20][2] = alpha;
			gauss_points[21][0] = alpha;
			gauss_points[21][1] = 0;
			gauss_points[21][2] = -alpha;
			gauss_points[22][0] = alpha;
			gauss_points[22][1] = 0;
			gauss_points[22][2] = 0;
			gauss_points[23][0] = alpha;
			gauss_points[23][1] = 0;
			gauss_points[23][2] = alpha;
			gauss_points[24][0] = alpha;
			gauss_points[24][1] = alpha;
			gauss_points[24][2] = -alpha;
			gauss_points[25][0] = alpha;
			gauss_points[25][1] = alpha;
			gauss_points[25][2] = 0;
			gauss_points[26][0] = alpha;
			gauss_points[26][1] = alpha;
			gauss_points[26][2] = alpha;

			gauss_weights[0] = c13;
			gauss_weights[1] = c12c2;
			gauss_weights[2] = c13;

			gauss_weights[3] = c12c2;
			gauss_weights[4] = c1c22;
			gauss_weights[5] = c12c2;

			gauss_weights[6] = c13;
			gauss_weights[7] = c12c2;
			gauss_weights[8] = c13;

			gauss_weights[9] = c12c2;
			gauss_weights[10] = c1c22;
			gauss_weights[11] = c12c2;

			gauss_weights[12] = c1c22;
			gauss_weights[13] = c23;
			gauss_weights[14] = c1c22;

			gauss_weights[15] = c12c2;
			gauss_weights[16] = c1c22;
			gauss_weights[17] = c12c2;

			gauss_weights[18] = c13;
			gauss_weights[19] = c12c2;
			gauss_weights[20] = c13;

			gauss_weights[21] = c12c2;
			gauss_weights[22] = c1c22;
			gauss_weights[23] = c12c2;

			gauss_weights[24] = c13;
			gauss_weights[25] = c12c2;
			gauss_weights[26] = c13;

			break;
		}
		case 64:
		{
			double gp[4], gw[4];

			gp[0] = -8.6113631159405257522394648889281e-01L;
			gp[1] = -3.3998104358485626480266575910324e-01L;
			gp[2] = -gp[1];
			gp[3] = -gp[0];

			gw[0] = 3.4785484513745385737306394922200e-01L;
			gw[1] = 6.5214515486254614262693605077800e-01L;
			gw[2] = gw[1];
			gw[3] = gw[0];

			for (int i = 0; i != 4; ++i)
			{
				for (int j = 0; j != 4; ++j)
				{
					for (int k = 0; k != 4; ++k)
					{
						int num = i * 16 + j * 4 + k;

						gauss_points[num][0] = gp[i];
						gauss_points[num][1] = gp[j];
						gauss_points[num][2] = gp[k];

						gauss_weights[num] = gw[i] * gw[j] * gw[k];
					}
				}
			}

			break;
		}
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 3D hexaHedron!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature on reference 3D prism element
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Prism integration formulation are constructed by 1D Gauss intergration and 2D triangular integration
	 *        Shape functions and points of integration of the finite element (v11370) (Code_Aster)
	 */
	inline void GaussPrism_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 0.0;
			gauss_points[0][1] = 1.0 / 3.0;
			gauss_points[0][2] = 1.0 / 3.0;

			gauss_weights[0] = 1.0;

			break;
		}
		case 6: //(P3-x,P2-yz)
		{
			double a = 1.0 / sqrt(3.0);
			double c23(2.0 / 3.0), c16(1.0 / 6.0);

			gauss_points[0][0] = -a;
			gauss_points[0][1] = c23;
			gauss_points[0][2] = c16;
			gauss_points[1][0] = -a;
			gauss_points[1][1] = c16;
			gauss_points[1][2] = c23;
			gauss_points[2][0] = -a;
			gauss_points[2][1] = c16;
			gauss_points[2][2] = c16;

			gauss_points[3][0] = a;
			gauss_points[3][1] = c23;
			gauss_points[3][2] = c16;
			gauss_points[4][0] = a;
			gauss_points[4][1] = c16;
			gauss_points[4][2] = c23;
			gauss_points[5][0] = a;
			gauss_points[5][1] = c16;
			gauss_points[5][2] = c16;

			gauss_weights[0] = c16;
			gauss_weights[1] = c16;
			gauss_weights[2] = c16;
			gauss_weights[3] = c16;
			gauss_weights[4] = c16;
			gauss_weights[5] = c16;

			break;
		}
		case 8: // P(3)-Note: this rule has a negative weight and may be unsuitable for some problems.
		{
			double a = 1.0 / sqrt(3.0);
			double c13(1.0 / 3.0), c23(2.0 / 3.0), c16(1.0 / 6.0);
			double w1(-27.0 / 96.0), w2(25.0 / 96.0);

			gauss_points[0][0] = -a;
			gauss_points[0][1] = c23;
			gauss_points[0][2] = c16;
			gauss_points[1][0] = -a;
			gauss_points[1][1] = c16;
			gauss_points[1][2] = c23;
			gauss_points[2][0] = -a;
			gauss_points[2][1] = c16;
			gauss_points[2][2] = c16;
			gauss_points[3][0] = -a;
			gauss_points[3][1] = c13;
			gauss_points[3][2] = c13;

			gauss_points[4][0] = a;
			gauss_points[4][1] = c23;
			gauss_points[4][2] = c16;
			gauss_points[5][0] = a;
			gauss_points[5][1] = c16;
			gauss_points[5][2] = c23;
			gauss_points[6][0] = a;
			gauss_points[6][1] = c16;
			gauss_points[6][2] = c16;
			gauss_points[7][0] = a;
			gauss_points[7][1] = c13;
			gauss_points[7][2] = c13;

			gauss_weights[0] = w2;
			gauss_weights[1] = w2;
			gauss_weights[2] = w2;
			gauss_weights[3] = w1;

			gauss_weights[4] = w2;
			gauss_weights[5] = w2;
			gauss_weights[6] = w2;
			gauss_weights[7] = w1;

			break;
		}
		case 21: //(P5)
		{
			double alpha = sqrt(3.0 / 5.0);
			double a = (6.0 + sqrt(15.0)) / 21.0;
			double b = (6.0 - sqrt(15.0)) / 21.0;
			double c13 = 1.0 / 3.0;
			double c1 = 5.0 / 9.0;
			double c2 = 8.0 / 9.0;
			double w1 = (155.0 + sqrt(15.0)) / 2400.0;
			double w2 = (155.0 - sqrt(15.0)) / 2400.0;

			gauss_points[0][0] = -alpha;
			gauss_points[0][1] = c13;
			gauss_points[0][2] = c13;
			gauss_points[1][0] = -alpha;
			gauss_points[1][1] = a;
			gauss_points[1][2] = a;
			gauss_points[2][0] = -alpha;
			gauss_points[2][1] = 1 - 2 * a;
			gauss_points[2][2] = a;
			gauss_points[3][0] = -alpha;
			gauss_points[3][1] = a;
			gauss_points[3][2] = 1 - 2 * a;
			gauss_points[4][0] = -alpha;
			gauss_points[4][1] = b;
			gauss_points[4][2] = b;
			gauss_points[5][0] = -alpha;
			gauss_points[5][1] = 1 - 2 * b;
			gauss_points[5][2] = b;
			gauss_points[6][0] = -alpha;
			gauss_points[6][1] = b;
			gauss_points[6][2] = 1 - 2 * b;

			gauss_points[7][0] = 0;
			gauss_points[7][1] = c13;
			gauss_points[7][2] = c13;
			gauss_points[8][0] = 0;
			gauss_points[8][1] = a;
			gauss_points[8][2] = a;
			gauss_points[9][0] = 0;
			gauss_points[9][1] = 1 - 2 * a;
			gauss_points[9][2] = a;
			gauss_points[10][0] = 0;
			gauss_points[10][1] = a;
			gauss_points[10][2] = 1 - 2 * a;
			gauss_points[11][0] = 0;
			gauss_points[11][1] = b;
			gauss_points[11][2] = b;
			gauss_points[12][0] = 0;
			gauss_points[12][1] = 1 - 2 * b;
			gauss_points[12][2] = b;
			gauss_points[13][0] = 0;
			gauss_points[13][1] = b;
			gauss_points[13][2] = 1 - 2 * b;

			gauss_points[14][0] = alpha;
			gauss_points[14][1] = c13;
			gauss_points[14][2] = c13;
			gauss_points[15][0] = alpha;
			gauss_points[15][1] = a;
			gauss_points[15][2] = a;
			gauss_points[16][0] = alpha;
			gauss_points[16][1] = 1 - 2 * a;
			gauss_points[16][2] = a;
			gauss_points[17][0] = alpha;
			gauss_points[17][1] = a;
			gauss_points[17][2] = 1 - 2 * a;
			gauss_points[18][0] = alpha;
			gauss_points[18][1] = b;
			gauss_points[18][2] = b;
			gauss_points[19][0] = alpha;
			gauss_points[19][1] = 1 - 2 * b;
			gauss_points[19][2] = b;
			gauss_points[20][0] = alpha;
			gauss_points[20][1] = b;
			gauss_points[20][2] = 1 - 2 * b;

			gauss_weights[0] = c1 * 9.0 / 80.0;
			gauss_weights[1] = c1 * w1;
			gauss_weights[2] = c1 * w1;
			gauss_weights[3] = c1 * w1;
			gauss_weights[4] = c1 * w2;
			gauss_weights[5] = c1 * w2;
			gauss_weights[6] = c1 * w2;

			gauss_weights[7] = c2 * 9.0 / 80.0;
			gauss_weights[8] = c2 * w1;
			gauss_weights[9] = c2 * w1;
			gauss_weights[10] = c2 * w1;
			gauss_weights[11] = c2 * w2;
			gauss_weights[12] = c2 * w2;
			gauss_weights[13] = c2 * w2;

			gauss_weights[14] = c1 * 9.0 / 80.0;
			gauss_weights[15] = c1 * w1;
			gauss_weights[16] = c1 * w1;
			gauss_weights[17] = c1 * w1;
			gauss_weights[18] = c1 * w2;
			gauss_weights[19] = c1 * w2;
			gauss_weights[20] = c1 * w2;

			break;
		}
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 3D prism!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}

	/**
	 * Get points and weights of Gauss quadrature on reference 3D pyramid element
	 * \param[in] num_point : num of gauss points
	 * \param[in|out] gauss_points : position of each gauss point
	 * \param[in|out] gauss_weights : weight of each gauss point
	 * \note: Pyramid rule: http://people.sc.fsu.edu/~jburkardt/m_src/pyramid_rule/pyramid_rule.html
	 */
	inline void GaussPyramid_ref(int num_point, Array2D<double> &gauss_points, Array1D<double> &gauss_weights)
	{
		switch (num_point)
		{
		case 1: //(P1)
		{
			gauss_points[0][0] = 0.0;
			gauss_points[0][1] = 0.0;
			gauss_points[0][2] = 0.25;

			gauss_weights[0] = 4.0 / 3.0;

			break;
		}
		case 8: //(P3)
		{
			double w1 = 0.1744105884401307 * 4.0 / 3.0;
			double w2 = 0.0755894115598690 * 4.0 / 3.0;

			gauss_points[0][0] = -0.5066163033497875;
			gauss_points[0][1] = -0.5066163033497875;
			gauss_points[0][2] = 0.1225148226554414;
			gauss_points[1][0] = 0.5066163033497875;
			gauss_points[1][1] = -0.5066163033497875;
			gauss_points[1][2] = 0.1225148226554414;
			gauss_points[2][0] = -0.5066163033497875;
			gauss_points[2][1] = 0.5066163033497875;
			gauss_points[2][2] = 0.1225148226554414;
			gauss_points[3][0] = 0.5066163033497875;
			gauss_points[3][1] = 0.5066163033497875;
			gauss_points[3][2] = 0.1225148226554414;

			gauss_points[4][0] = -0.2631840555697136;
			gauss_points[4][1] = -0.2631840555697136;
			gauss_points[4][2] = 0.5441518440112253;
			gauss_points[5][0] = 0.2631840555697136;
			gauss_points[5][1] = -0.2631840555697136;
			gauss_points[5][2] = 0.5441518440112253;
			gauss_points[6][0] = -0.2631840555697136;
			gauss_points[6][1] = 0.2631840555697136;
			gauss_points[6][2] = 0.5441518440112253;
			gauss_points[7][0] = 0.2631840555697136;
			gauss_points[7][1] = 0.2631840555697136;
			gauss_points[7][2] = 0.5441518440112253;

			gauss_weights[0] = w1;
			gauss_weights[1] = w1;
			gauss_weights[2] = w1;
			gauss_weights[3] = w1;
			gauss_weights[4] = w2;
			gauss_weights[5] = w2;
			gauss_weights[6] = w2;
			gauss_weights[7] = w2;

			break;
		}
		case 27: //(P5)
		{
			double a(0), b(0);
			double w1(0), w2(0), w3(0);

			a = 0.7180557413198889;
			b = 0.0729940240731498;

			gauss_points[0][0] = -a;
			gauss_points[0][1] = -a;
			gauss_points[0][2] = b;
			gauss_points[1][0] = 0;
			gauss_points[1][1] = -a;
			gauss_points[1][2] = b;
			gauss_points[2][0] = a;
			gauss_points[2][1] = -a;
			gauss_points[2][2] = b;
			gauss_points[3][0] = -a;
			gauss_points[3][1] = 0;
			gauss_points[3][2] = b;
			gauss_points[4][0] = 0;
			gauss_points[4][1] = 0;
			gauss_points[4][2] = b;
			gauss_points[5][0] = a;
			gauss_points[5][1] = 0;
			gauss_points[5][2] = b;
			gauss_points[6][0] = -a;
			gauss_points[6][1] = a;
			gauss_points[6][2] = b;
			gauss_points[7][0] = 0;
			gauss_points[7][1] = a;
			gauss_points[7][2] = b;
			gauss_points[8][0] = a;
			gauss_points[8][1] = a;
			gauss_points[8][2] = b;

			a = 0.5058087078539250;
			b = 0.3470037660383519;

			gauss_points[9][0] = -a;
			gauss_points[9][1] = -a;
			gauss_points[9][2] = b;
			gauss_points[10][0] = 0;
			gauss_points[10][1] = -a;
			gauss_points[10][2] = b;
			gauss_points[11][0] = a;
			gauss_points[11][1] = -a;
			gauss_points[11][2] = b;
			gauss_points[12][0] = -a;
			gauss_points[12][1] = 0;
			gauss_points[12][2] = b;
			gauss_points[13][0] = 0;
			gauss_points[13][1] = 0;
			gauss_points[13][2] = b;
			gauss_points[14][0] = a;
			gauss_points[14][1] = 0;
			gauss_points[14][2] = b;
			gauss_points[15][0] = -a;
			gauss_points[15][1] = a;
			gauss_points[15][2] = b;
			gauss_points[16][0] = 0;
			gauss_points[16][1] = a;
			gauss_points[16][2] = b;
			gauss_points[17][0] = a;
			gauss_points[17][1] = a;
			gauss_points[17][2] = b;

			a = 0.2285043056539673;
			b = 0.7050022098884984;

			gauss_points[18][0] = -a;
			gauss_points[18][1] = -a;
			gauss_points[18][2] = b;
			gauss_points[19][0] = 0;
			gauss_points[19][1] = -a;
			gauss_points[19][2] = b;
			gauss_points[20][0] = a;
			gauss_points[20][1] = -a;
			gauss_points[20][2] = b;
			gauss_points[21][0] = -a;
			gauss_points[21][1] = 0;
			gauss_points[21][2] = b;
			gauss_points[22][0] = 0;
			gauss_points[22][1] = 0;
			gauss_points[22][2] = b;
			gauss_points[23][0] = a;
			gauss_points[23][1] = 0;
			gauss_points[23][2] = b;
			gauss_points[24][0] = -a;
			gauss_points[24][1] = a;
			gauss_points[24][2] = b;
			gauss_points[25][0] = 0;
			gauss_points[25][1] = a;
			gauss_points[25][2] = b;
			gauss_points[26][0] = a;
			gauss_points[26][1] = a;
			gauss_points[26][2] = b;

			w1 = 0.0363741576539090 * 4.0 / 3.0;
			w2 = 0.0581986522462544 * 4.0 / 3.0;
			w3 = 0.0931178435940069 * 4.0 / 3.0;

			gauss_weights[0] = w1;
			gauss_weights[1] = w2;
			gauss_weights[2] = w1;
			gauss_weights[3] = w2;
			gauss_weights[4] = w3;
			gauss_weights[5] = w2;
			gauss_weights[6] = w1;
			gauss_weights[7] = w2;
			gauss_weights[8] = w1;

			w1 = 0.0338533030694135 * 4.0 / 3.0;
			w2 = 0.0541652849110615 * 4.0 / 3.0;
			w3 = 0.0866644558576984 * 4.0 / 3.0;

			gauss_weights[9] = w1;
			gauss_weights[10] = w2;
			gauss_weights[11] = w1;
			gauss_weights[12] = w2;
			gauss_weights[13] = w3;
			gauss_weights[14] = w2;
			gauss_weights[15] = w1;
			gauss_weights[16] = w2;
			gauss_weights[17] = w1;

			w1 = 0.0069330331038381 * 4.0 / 3.0;
			w2 = 0.0110928529661410 * 4.0 / 3.0;
			w3 = 0.0177485647458256 * 4.0 / 3.0;

			gauss_weights[18] = w1;
			gauss_weights[19] = w2;
			gauss_weights[20] = w1;
			gauss_weights[21] = w2;
			gauss_weights[22] = w3;
			gauss_weights[23] = w2;
			gauss_weights[24] = w1;
			gauss_weights[25] = w2;
			gauss_weights[26] = w1;

			break;
		}
		default:
		{
			std::cout << "Error: Gauss quadrature rule is not supported for 3D pyramid!" << std::endl;
			std::cout << "File: " << __FILE__ << std::endl;
			std::cout << "Line: " << __LINE__ << std::endl;
			exit(1);
		}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
#endif
