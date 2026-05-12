#ifndef __EQUATIONS__HPP__
#define __EQUATIONS__HPP__

#include "../third_party/scTools/scTools.h"

enum RIEMANNFLUXTYPE
{
	LLF = 0,
	HLLC = 1
};
enum SCHEMETYPE
{
	WENO = 0,
	WENOZ = 1,
};

//	方程的性质
enum INITIALVALUETYPE
{
	CONSTANT = 0,
	POLYNOMIAL = 1
};
enum TESTCASE
{
	SMOOTH = 0,
	Lax_problem = 1,
	Sod_problem = 2,
	SHU_OSHER_PROBLEM = 3,
	BLAST_WAVE = 4
};
enum BOUNDARYTYPE
{
	PERIOD = 0,
	NEUMANN = 1,
	REFLECTIVE = 2
};

class BaseEquation
{
public:
	double xL, xR, outputime;
	BOUNDARYTYPE boundaryCondition;
	INITIALVALUETYPE initialCondition;
	bool u_exact_exist = false;
	std::function<double(double, double)> theVarExact;
	std::function<double(Array1D<double>)> theVarUh;

	virtual ~BaseEquation() {}
	virtual int getVarNum() = 0;
	virtual void getPhyFlux(const Array1D<double> Uh, Array1D<double> &Flux) = 0;
	virtual double getMaxEigenValue(const Array1D<double> Uh) = 0;
	virtual void getLEigenMatrix(const Array1D<double> Uh, Array2D<double> &eigMatrix) = 0;
	virtual void getREigenMatrix(const Array1D<double> Uh, Array2D<double> &eigMatrix) = 0;
	virtual void getLLFRiemannFlux(const Array1D<double> UL, const Array1D<double> UR, Array1D<double> &Flux) = 0;
	virtual void setEquationParameters(TESTCASE type) = 0;
	virtual void getU0(double xP, Array1D<double> &U) = 0;
	virtual int getVitalVarNum() = 0;
	virtual void getVitalVarName(Array1D<std::string> &VitalVarName) = 0;
	virtual void getVitalVarVal(const Array1D<double> Uh, Array1D<double> &VitalVar) = 0;
};

class EulerEquation : public BaseEquation
{
public:
	// ------------- (1) 方程形式与特征值，特征矩阵 -------------//
	// (1-0) 变量
	double gamma;

	// (1-1) 方程数量
	inline int getVarNum() override
	{
		return 3;
	}

	// (1-2) PhyFlux
	inline void getPhyFlux(const Array1D<double> Uh, Array1D<double> &Flux) override
	{
		double rho(0), u(0), E(0), pre(0);

		// U
		rho = Uh[0];
		u = Uh[1] / Uh[0];
		E = Uh[2];
		pre = (gamma - 1.0) * (E - 0.5 * rho * u * u);

		// Flux
		Flux[0] = rho * u;
		Flux[1] = rho * u * u + pre;
		Flux[2] = u * (E + pre);
	}

	// (1-3) eigen Values
	inline double getMaxEigenValue(const Array1D<double> Uh) override
	{
		double rho(0), pre(0), c(0), u(0), E(0);
		rho = Uh[0];
		u = Uh[1] / Uh[0];
		E = Uh[2];
		pre = (gamma - 1.0) * (E - 0.5 * rho * u * u);
		c = sqrt(gamma * pre / rho);

		return fabs(u) + c;
	}

	// (1-4) Left Eigen Matrix
	inline void getLEigenMatrix(const Array1D<double> Uh, Array2D<double> &eigMatrix) override
	{
		double rho(0), u(0), E(0), pre(0), c(0), B1(0), B2(0);

		// U
		rho = Uh[0];
		u = Uh[1] / Uh[0];
		E = Uh[2];
		pre = (gamma - 1.0) * (E - 0.5 * rho * u * u);
		c = sqrt(gamma * pre / rho);
		B1 = (gamma - 1) / c / c;
		B2 = B1 * u * u / 2;

		eigMatrix[0][0] = 0.5 * (B2 + u / c);
		eigMatrix[1][0] = 1 - B2;
		eigMatrix[2][0] = 0.5 * (B2 - u / c);

		eigMatrix[0][1] = -(B1 * u + 1.0 / c) / 2;
		eigMatrix[1][1] = B1 * u;
		eigMatrix[2][1] = -(B1 * u - 1.0 / c) / 2;

		eigMatrix[0][2] = 0.5 * B1;
		eigMatrix[1][2] = -B1;
		eigMatrix[2][2] = 0.5 * B1;
	}

	// (1-5) Right Eigen Matrix
	inline void getREigenMatrix(const Array1D<double> Uh, Array2D<double> &eigMatrix) override
	{
		double rho(0), u(0), E(0), pre(0), c(0), H(0);

		// U
		rho = Uh[0];
		u = Uh[1] / Uh[0];
		E = Uh[2];
		pre = (gamma - 1.0) * (E - 0.5 * rho * u * u);
		c = sqrt(gamma * pre / rho);
		H = (E + pre) / rho;

		eigMatrix[0][0] = 1.0;
		eigMatrix[1][0] = u - c;
		eigMatrix[2][0] = H - u * c;

		eigMatrix[0][1] = 1.0;
		eigMatrix[1][1] = u;
		eigMatrix[2][1] = 0.5 * u * u;

		eigMatrix[0][2] = 1.0;
		eigMatrix[1][2] = u + c;
		eigMatrix[2][2] = H + u * c;
	}

	// (1-6) Num flux - (a) LLF
	inline void getLLFRiemannFlux(const Array1D<double> UL, const Array1D<double> UR, Array1D<double> &Flux) override
	{
		// LLF approximate Riemann flux
		int varNum = getVarNum();
		double ws(0), ws_L(0), ws_R(0);
		Array1D<double> FUL(varNum, 0.0), FUR(varNum, 0.0);

		getPhyFlux(UL, FUL);
		getPhyFlux(UR, FUR);
		ws_L = getMaxEigenValue(UL);
		ws_R = getMaxEigenValue(UR);
		ws = max(ws_L, ws_R);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
			Flux[r] = 0.5 * (FUL[r] + FUR[r] - ws * (UR[r] - UL[r]));
	}

	// -------- (2) 方程的初边值条件与参数设置 ----------- //
	// (2-1) 需要用到的物理变量
	std::function<double(double)> u0, pre0, rho0;

	// (2-2) 根据算例进行方程设置
	void setEquationParameters(TESTCASE type) override
	{
		switch (type)
		{
		case SMOOTH:
			// 参数设置
			this->xL = 0.0;
			this->xR = 1.0;
			this->outputime = 2;
			this->boundaryCondition = PERIOD;
			this->initialCondition = POLYNOMIAL;
			gamma = 1.4;

			// 初边值条件
			rho0 = [](double x)
			{ return 1 + 0.2 * sin(2 * M_PI * x); };
			u0 = [](double x)
			{ return 1; };
			pre0 = [](double x)
			{ return 1; };

			// 函数解析解
			this->u_exact_exist = true;
			this->theVarExact = [](double x, double t)
			{ return 1 + 0.2 * sin(2 * M_PI * (x - t)); };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;
		case Lax_problem:
			// 参数设置
			this->xL = 0.0;
			this->xR = 1.0;
			this->outputime = 0.14;
			this->boundaryCondition = NEUMANN;
			this->initialCondition = CONSTANT;
			gamma = 1.4;

			// 初边值条件
			rho0 = [](double x)
			{ return (x < 0.5) ? 0.445 : 0.5; };
			u0 = [](double x)
			{ return (x < 0.5) ? 0.698 : 0; };
			pre0 = [](double x)
			{ return (x < 0.5) ? 3.528 : 0.571; };

			// 函数解析解
			this->u_exact_exist = false;
			this->theVarExact = [](double x, double t)
			{
				std::cout << "该算例无解析解，请检查。" << std::endl;
				std::cin.get();
				exit(1);
				return 0.0;
			};
			break;
		case Sod_problem:
			// 参数设置
			this->xL = 0.0;
			this->xR = 1.0;
			this->outputime = 0.2;
			this->boundaryCondition = NEUMANN;
			this->initialCondition = CONSTANT;
			gamma = 1.4;

			// 初边值条件
			rho0 = [](double x)
			{ return (x < 0.5) ? 1 : 0.125; };
			u0 = [](double x)
			{ return 0; };
			pre0 = [](double x)
			{ return (x < 0.5) ? 1 : 0.1; };

			// 函数解析解
			this->u_exact_exist = false;
			this->theVarExact = [](double x, double t)
			{
				std::cout << "该算例无解析解，请检查。" << std::endl;
				std::cin.get();
				exit(1);
				return 0.0;
			};
			break;
		case SHU_OSHER_PROBLEM:
			// 参数设置
			xL = -5.0;
			xR = 5.0;
			outputime = 1.8;
			boundaryCondition = NEUMANN;
			initialCondition = CONSTANT;
			gamma = 1.4;

			// 初边值条件
			rho0 = [](double x)
			{ return (x < -4) ? 3.857143 : 1 + 0.2 * sin(5 * x); };
			u0 = [](double x)
			{ return (x < -4) ? 2.629369 : 0; };
			pre0 = [](double x)
			{ return (x < -4) ? 10.333333 : 1; };

			// 函数解析解
			this->u_exact_exist = false;
			this->theVarExact = [](double x, double t)
			{
				std::cout << "该算例无解析解，请检查。" << std::endl;
				std::cin.get();
				exit(1);
				return 0.0;
			};
			break;
		case BLAST_WAVE:
			this->xL = 0.0;
			this->xR = 1.0;
			this->outputime = 0.038;
			this->boundaryCondition = REFLECTIVE;
			this->initialCondition = CONSTANT;
			gamma = 1.4;

			rho0 = [](double x)
			{ return 1.0; };
			u0 = [](double x)
			{ return 0.0; };
			pre0 = [](double x)
			{
				if (x < 0.1) return 1000.0;
				else if (x < 0.9) return 0.01;
				else return 100.0;
			};

			this->u_exact_exist = false;
			this->theVarExact = [](double x, double t)
			{
				std::cout << "该算例无解析解，请检查。" << std::endl;
				std::cin.get();
				exit(1);
				return 0.0;
			};
			break;
		default:
			throw std::invalid_argument("Invalid example number");
		}
	}

	// (2-3) 将物理变量转换成双曲型方程 Ut+f(U)x = 0 的 U
	inline void getU0(double xP, Array1D<double> &U) override
	{
		double rho = rho0(xP);
		double u = u0(xP);
		double pre = pre0(xP);

		double E = pre / (gamma - 1.0) + 0.5 * rho * u * u;

		// Flux
		U[0] = rho;
		U[1] = rho * u;
		U[2] = E;
	}

	// ---------- (3) 输出 --------------//
	// (3-1) 输出变量的数量
	inline int getVitalVarNum() override
	{
		return 5;
	}

	// (3-2) 输出变量的名字
	inline void getVitalVarName(Array1D<std::string> &VitalVarName) override
	{
		VitalVarName[0] = "rho";
		VitalVarName[1] = "u";
		VitalVarName[2] = "pre";
		VitalVarName[3] = "c";
		VitalVarName[4] = "E";
	}

	// (3-3) 输出变量的值
	inline void getVitalVarVal(const Array1D<double> Uh, Array1D<double> &VitalVar) override
	{
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double E = Uh[2];
		double pre = (gamma - 1.0) * (E - 0.5 * rho * u * u);
		double c = sqrt(gamma * pre / rho);

		VitalVar[0] = rho;
		VitalVar[1] = u;
		VitalVar[2] = pre;
		VitalVar[3] = c;
		VitalVar[4] = E;
	}
};
#endif