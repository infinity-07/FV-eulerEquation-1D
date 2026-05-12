#ifndef __WENOFV__HPP__
#define __WENOFV__HPP__

#include <iostream>
#include <fstream>
#include "../third_party/scTools/scTools.h"

#include "equations.hpp"
#include <map>
#include <random>

// WENO reconstruction types
enum RECONSTRUCTIONTYPE
{
    CONVERVATIVE = 0,
    CHARACTERISTIC = 1,
};

class Cgrid
{
public:
    // 成员变量
    double m_xL;
    double m_xR;
    double m_xC;
    double m_xD;
};

class minusAndPlusFlux
{
public:
    Array1D<double> flux_minus;
    Array1D<double> flux_plus;
};

class CWENOFV
{
public:
    CWENOFV() {}
    ~CWENOFV();

    CWENOFV(EulerEquation *equation)
    {
        this->equation = equation;
    }

    EulerEquation *equation;

public:
    // Variables
    int m_elemNum;
    int m_varNum;
    double m_globalXL;
    double m_globalXR;
    double m_cfl;
    double m_deltaX;
    double m_outputTime;
    double m_restartTime;
    double m_now;
    int m_count;

    int m_ghostNum, m_startElemX, m_endElemX;

    // Arrays
    Array2D<double> m_Uh;
    Array1D<Cgrid> m_grids;
    Array1D<minusAndPlusFlux> m_flux;
    Array1D<Array2D<double>> m_cellFlux;
    Array1D<double> m_troubledCell;
    Array1D<double> m_badcell;
Array1D<double> m_lambdaVec;

    RECONSTRUCTIONTYPE m_reconstructionType;

    // Timer
    CTime m_timer;

    // Scheme and Testcase
    SCHEMETYPE m_scheme;
    TESTCASE m_testcase;

    // Methods
    void initializeSolver(std::map<std::string, std::string> option);

    void initializeAve();

    void run();
    void RunRK();

    double calculateDeltaT();
    void assembleRHS(Array2D<double> &L);
    void getFlux_conservative();
    void getFlux_characteristic();

    void getL(Array2D<double> &L);

    // Output
    void outputAve(std::string prefix);
    void outputCharacteristicAve(std::string prefix);
    void outputAccuracy(string prefix);

};

///////////////////////////////////////////////////////////////////////
#endif