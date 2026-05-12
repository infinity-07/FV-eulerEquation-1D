#include "weno_fv.hpp"
#include <cmath>
#include "equations.hpp"
#include "weno_reconstructor.hpp"
#include <iomanip>
#include <map>
#include <filesystem>
#include "MyTimer.hpp"

CWENOFV::~CWENOFV()
{
}

void CWENOFV::initializeSolver(std::map<std::string, std::string> option)
{
    m_testcase = (TESTCASE)std::stoi(option["TESTCASE"]);
    m_scheme = (SCHEMETYPE)std::stoi(option["SCHEME"]);
    m_cfl = std::stod(option["CFL"]);
    m_elemNum = std::stoi(option["ELEMNUM"]);

    equation->setEquationParameters(m_testcase);

    m_ghostNum = 3;

    /////////////////////////////////////////////////
    std::cout << "Initializing variables" << std::endl;
    m_varNum = equation->getVarNum();
    m_outputTime = equation->outputime;
    m_globalXL = equation->xL;
    m_globalXR = equation->xR;
    m_deltaX = (m_globalXR - m_globalXL) / m_elemNum;

    // Generate grid
    std::cout << "Generating grid..." << std::endl;
    m_grids.Resize(m_elemNum + 2 * m_ghostNum);
    for (int e = 0; e < m_elemNum + 2 * m_ghostNum; ++e)
    {
        m_grids[e].m_xL = m_globalXL + e * m_deltaX - m_ghostNum * m_deltaX;
        m_grids[e].m_xR = m_globalXL + (e + 1) * m_deltaX - m_ghostNum * m_deltaX;

        m_grids[e].m_xC = (m_grids[e].m_xL + m_grids[e].m_xR) / 2.0;
        m_grids[e].m_xD = m_grids[e].m_xR - m_grids[e].m_xL;
    }

    m_reconstructionType = (RECONSTRUCTIONTYPE)std::stoi(option["RECONSTRUCTIONTYPE"]);
    m_usePositivityLimiter = option.count("POSITIVITY_LIMITER") && std::stoi(option["POSITIVITY_LIMITER"]) != 0;

    // Initialize cell flux arrays
    m_cellFlux.Resize(m_elemNum + 2 * m_ghostNum);
    for (int ei = 0; ei != m_elemNum + 2 * m_ghostNum; ei++)
        m_cellFlux[ei].Resize(2, m_varNum);

    m_Uh.Resize(m_elemNum + 2 * m_ghostNum, m_varNum);

    m_startElemX = m_ghostNum;
    m_endElemX = m_elemNum + m_ghostNum;
}

void CWENOFV::initializeAve()
{
    const int gpNum = 4;
    Array1D<double> gpoints_ref(gpNum), gweights_ref(gpNum);
    sc_math::GaussLegendre_ref(gpNum, gpoints_ref, gweights_ref);

    std::cout << "Initialize cell averages..." << std::endl;
    double gPoint, gWeight;
    m_Uh.setZero();

    for (int e = m_startElemX; e != m_endElemX; ++e)
    {
        for (int gp = 0; gp != gpNum; ++gp)
        {
            gPoint = m_grids[e].m_xC + 0.5 * m_grids[e].m_xD * gpoints_ref[gp];
            gWeight = 0.5 * m_grids[e].m_xD * gweights_ref[gp];

            Array1D<double> Conserved_var(m_varNum);
            equation->getU0(gPoint, Conserved_var);

            for (int r = 0; r != m_varNum; ++r)
                m_Uh[e][r] += Conserved_var[r] * gWeight;
        }
        for (int r = 0; r != m_varNum; ++r)
            m_Uh[e][r] = m_Uh[e][r] / m_grids[e].m_xD;
    }
}

void CWENOFV::run()
{
    m_now = 0;
    m_count = 0;
    Timer progressTimer;
    std::cout << "Start iteration..." << std::endl;
    progressTimer.start();
    while (fabs(m_now - m_outputTime) > 1e-9)
    {
        RunRK();
        progressTimer.updateProgress(m_now / equation->outputime);
    }
    double executionTime = progressTimer.stop();
    std::cout << "\nTotal execution time: " << executionTime << " seconds\n";
    std::cout << "\nIteration complete" << std::endl;
}

void CWENOFV::RunRK()
{
    // Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001).
    // Strong Stability-Preserving High-Order Time Discretization Methods.
    // SIAM Review, 43(1), 89-112. https://doi.org/10.1137/s003614450036757x

    Array2D<double> L(m_elemNum + 2 * m_ghostNum, m_varNum);
    Array2D<double> m_Un(m_elemNum + 2 * m_ghostNum, m_varNum);

    const double deltaT = calculateDeltaT();

    for (int e = m_startElemX; e != m_endElemX; ++e)
        for (int r = 0; r != m_varNum; ++r)
            m_Un[e][r] = m_Uh[e][r];

    for (int step = 0; step != 3; step++)
    {
        assembleRHS(L);

        double a, b, c;
        switch (step)
        {
        case 0:
            a = 1.0;  b = 0.0;        c = 1.0;        break;
        case 1:
            a = 3.0 / 4.0; b = 1.0 / 4.0; c = 1.0 / 4.0; break;
        case 2:
            a = 1.0 / 3.0; b = 2.0 / 3.0; c = 2.0 / 3.0; break;
        default:
            std::cout << "Error: Invalid RK3 step" << std::endl;
            std::cin.get();
            exit(1);
        }

        for (int e = m_startElemX; e != m_endElemX; ++e)
            for (int r = 0; r != m_varNum; ++r)
                m_Uh[e][r] = a * m_Un[e][r] + b * m_Uh[e][r] + c * deltaT * L[e][r];
    }

    m_now = m_now + deltaT;
    m_count++;
}

void CWENOFV::assembleRHS(Array2D<double> &L)
{
    switch (m_reconstructionType)
    {
    case CONVERVATIVE:
        getFlux_conservative();
        break;
    case CHARACTERISTIC:
        getFlux_characteristic();
        break;
    default:
        std::cout << "Unsupported reconstruction type for EulerEquation" << std::endl;
        std::cin.get();
        break;
    }

    if (m_usePositivityLimiter)
        applyPositivityLimiter();

    getL(L);
}

void CWENOFV::getL(Array2D<double> &L)
{
    Array1D<double> flux_right(m_varNum);
    Array1D<double> flux_left(m_varNum);
    Array1D<double> u_right_minus(m_varNum);
    Array1D<double> u_right_plus(m_varNum);
    Array1D<double> u_left_minus(m_varNum);
    Array1D<double> u_left_plus(m_varNum);

    if (equation->boundaryCondition == PERIOD)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int gp = 0; gp != 2; gp++)
            {
                m_cellFlux[m_startElemX - 1][gp][r] = m_cellFlux[m_endElemX - 1][gp][r];
                m_cellFlux[m_endElemX][gp][r] = m_cellFlux[m_startElemX][gp][r];
            }
    }
    else if (equation->boundaryCondition == NEUMANN)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int gp = 0; gp != 2; gp++)
            {
                m_cellFlux[m_startElemX - 1][gp][r] = m_cellFlux[m_startElemX][gp][r];
                m_cellFlux[m_endElemX][gp][r] = m_cellFlux[m_endElemX - 1][gp][r];
            }
    }
    else // REFLECTIVE
    {
        for (int r = 0; r < m_varNum; r++)
        {
            const double sign = (r == 1) ? -1.0 : 1.0;
            for (int gp = 0; gp != 2; gp++)
            {
                m_cellFlux[m_startElemX - 1][gp][r] = sign * m_cellFlux[m_startElemX][1 - gp][r];
                m_cellFlux[m_endElemX][gp][r]        = sign * m_cellFlux[m_endElemX - 1][1 - gp][r];
            }
        }
    }

    for (int e = m_startElemX; e < m_endElemX; e++)
    {
        for (int r = 0; r != m_varNum; r++)
        {
            u_right_minus[r] = m_cellFlux[e][1][r];
            u_right_plus[r]  = m_cellFlux[e + 1][0][r];
            u_left_minus[r]  = m_cellFlux[e - 1][1][r];
            u_left_plus[r]   = m_cellFlux[e][0][r];
        }

        equation->getLLFRiemannFlux(u_left_minus, u_left_plus, flux_left);
        equation->getLLFRiemannFlux(u_right_minus, u_right_plus, flux_right);

        for (int r = 0; r < m_varNum; r++)
            L[e][r] = -(flux_right[r] - flux_left[r]) / m_grids[e].m_xD;
    }
}

void CWENOFV::getFlux_conservative()
{
    Array1D<double> uavemm(m_varNum), uavem(m_varNum), uave(m_varNum), uavep(m_varNum), uavepp(m_varNum);
    Array1D<double> u_minus(m_varNum), u_plus(m_varNum);
    Array2D<double> tempAve(7, m_varNum);

    if (equation->boundaryCondition == PERIOD)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[e][r] = m_Uh[m_endElemX - m_ghostNum + e][r];
                m_Uh[m_endElemX + e][r] = m_Uh[m_ghostNum + e][r];
            }
    }
    else if (equation->boundaryCondition == NEUMANN)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[e][r] = m_Uh[m_startElemX][r];
                m_Uh[m_endElemX + e][r] = m_Uh[m_endElemX - 1][r];
            }
    }
    else // REFLECTIVE
    {
        for (int r = 0; r < m_varNum; r++)
        {
            const double sign = (r == 1) ? -1.0 : 1.0;
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[m_startElemX - 1 - e][r] = sign * m_Uh[m_startElemX + e][r];
                m_Uh[m_endElemX + e][r]        = sign * m_Uh[m_endElemX - 1 - e][r];
            }
        }
    }

    for (int ei = m_startElemX; ei < m_endElemX; ei++)
    {
        for (int tempIndex = 0; tempIndex < 7; tempIndex++)
            for (int r = 0; r < m_varNum; r++)
                tempAve[tempIndex][r] = m_Uh[ei - 3 + tempIndex][r];

        for (int r = 0; r != m_varNum; r++)
        {
            uavemm[r] = tempAve[1][r];
            uavem[r]  = tempAve[2][r];
            uave[r]   = tempAve[3][r];
            uavep[r]  = tempAve[4][r];
            uavepp[r] = tempAve[5][r];
        }

        for (int r = 0; r < m_varNum; r++)
        {
            switch (m_scheme)
            {
            case WENO:
                u_minus[r] = weno(uavemm[r], uavem[r], uave[r], uavep[r], uavepp[r]);
                break;
            case WENOZ:
                u_minus[r] = wenoz(uavemm[r], uavem[r], uave[r], uavep[r], uavepp[r]);
                break;
            default:
                std::cout << "没有可执行方案" << std::endl;
                std::cin.get();
                break;
            }
        }

        for (int r = 0; r < m_varNum; r++)
        {
            switch (m_scheme)
            {
            case WENO:
                u_plus[r] = weno(uavepp[r], uavep[r], uave[r], uavem[r], uavemm[r]);
                break;
            case WENOZ:
                u_plus[r] = wenoz(uavepp[r], uavep[r], uave[r], uavem[r], uavemm[r]);
                break;
            default:
                std::cout << "没有可执行方案" << std::endl;
                std::cin.get();
                break;
            }
        }

        for (int r = 0; r != m_varNum; r++)
        {
            m_cellFlux[ei][0][r] = u_plus[r];
            m_cellFlux[ei][1][r] = u_minus[r];
        }
    }
}

void CWENOFV::getFlux_characteristic()
{
    Array1D<double> uavemm(m_varNum), uavem(m_varNum), uave(m_varNum), uavep(m_varNum), uavepp(m_varNum);
    Array1D<double> uavemm_c(m_varNum), uavem_c(m_varNum), uave_c(m_varNum), uavep_c(m_varNum), uavepp_c(m_varNum);
    Array1D<double> u_minus(m_varNum), u_plus(m_varNum);
    Array2D<double> eigMatrixL(m_varNum, m_varNum), eigMatrixR(m_varNum, m_varNum);
    Array2D<double> tempAve(7, m_varNum);

    if (equation->boundaryCondition == PERIOD)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[e][r] = m_Uh[m_endElemX - m_ghostNum + e][r];
                m_Uh[m_endElemX + e][r] = m_Uh[m_ghostNum + e][r];
            }
    }
    else if (equation->boundaryCondition == NEUMANN)
    {
        for (int r = 0; r < m_varNum; r++)
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[e][r] = m_Uh[m_startElemX][r];
                m_Uh[m_endElemX + e][r] = m_Uh[m_endElemX - 1][r];
            }
    }
    else // REFLECTIVE
    {
        for (int r = 0; r < m_varNum; r++)
        {
            const double sign = (r == 1) ? -1.0 : 1.0;
            for (int e = 0; e < m_ghostNum; e++)
            {
                m_Uh[m_startElemX - 1 - e][r] = sign * m_Uh[m_startElemX + e][r];
                m_Uh[m_endElemX + e][r]        = sign * m_Uh[m_endElemX - 1 - e][r];
            }
        }
    }

    for (int ei = m_startElemX; ei < m_endElemX; ei++)
    {
        for (int tempIndex = 0; tempIndex < 7; tempIndex++)
            for (int r = 0; r < m_varNum; r++)
                tempAve[tempIndex][r] = m_Uh[ei - 3 + tempIndex][r];

        for (int r = 0; r != m_varNum; r++)
        {
            uavemm[r] = tempAve[1][r];
            uavem[r]  = tempAve[2][r];
            uave[r]   = tempAve[3][r];
            uavep[r]  = tempAve[4][r];
            uavepp[r] = tempAve[5][r];
        }

        equation->getLEigenMatrix(uave, eigMatrixL);
        equation->getREigenMatrix(uave, eigMatrixR);
        uavemm_c = eigMatrixL * uavemm;
        uavem_c  = eigMatrixL * uavem;
        uave_c   = eigMatrixL * uave;
        uavep_c  = eigMatrixL * uavep;
        uavepp_c = eigMatrixL * uavepp;

        for (int r = 0; r < m_varNum; r++)
        {
            switch (m_scheme)
            {
            case WENO:
                u_minus[r] = weno(uavemm_c[r], uavem_c[r], uave_c[r], uavep_c[r], uavepp_c[r]);
                break;
            case WENOZ:
                u_minus[r] = wenoz(uavemm_c[r], uavem_c[r], uave_c[r], uavep_c[r], uavepp_c[r]);
                break;
            default:
                std::cout << "没有可执行方案" << std::endl;
                std::cin.get();
                break;
            }
        }

        for (int r = 0; r < m_varNum; r++)
        {
            switch (m_scheme)
            {
            case WENO:
                u_plus[r] = weno(uavepp_c[r], uavep_c[r], uave_c[r], uavem_c[r], uavemm_c[r]);
                break;
            case WENOZ:
                u_plus[r] = wenoz(uavepp_c[r], uavep_c[r], uave_c[r], uavem_c[r], uavemm_c[r]);
                break;
            default:
                std::cout << "没有可执行方案" << std::endl;
                std::cin.get();
                break;
            }
        }

        u_minus = eigMatrixR * u_minus;
        u_plus  = eigMatrixR * u_plus;

        for (int r = 0; r != m_varNum; r++)
        {
            m_cellFlux[ei][0][r] = u_plus[r];
            m_cellFlux[ei][1][r] = u_minus[r];
        }
    }
}

double CWENOFV::calculateDeltaT()
{
    double timestep(1.0), tmp(0);
    Array1D<double> Uave(m_varNum);
    double eigenvalue(0);

    Uave.setZero();
    for (int e = m_startElemX; e != m_endElemX; ++e)
    {
        for (int r = 0; r != m_varNum; ++r)
            Uave[r] = m_Uh[e][r];

        eigenvalue = equation->getMaxEigenValue(Uave);
        tmp = m_cfl * m_deltaX / eigenvalue;

        if (tmp > timestep)
            continue;
        else
            timestep = tmp;
    }

    if (m_outputTime - m_now < timestep)
        timestep = m_outputTime - m_now;

    if (timestep != timestep)
    {
        std::cout << "Error: Timestep is not a number..." << std::endl;
        std::cin.get();
        exit(1);
    }

    return timestep;
}

void CWENOFV::applyPositivityLimiter()
{
    // Zhang, X., & Shu, C.-W. (2011).
    // Maximum-principle-satisfying and positivity-preserving high-order schemes
    // for conservation laws: survey and new developments.
    // Proceedings of the Royal Society A, 467(2134), 2752-2776.

    const double eps   = 1e-4;
    const double gamma = equation->gamma;

    for (int ei = m_startElemX; ei < m_endElemX; ei++)
    {
        const double rho_bar = m_Uh[ei][0];
        const double m_bar   = m_Uh[ei][1];
        const double E_bar   = m_Uh[ei][2];
        const double p_bar   = (gamma - 1.0) * (E_bar - 0.5 * m_bar * m_bar / rho_bar);

        if (rho_bar < eps || p_bar < eps) continue;

        // === Step 1: Density limiter ===
        // Scale density reconstruction so that rho >= eps at all face points,
        // while preserving the cell average.
        const double rho_min = std::min(m_cellFlux[ei][0][0], m_cellFlux[ei][1][0]);

        if (rho_min < eps)
        {
            const double theta1 = std::min(1.0, (rho_bar - eps) / (rho_bar - rho_min));
            for (int gp = 0; gp < 2; gp++)
                m_cellFlux[ei][gp][0] = theta1 * (m_cellFlux[ei][gp][0] - rho_bar) + rho_bar;
        }

        // === Step 2: Pressure limiter ===
        // For each face point a, if p(q̂_a) < eps, find the convex combination
        // parameter t such that p((1-t)*w_bar + t*q̂_a) = eps by solving a
        // quadratic, then take theta2 = min_a(t). Scale the full state vector.
        double theta2 = 1.0;

        for (int gp = 0; gp < 2; gp++)
        {
            const double rho_a = m_cellFlux[ei][gp][0];
            const double mom_a = m_cellFlux[ei][gp][1];
            const double E_a   = m_cellFlux[ei][gp][2];
            const double p_a   = (gamma - 1.0) * (E_a - 0.5 * mom_a * mom_a / rho_a);

            if (p_a >= eps) continue;

            // s(t) = (1-t)*w_bar + t*q̂_a, find t in (0,1] where p(s(t)) = eps.
            // Expanding p(s(t))*rho(s(t)) = eps*rho(s(t)) gives the quadratic:
            // A*t^2 + B*t + C = 0,  C = rho_bar*(p_bar - eps) > 0
            const double Drho = rho_a - rho_bar;
            const double Dm   = mom_a - m_bar;
            const double DE   = E_a   - E_bar;

            const double A = (gamma - 1.0) * (DE * Drho - 0.5 * Dm * Dm);
            const double B = (gamma - 1.0) * (DE * rho_bar + E_bar * Drho - m_bar * Dm)
                             - eps * Drho;
            const double C = rho_bar * (p_bar - eps);

            double t_eps;
            if (std::fabs(A) < 1e-14 * (std::fabs(B) + std::fabs(C)))
            {
                t_eps = -C / B;
            }
            else
            {
                double disc = B * B - 4.0 * A * C;
                if (disc < 0.0) disc = 0.0;
                const double sqrtD = std::sqrt(disc);
                const double t1 = (-B - sqrtD) / (2.0 * A);
                const double t2 = (-B + sqrtD) / (2.0 * A);
                t_eps = 1.0;
                if (t1 > 1e-14 && t1 <= 1.0 + 1e-14) t_eps = std::min(t_eps, t1);
                if (t2 > 1e-14 && t2 <= 1.0 + 1e-14) t_eps = std::min(t_eps, t2);
            }
            theta2 = std::min(theta2, std::max(0.0, std::min(t_eps, 1.0)));
        }

        if (theta2 < 1.0)
        {
            for (int gp = 0; gp < 2; gp++)
                for (int r = 0; r < m_varNum; r++)
                    m_cellFlux[ei][gp][r] = theta2 * (m_cellFlux[ei][gp][r] - m_Uh[ei][r]) + m_Uh[ei][r];
        }
    }
}

void CWENOFV::outputAve(string prefix)
{
    int vitalVarCount = equation->getVitalVarNum();
    Array1D<double> vitalVars(vitalVarCount);
    Array1D<string> vitalVarNames(vitalVarCount);
    equation->getVitalVarName(vitalVarNames);

    string filename = "./output/average_" + prefix + ".plt";
    ofstream outputFile(filename);

    cout << "Starting to output cell averages..." << endl;

    if (!outputFile.is_open())
    {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    outputFile << "TITLE=FVSolution\n"
               << "VARIABLES=X ";
    for (int i = 0; i < vitalVarCount; i++)
        outputFile << vitalVarNames[i] << " ";
    outputFile << "\nZONE T=TA\n";

    for (int elem = m_startElemX; elem < m_endElemX; ++elem)
    {
        Array1D<double> elementSolution(m_varNum);
        for (int var = 0; var < m_varNum; ++var)
            elementSolution[var] = m_Uh[elem][var];

        equation->getVitalVarVal(elementSolution, vitalVars);

        outputFile << m_grids[elem].m_xC << " ";
        for (int var = 0; var < vitalVarCount; ++var)
            outputFile << setprecision(15) << setw(20) << vitalVars[var] << " ";
        outputFile << "\n";
    }

    outputFile.close();
    cout << "Average values output completed." << endl;
}

void CWENOFV::outputAccuracy(string prefix)
{
    int gpNum = 4;
    Array1D<double> gpoints_ref(gpNum), gweights_ref(gpNum);
    sc_math::GaussLegendre_ref(gpNum, gpoints_ref, gweights_ref);

    Array1D<double> Uexact(m_elemNum + 2 * m_ghostNum);
    Array1D<double> Uerr(m_elemNum + 2 * m_ghostNum);

    double err_1(0), err_2(0), err_inf(0);
    double gPoint, gWeight;

    Uexact.setZero();
    for (int e = m_startElemX; e != m_endElemX; ++e)
    {
        for (int gp = 0; gp != gpNum; ++gp)
        {
            gPoint = m_grids[e].m_xC + 0.5 * m_grids[e].m_xD * gpoints_ref[gp];
            gWeight = 0.5 * m_grids[e].m_xD * gweights_ref[gp];

            const double Uref = equation->theVarExact(gPoint, m_now);
            Uexact[e] += Uref * gWeight / m_grids[e].m_xD;
        }
    }

    Array1D<double> Uave(m_varNum);
    for (int e = m_startElemX; e != m_endElemX; ++e)
    {
        for (int r = 0; r != m_varNum; ++r)
            Uave[r] = m_Uh[e][r];
        Uerr[e] = fabs(Uexact[e] - equation->theVarUh(Uave));

        err_1 += Uerr[e];
        err_2 += pow(Uerr[e], 2);
        if (Uerr[e] > err_inf)
            err_inf = Uerr[e];
    }

    err_1 = err_1 / m_elemNum;
    err_2 = sqrt(err_2 / m_elemNum);

    std::fstream fileout;
    string elemNum = sc_common::intToString(m_elemNum);
    string filename = "./output/accuracy_" + prefix + ".csv";

    std::ifstream fileExists(filename.c_str());
    if (fileExists)
    {
        fileout.open(filename.c_str(), ios::out | ios::app);
    }
    else
    {
        fileout.open(filename.c_str(), ios::out);
        fileout << "elemNum, Linf-norm (z1), L1-norm (z1), L2-norm (z1)" << std::endl;
    }

    fileout << elemNum << ",";
    fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_inf << ", ";
    fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_1 << ", ";
    fileout << setprecision(15) << setw(20) << setiosflags(ios::scientific) << err_2;
    fileout << std::endl;
    fileout.close();
}
