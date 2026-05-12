#include "../third_party/scTools/scTools.h"

#include <iostream>
#include <fstream>
#include "weno_fv.hpp"
#include "config_reader.hpp"

int main()
{
    EulerEquation equation;
    string inputfile = "./input/config.cfg";

    // ---------------- Read Configuration File --------------//
    std::cout << "Starting to read configuration file..." << std::endl;
    CConfig config(inputfile);
    std::map<std::string, std::string> option;
    config.read(option);

    // ------------ Pass the equation and initial boundary conditions to the solver -------------- //
    CWENOFV solver(&equation);
    solver.initializeSolver(option);
    solver.initializeAve();

    solver.run();
    solver.outputAve("final");

    return 0;
}
