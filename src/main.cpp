// Calculate the eigenvalues and eigenvectors of the W transition matrix
// The matrix was constructed using eqs. from Clark et al. (2013) supplementary material
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;

int main(int argc, char **argv){
    // Initial transition rates
    double k_1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0,\
        ko_dNv0, ko_bNv0, ko_dKv0, ko_bKv0, k_31,\
        k_32, k_4f, k_4r, ki_dN1v0, ki_bN1v0,\
        ki_dN, ki_bN, ki_dK, ki_bk;

    k_1 = std::stod(argv[1]);
    k_2fv0 = std::stod(argv[2]);
    k_2rv0 = std::stod(argv[3]);
    ko_dN1v0 = std::stod(argv[4]);
    ko_bN1v0 = std::stod(argv[5]);
    ko_dNv0 = std::stod(argv[6]);
    ko_bNv0 = std::stod(argv[7]);
    ko_dKv0 = std::stod(argv[8]);
    ko_bKv0 = std::stod(argv[9]);
    k_31 = std::stod(argv[10]);
    k_32 = std::stod(argv[11]);
    k_4f = std::stod(argv[12]);
    k_4r = std::stod(argv[13]);
    ki_dN1v0 = std::stod(argv[14]);
    ki_bN1v0 = std::stod(argv[15]);
    ki_dN = std::stod(argv[16]);
    ki_bN = std::stod(argv[17]);
    ki_dK = std::stod(argv[18]);
    ki_bk = std::stod(argv[19]);

    //Other parameteres
    double T = 297.15; // Temperature in K
    double V = -79e-3; // Transmembrane potential in V (V_in - V_out)
    double F = 9.648533212e4; // Faraday constant in C Mol^-1
    double R = 8.314462618; // Ideal gas constant in J K^-1 Mol^-1
    double FV_RT = F*V/(R*T); //Exponent of the Boltzmann factor for transition rates

    // Transition rates taking electrical potential into account
    double k_2f = k_2fv0*exp(0.1*FV_RT);
    double k_2r = k_2rv0*exp(-0.1*FV_RT);
    double ko_dN1 = ko_dN1v0*exp(0.65*FV_RT);
    double ko_bN1 = ko_bN1v0*exp(-0.65*FV_RT);
    double ko_dN = ko_dNv0*exp(0.185*FV_RT);
    double ko_bN = ko_bNv0*exp(-0.185*FV_RT);
    double ko_dK = ko_dKv0*exp(0.185*FV_RT);
    double ko_bK = ko_bKv0*exp(-0.185*FV_RT);
    double ki_dN1 = ki_dN1v0*exp(-0.25*FV_RT);
    double ki_bN1 = ki_bN1v0*exp(0.25*FV_RT);

    // Ion concentrations (mMol)
    double c_Na_out = 140;
    double c_Na_in = 15;
    double c_K_out = 4;
    double c_K_in = 120;

    // W transition matrix
    MatrixXd W = MatrixXd::Zero(19,19);
    W(0,0) = -(k_1+ki_dN1); W(0,13) = ki_bN1*c_Na_in;
    W(1,0) = k_1; W(1,1) = -k_2f; W(1,2) = k_2r;
    W(2,1) = k_2f; W(2,2) = -(k_2r+ko_dN1); W(2,3)=ko_bN1*c_Na_out;
    W(3,2) = ko_dN1; W(3,3) = -(ko_bN1*c_Na_out+k_31+2*ko_dN); W(3,4) = ko_bN*c_Na_out;
    W(4,3) = 2*ko_dN; W(4,4) = -(ko_bN*c_Na_out+ko_bK*c_K_out+ko_dN); W(4,5) = 2*ko_bN*c_Na_out; W(4,15) = ko_dK;
    W(5,4) = ko_dN; W(5,5) = -(2*ko_bN*c_Na_out+2*ko_bK*c_K_out); W(5,6) = ko_dK;


    return 0;
}
