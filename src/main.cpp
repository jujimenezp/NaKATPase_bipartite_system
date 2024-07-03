// Calculate the eigenvalues and eigenvectors of the W transition matrix
// The matrix was constructed using eqs. from Clark et al. (2013) supplementary material
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

using Eigen::MatrixXd;
using Eigen::EigenSolver;

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
    double T = std::stod(argv[20]); // Temperature in K
    double V = std::stod(argv[21]); // Transmembrane potential in V (V_in - V_out)
    double F = std::stod(argv[22]); // Faraday constant in C Mol^-1
    double R = std::stod(argv[23]); // Ideal gas constant in J K^-1 Mol^-1
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
    W(6,5) = 2*ko_bK*c_K_out; W(6,6) = -(ko_dK+ko_bN*c_Na_out+ko_bK*c_K_out); W(6,7) = 2*ko_dK; W(6,17) = ko_dN;
    W(7,6) = ko_bK*c_K_out; W(7,7) = -(2*ko_dK+k_32);
    W(8,7) = k_32; W(8,8) = -k_4f; W(8,9) = k_4r;
    W(9,8) = k_4f; W(9,9) = -(k_4r+2*ki_dK); W(9,10) = ki_bk*c_K_in;
    W(10,9) = 2*ki_dK; W(10,10) = -(ki_bk*c_K_in+ki_bN*c_Na_in+ki_dK); W(10,11) = 2*ki_bk*c_K_in; W(10,18) = ki_dN;
    W(11,10) = ki_dK; W(11,11) = -(2*ki_bk*c_K_in+2*ki_bN*c_Na_in); W(11,12) = ki_dN;
    W(12,11) = 2*ki_bN*c_Na_in; W(12,12) = -(ki_dN+ki_bk*c_K_in+ki_bN*c_Na_in); W(12,13) = 2*ki_dN; W(12,16) = ki_dK;
    W(13,0) = ki_dN1; W(13,12) = ki_bN*c_Na_in; W(13,13) = -(2*ki_dN+ki_bN1*c_Na_in); W(13,14) = k_4f;
    W(14,3) = k_31; W(14,14) = -k_4f;
    W(15,4) = ko_bK*c_K_out; W(15,15) = -ko_dK;
    W(16,12) = ki_bk*c_K_in; W(16,16) = -ki_dK;
    W(17,6) = ko_bN*c_Na_out; W(17,17) = -ko_dN;
    W(18,10) = ki_bN*c_Na_in; W(18,18) = -ki_dN;

    std::cout << W << std::endl;

    // Eigenvalues and eigenvectors
    EigenSolver<MatrixXd> solver(W);
    double eigenvalue_ss = solver.eigenvalues()[13].real();
    Eigen::VectorXd eigenvector_ss = (solver.eigenvectors().col(13)).real();
    std::cout << "The eigenvalue for the steady state of W is: " << std::endl << eigenvalue_ss << std::endl;
    std::cout << "The eigenvector for the steady state of W is: " << std::endl << eigenvector_ss << std::endl;

    MatrixXd diff = (W * eigenvector_ss) - (eigenvalue_ss * eigenvector_ss );
    std::cout << "Difference of eigenvalue*W - eigenvalue*eigenvector" << std::endl
              << diff <<std::endl;

    std::cout << "Sum of elements of eigenvector" << std::endl
              << eigenvector_ss.sum() << std::endl;

    return 0;
}
