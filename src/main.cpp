// Calculate the eigenvalues and eigenvectors of the W matrix
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;

int main(int argc, char **argv){
    // Transition rates
    double k1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0,\
        ko_dNv0, ko_bNv0, Ko_dKv0, ko_bKv0, k_31,\
        k_32, k_4f, k_4r, ki_dN1v0, ki_bN1v0,\
        ki_dN, ki_bN, ki_dK, ki_bk;

    k1=std::stod(argv[1]);
    k_2fv0=std::stod(argv[2]);
    k_2rv0=std::stod(argv[3]);
    ko_dN1v0=std::stod(argv[4]);
    ko_bN1v0=std::stod(argv[5]);
    ko_dNv0=std::stod(argv[6]);
    ko_bNv0=std::stod(argv[7]);
    Ko_dKv0=std::stod(argv[8]);
    ko_bKv0=std::stod(argv[9]);
    k_31=std::stod(argv[10]);
    k_32=std::stod(argv[11]);
    k_4f=std::stod(argv[12]);
    k_4r=std::stod(argv[13]);
    ki_dN1v0=std::stod(argv[14]);
    ki_bN1v0=std::stod(argv[15]);
    ki_dN=std::stod(argv[16]);
    ki_bN=std::stod(argv[17]);
    ki_dK=std::stod(argv[18]);
    ki_bk=std::stod(argv[19]);
    std::cout << k1 <<std::endl;

    MatrixXd W = MatrixXd::Zero(19,19);

    return 0;
}
