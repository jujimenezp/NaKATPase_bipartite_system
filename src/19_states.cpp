#include "W_matrix.hpp"

int main(int argc, char **argv){
    W_matrix mat(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]),\
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]));

    solver solv;
    solv.initialize(mat.W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(mat.W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(mat.W);
    std::cout << mat.W << std::endl;
    std::cout << "eigenvalues: \n" << eigenvalues << std::endl;
    std::cout << "eigenvectors: \n" << eigenvectors << std::endl;

}
