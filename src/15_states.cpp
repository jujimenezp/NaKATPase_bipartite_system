#include "W_matrix.hpp"

int main(int argc, char **argv){

    std::ofstream output_file("results/15_states.dat");
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]),\
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]));

    W.delete_state(18);
    W.delete_state(17);
    W.delete_state(16);
    W.delete_state(15);
    output_file << W << std::endl;

    solver solv;
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);
    int i = solv.steady_state_index(eigenvalues, 3e-12);
    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 7, 8);
    double J_E2K2_out = solv.get_current(W, eigenvectors.col(i), 8, 9);
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 0, 1);
    double J_E2Na2_in = solv.get_current(W, eigenvectors.col(i), 3, 14);
    
    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (main cycle): " << J_E1PNa3_in << std::endl
                << "Current from [E2PNa+2] to [E2(Na+)_2] (secondary cycle): " << J_E2Na2_in << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main cycle): " << J_E2K2_in << std::endl
                << "Current from [E2(K+)_2] to [E1K+2] (main cycle): " << J_E2K2_out << std::endl;

    output_file.close();
    return 0;
}
