#include "W_matrix.hpp"

int main(int argc, char **argv){

    std::ofstream output_file("results/14_states.dat");
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]), \
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]), std::stod(argv[28]), \
               std::stod(argv[29]),std::stod(argv[30]),std::stod(argv[31]));

    W.delete_state(18);
    W.delete_state(17);
    W.delete_state(16);
    W.delete_state(15);
    W.delete_state(14);
    output_file << W << std::endl;

    solver solv;
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);
    std::cout << "\nEigenvalues: \n" << eigenvalues << std::endl;
    //eigenvectors = eigenvectors*628.0687752658882;
    int i = solv.steady_state_index(eigenvalues, 1e-11);

    // Finding the steady state eigenvalue
    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 7, 8);
    double J_E2K2_out = solv.get_current(W, eigenvectors.col(i), 8, 9);
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 0, 1);
    double J_E1Na2_in = solv.get_current(W, eigenvectors.col(i), 12, 13);

    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (main cycle): " << J_E1PNa3_in << std::endl
                << "Current from [E1Na+] to [E1Na+2] (main cycle): " << J_E1Na2_in << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main cycle): " << J_E2K2_in << std::endl
                << "Current from [E2(K+)_2] to [E1K+2] (main cycle): " << J_E2K2_out << std::endl;

    // Work done in the 3Na-2K path
    double work_3Na_2K = solv.Work_3Na_2K(W, eigenvectors.col(i));
    output_file << "\nWork rate in the 3Na_2K path: " << work_3Na_2K << std::endl;
    double energy_3Na_2K = solv.Energy_3Na_2K(W, eigenvectors.col(i));
    double Q = -energy_3Na_2K-work_3Na_2K;
    output_file << "Energy rate by ATP hydrolysis: " << energy_3Na_2K << std::endl
                << "Heat rate through the 3Na_2K path:" << Q << std::endl;
    double entropy_sys_3Na_2K = solv.System_entropy(eigenvectors.col(i));
    output_file << "System entropy in the steady state: " << entropy_sys_3Na_2K << std::endl;
    output_file << "Environment entropy: " << -(Q/W.T)/(solv.kB*log(2)) << std::endl;

     // Efficiency
    double eff;
    eff = solv.Efficiency_3Na_2K(W, eigenvectors.col(i));
    std::cout << "\nEfficiency: " << eff << std::endl;

    // Information of bipartite system
    double I_dot_X = solv.Idot_X(W, eigenvectors.col(i));
    std::cout << "\nInformation of bipartite system\ndI_X/dt = " << I_dot_X << std::endl;

    output_file.close();
    return 0;
}
