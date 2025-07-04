#include "W_matrix.hpp"

int main(int argc, char **argv){

    std::cout.precision(7);

    std::ofstream output_file("results/14_states.dat");
    output_file.precision(15);

    // Create transition matrix with coefficients parsed by bin/driver.py
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]), \
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]), std::stod(argv[28]), \
               std::stod(argv[29]),std::stod(argv[30]),std::stod(argv[31]), "eV", 1, 1, 1, 1, \
               1, 1, 1, std::stod(argv[32]));

    // Dead-end states and secondary path states deleted
    W.delete_state(18);
    W.delete_state(17);
    W.delete_state(16);
    W.delete_state(15);
    W.delete_state(14);
    output_file << W << std::endl;

    // Initialize solver
    solver solv(W.cols()); // "J" for Joules or "eV" for electronvolts, size of currents matrix
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);
    std::cout << "\nEigenvalues: \n" << eigenvalues << std::endl;

    int i = solv.steady_state_index(eigenvalues, 1e-10);

    // Storing currents for main cycle
    solv.get_main_cycle_currents(W, eigenvectors.col(i));

    // Finding the steady state eigenvalue
    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    // Define bipartite subsystems
    std::unordered_map<int, std::vector<int>> X, Y;
    X[0] = {9,10,11,12,13,0};
    X[1] = {1};
    X[2] = {2,3,4,5,6,7};
    X[3] = {8};
    Y[0] = {0,1,2};
    Y[1] = {13,3};
    Y[2] = {12,4};
    Y[3] = {11,5};
    Y[4] = {10,6};
    Y[5] = {7,8,9};
    // X[0] = {9,10,11,12,13,0};
    // X[1] = {1};
    // X[2] = {2,3,4,5,6,7};
    // X[3] = {8};
    // Y[0] = {0,1,2};
    // Y[1] = {13,3};
    // Y[2] = {12,4};
    // Y[3] = {11,5};
    // Y[4] = {10,6};
    // Y[5] = {7,8,9};
    solv.bipartite_subsystems(W, eigenvectors.col(i), X, Y);
    solv.subsystems_currents(W, eigenvectors.col(i));

    // Different currents calculated
    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 8, 7);
    double J_E2K2_out = solv.get_current(W, eigenvectors.col(i), 9, 8);
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 1, 0);
    double J_E1Na2_in = solv.get_current(W, eigenvectors.col(i), 13, 12);

    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (main cycle): " << J_E1PNa3_in << " 1/s" << std::endl
                << "Current from [E1Na+] to [E1Na+2] (main cycle): " << J_E1Na2_in << " 1/s" << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main cycle): " << J_E2K2_in << " 1/s" << std::endl
                << "Current from [E2(K+)_2] to [E1K+2] (main cycle): " << J_E2K2_out << " 1/s" << std::endl;

    // Work and heat rates in the 3Na-2K path
    double work_3Na_2K = -(solv.Work_3Na_2K(W, eigenvectors.col(i)) + solv.Energy_3Na_2K(W, eigenvectors.col(i)))/(W.T*W.kB*J_E1PNa3_in);
    output_file << "\nWork rate through the 3Na_2K path: " << work_3Na_2K << " kBT/cycle" << std::endl;
    double Qdot = solv.Qdot(W, eigenvectors.col(i))/(W.T*W.kB*J_E1PNa3_in);
    output_file << "Heat rate through the 3Na_2K path: " << Qdot << " kBT/cycle" << std::endl;

    // Work and heat rates for the subsystems
    double Qdot_x = solv.Qdot_X(W, eigenvectors.col(i))/(W.T*W.kB*J_E1PNa3_in);
    double Qdot_y = solv.Qdot_Y(W, eigenvectors.col(i))/(W.T*W.kB*J_E1PNa3_in);
    output_file << "\nHeat rate for bipartite system:" << std::endl
                << "Qdot_x = " << Qdot_x << " kBT/cycle" << std::endl
                << "Qdot_y = " << Qdot_y << " kBT/cycle" << std::endl;

    double Wdot_x = solv.Wdot_X(W, eigenvectors.col(i))/(W.T*W.kB*J_E1PNa3_in);
    double Wdot_y = solv.Wdot_Y(W, eigenvectors.col(i))/(W.T*W.kB*J_E1PNa3_in);
    output_file << "\nWork rate for bipartite system:" << std::endl
                << "Wdot_x = " << Wdot_x << " kBT/cycle" << std::endl
                << "Wdot_y = " << Wdot_y << " kBT/cycle" << std::endl;

    // Entropy production
    double entropy_sys_3Na_2K = solv.System_entropy(W, eigenvectors.col(i));
    output_file << "\nEnvironment entropy rate: " << -Qdot/(W.kB*W.T) <<  "kB T/s"<< std::endl;
    output_file << "System entropy in the steady state: " << entropy_sys_3Na_2K << "   1/cycle" << std::endl;
    double sdot_x = solv.Sdot_X(W,eigenvectors.col(i));
    output_file << "Marginal entropy of X subsystem: " << sdot_x << "  1/cycle" << std::endl;

     // Efficiency
    double eff;
    eff = solv.Efficiency_3Na_2K(W, eigenvectors.col(i));
    output_file << "\nEfficiency: " << eff << std::endl;

    // Information of bipartite system
    double I_dot_X = solv.Idot_X(W, eigenvectors.col(i))/std::log2(M_E);
    double I_dot_Y = solv.Idot_Y(W, eigenvectors.col(i))/std::log2(M_E);
    output_file << "\nInformation of bipartite system\ndI_X/dt = " << I_dot_X << " nats/cycle" << std::endl
              << "dI_Y/dt = " << I_dot_Y << " nats/cycle" << std::endl;

    output_file.close();
    std::cout << "Current: " << J_E1PNa3_in << std::endl;

    std::ofstream results_file("results/concentration_range_KIn.dat", std::ofstream::app);
    results_file << W.c_K_in <<"\t"<< J_E1PNa3_in <<"\t"<< work_3Na_2K <<"\t"<< Qdot <<"\t"<< Qdot_x <<"\t"<< Qdot_y \
                 << "\t" << Wdot_x <<"\t"<< Wdot_y <<"\t"<< I_dot_X <<"\t"<< I_dot_Y <<"\t" << eff << std::endl;
    results_file.close();

    return 0;
}
