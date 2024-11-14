#include "W_matrix.hpp"

// Find concentrations for steady state of Clarke et al.
Eigen::VectorXd Clarke_concentrations(const Eigen::MatrixXd &W);


int main(int argc, char **argv){

    std::ofstream output_file("results/19_states.dat");
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]),\
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]), std::stod(argv[28]), \
               std::stod(argv[29]),std::stod(argv[30]),std::stod(argv[31]), "eV", 1, 1, 1, 1, \
               1, 1, 1);

    output_file << W << std::endl;

    //Initialize solver
    solver solv(W.cols());
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);

    //Scaling according to Clarke et al.
    // Eigen::VectorXd P = Clarke_concentrations(W);
    // std::cout << "Clarke normalization: " << P.sum() << std::endl;
    // eigenvectors = eigenvectors*P.sum();

    std::cout << "\nEigenvalues: \n" << eigenvalues << std::endl;

    // Finding the steady state eigenvalue
    int i = solv.steady_state_index(eigenvalues, 1e-11);

    Eigen::VectorXd P_clarke  {{0.00363, 0.00469855, 9.25632e-5, 0.0005, 0.00053404, 4.88017e-5, 0.000493945, \
        0.00211696, 524.696, 85.8581, 14.3097, 0.596226, 2.23561, 2.09576, 2.22222e-5, \
        0.00503933, 26.8273, 0.000541375, 26.8306}};

    eigenvectors.col(i) = P_clarke;

    // Storing currents for main cycle
    solv.get_main_cycle_currents(W, eigenvectors.col(i));

    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    // Define bipartite subsystems
    std::unordered_map<int, std::vector<int>> X, Y;
    X[0] = {9,10,11,12,13,0};
    X[1] = {1};
    X[2] = {2,3,4,5,6,7};
    X[3] = {8,18};
    Y[0] = {0,1,2};
    Y[1] = {13,3,18};
    Y[2] = {12,4};
    Y[3] = {11,5};
    Y[4] = {10,6};
    Y[5] = {7,8,9};
    solv.bipartite_subsystems(W, eigenvectors.col(i), X, Y);

    // Cycles currents
    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 8, 7);
    double J_E2K2_out = solv.get_current(W, eigenvectors.col(i), 9, 8);
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 1, 0);
    double J_E2Na2_in = solv.get_current(W, eigenvectors.col(i), 14, 3);
    double J_E2PK_in = solv.get_current(W, eigenvectors.col(i), 6, 5);
    double J_E1Na2_in = solv.get_current(W, eigenvectors.col(i), 13, 12);

    // Dead-ends currents
    double J_E2PNaK_in = solv.get_current(W, eigenvectors.col(i), 4, 15);
    double J_E1NaK_in = solv.get_current(W, eigenvectors.col(i), 12, 16);
    double J_E2PKNa_in = solv.get_current(W, eigenvectors.col(i), 6,17);
    double J_E1KNa_in = solv.get_current(W, eigenvectors.col(i), 10, 18);

    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (Total): " << J_E1PNa3_in << std::endl
                << "Current from [E2PNa+2] to [E2(Na+)_2] (secondary path): " << J_E2Na2_in << std::endl
                << "Current from [E2P] to [E2PK+] (main path): " << J_E2PK_in << std::endl
                << "Current from [E1Na+] to [E1Na+2] (main path): " << J_E1Na2_in << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main path): " << J_E2K2_in << std::endl
                << "Current from [E2(K+)_2] to [E1K+2] (main path): " << J_E2K2_out << std::endl;
    
    output_file << "\nCurrent from [E2PNa+] to [E2PNa+K+] (dead-end): " << J_E2PNaK_in << std::endl
                << "Current from [E1Na+] to [E1Na+K+] (dead-end): " << J_E1NaK_in << std::endl
                << "Current from [E2PK+] to [E2PK+Na+] (dead-end): " << J_E2PKNa_in << std::endl
                << "Current from [E1K+] to [E1K+Na+] (dead-end): " << J_E1KNa_in << std::endl;

    // Work and heat rates in the 3Na-2K path
    double work_3Na_2K = solv.Work_3Na_2K(W, eigenvectors.col(i)) + solv.Energy_3Na_2K(W, eigenvectors.col(i));
    output_file << "\nWork rate through the 3Na_2K path: " << work_3Na_2K/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;
    double Qdot = solv.Qdot(W, eigenvectors.col(i));
    output_file << "Heat rate through the 3Na_2K path: " << Qdot/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

    // Work and heat rates for the subsystems
    double Qdot_x = solv.Qdot_X(W, eigenvectors.col(i));
    double Qdot_y = solv.Qdot_Y(W, eigenvectors.col(i));
    output_file << "\nHeat rate for bipartite system:" << std::endl
                << "Qdot_x = " << Qdot_x/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl
                << "Qdot_y = " << Qdot_y/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

    double Wdot_x = solv.Wdot_X(W, eigenvectors.col(i));
    double Wdot_y = solv.Wdot_Y(W, eigenvectors.col(i));
    output_file << "\nWork rate for bipartite system:" << std::endl
                << "Wdot_x = " << Wdot_x/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl
                << "Wdot_y = " << Wdot_y/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

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
    double I_dot_X = solv.Idot_X(W, eigenvectors.col(i));
    double I_dot_Y = solv.Idot_Y(W, eigenvectors.col(i));
    output_file << "\nInformation of bipartite system\ndI_X/dt = " << I_dot_X/std::log2(M_E) << " nats/cycle" << std::endl
                << "dI_Y/dt = " << I_dot_Y/std::log2(M_E) << " nats/cycle" << std::endl;


    output_file.close();


    // //Finding the eigenstate of Clarke, et al.

    // std::cout << "Determining the eigenstate of Clarke et al.\n"
    //     << std::setw(8) << "P(i)"<<std::setw(15)<<"Clarke et al"<<std::setw(20)<<"This simulation"<< std::setw(15)<<"Discrepancy" << std::endl;

    // //Normalization
    // //eigenvectors.col(i)= eigenvectors.col(i)*P.sum();

    // for(int j=0; j < 19; j++){
    //     std::cout << std::setw(8) << "P("+std::to_string(j)+") = " << std::setw(15) << P(j) << std::setw(20) << eigenvectors(j,i) <<std::setw(15)<< (P(j)-eigenvectors(j,i))/P(j) << std::endl;
    // }

    // std::cout << std::setw(8) << "Sum = " << std::setw(15) << P.sum() << std::setw(20) << eigenvectors.col(i).sum() << std::endl;

    return 0;
}

Eigen::VectorXd Clarke_concentrations(const Eigen::MatrixXd &W){
    Eigen::VectorXd P(19);
    double C1 = 0.002;
    double C2 = 0.724;
    double J = C1+C2;

    P(0) = J/W(1,0);
    P(3) = C1/W(14,3);
    //P[7] = C2/W(8,7);
    P(14) = C1/W(13,14);
    P(13) = (J+W(13,0)*P[0])/W(0,13);

    for(int j=12; j > 3; j--) {
        P(j) = (C2+W(j,j+1)*P[j+1])/W(j+1,j);
    }
    for (int j=2; j>0; j--) {
        P(j) = (J+W(j,j+1)*P[j+1])/W(j+1,j);
    }

    P(15) = W(15,4)*P(4)/W(4,15);
    P(16) = W(16,12)*P(12)/W(12,16);
    P(17) = W(17,6)*P(6)/W(6,17);
    P(18) = W(18,10)*P(10)/W(10,18);
    return P;
}
