#include "W_matrix.hpp"

int main(int argc, char **argv){

    std::cout.precision(7);

    std::ofstream output_file("results/18_states.dat");
    output_file.precision(15);

    // Create transition matrix with coefficients parsed by bin/driver.py
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]), \
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]), std::stod(argv[28]), \
               std::stod(argv[29]),std::stod(argv[30]),std::stod(argv[31]), "eV");

    // Dead-end states and secondary path states deleted
    W.delete_state(14);
    output_file << W << std::endl;

    // Initialize solver
    solver solv(W.cols()); // "J" for Joules or "eV" for electronvolts, size of currents matrix
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);
    std::cout << "\nEigenvalues: \n" << eigenvalues << std::endl;

    int i = solv.steady_state_index(eigenvalues, 1e-11);

    // Storing currents for main cycle
    solv.get_main_cycle_currents(W, eigenvectors.col(i));

    // Finding the steady state eigenvalue
    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    // Different currents calculated
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 1, 0);
    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 8, 7);
    // Dead-ends currents
    double J_E2PNaK_in = solv.get_current(W, eigenvectors.col(i), 4, 14);
    double J_E1NaK_in = solv.get_current(W, eigenvectors.col(i), 12, 15);
    double J_E2PKNa_in = solv.get_current(W, eigenvectors.col(i), 6,16);
    double J_E1KNa_in = solv.get_current(W, eigenvectors.col(i), 10, 17);


    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (main cycle): " << J_E1PNa3_in << " 1/s" << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main cycle): " << J_E2K2_in << " 1/s" << std::endl;
    output_file << "\nCurrent from [E2PNa+] to [E2PNa+K+] (dead-end)" << J_E2PNaK_in << std::endl
                << "Current from [E1Na+] to [E1Na+K+] (dead-end)" << J_E1NaK_in << std::endl
                << "Current from [E2PK+] to [E2PK+Na+] (dead-end)" << J_E2PKNa_in << std::endl
                << "Current from [E1K+] to [E1K+Na+] (dead-end)" << J_E1KNa_in << std::endl;

    // Work and heat rates in the 3Na-2K path
    double work_3Na_2K = solv.Work_3Na_2K_18states(W, eigenvectors.col(i)) + solv.Energy_3Na_2K_18states(W, eigenvectors.col(i));
    output_file << "\nWork rate through the 3Na_2K path: " << work_3Na_2K/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;
    double Qdot = solv.Qdot_18states(W, eigenvectors.col(i));
    output_file << "Heat rate through the 3Na_2K path: " << Qdot/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

    // Work and heat rates for the subsystems
    double Qdot_x = solv.Qdot_X_18states(W, eigenvectors.col(i));
    double Qdot_y = solv.Qdot_Y_18states(W, eigenvectors.col(i));
    output_file << "\nHeat rate for bipartite system:" << std::endl
                << "Qdot_x = " << Qdot_x/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl
                << "Qdot_y = " << Qdot_y/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

    double Wdot_x = solv.Wdot_X_18states(W, eigenvectors.col(i));
    double Wdot_y = solv.Wdot_Y_18states(W, eigenvectors.col(i));
    output_file << "\nWork rate for bipartite system:" << std::endl
                << "Wdot_x = " << Wdot_x/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl
                << "Wdot_y = " << Wdot_y/(W.T*W.kB*J_E1PNa3_in) << " kBT/cycle" << std::endl;

    // Entropy production
    double entropy_sys_3Na_2K = solv.System_entropy(W, eigenvectors.col(i));
    output_file << "\nEnvironment entropy rate: " << -Qdot/(W.kB*W.T) <<  " kB T/s"<< std::endl;
    output_file << "System entropy in the steady state: " << entropy_sys_3Na_2K << "   1/cycle" << std::endl;
    double sdot_x = solv.Sdot_X(W,eigenvectors.col(i));
    output_file << "Marginal entropy of X subsystem: " << sdot_x << "  1/cycle" << std::endl;

     // Efficiency
    double eff;
    eff = solv.Efficiency_3Na_2K(W, eigenvectors.col(i));
    output_file << "\nEfficiency: " << eff << std::endl;

    // Information of bipartite system
    double I_dot_X = solv.Idot_X_18states(W, eigenvectors.col(i));
    double I_dot_Y = solv.Idot_Y_18states(W, eigenvectors.col(i));
    output_file << "\nInformation of bipartite system\ndI_X/dt = " << I_dot_X/(J_E1PNa3_in*std::log2(M_E)) << " nats/cycle" << std::endl
              << "dI_Y/dt = " << I_dot_Y/(J_E1PNa3_in*std::log2(M_E)) << " nats/cycle" << std::endl;

    output_file.close();
    return 0;
}


double solver::Work_3Na_2K_18states(const W_matrix &W, const Eigen::VectorXd &v)
const{
    double work=0;

    double J_14_4 = get_current(W, v, 14, 4);
    double J_15_12 = get_current(W, v, 15, 12);
    double J_16_6 = get_current(W, v, 16, 6);
    double J_17_10 = get_current(W, v, 17, 10);
    std::cout << "Dead_end currents:\n"
              << J_14_4 <<"\n"
              << J_15_12 << "\n"
              << J_16_6 << "\n"
              << J_17_10 << std::endl;
    //E2PNa+3 -> E2PNa+2
    work += J(3,2) * std::log(W.c_Na_out);

    //E2PNa+2 -> E2PNa+
    work += J(4,3) * std::log(W.c_Na_out);

    //E2PNa+ -> E2P
    work += J(5,4) * std::log(W.c_Na_out);

    //E2P -> E2PK+
    work -= J(6,5) * std::log(W.c_K_out);

    //E2PK+ -> E2PK+2
    work -= J(7,6) * std::log(W.c_K_out);

    //E1K+2 -> E1K+
    work += J(10,9) * std::log(W.c_K_in);

    //E1K+ -> E1
    work += J(11,10) * std::log(W.c_K_in);

    //E1 -> E1Na+
    work -= J(12,11) * std::log(W.c_Na_in);

    //E1Na+ -> E1Na+2
    work -= J(13,12) * std::log(W.c_Na_in);

    //E1Na+2 -> E1Na+3
    work -= J(0,13) * std::log(W.c_Na_in);

    //Dead-end states work
    work -= J_14_4*std::log(W.c_K_out) + J_15_12*std::log(W.c_K_in) + J_16_6*std::log(W.c_Na_out) + J_17_10*std::log(W.c_Na_in);

    work *= W.kB*W.T;

    // Effect of the transmembrane potential
    work -= J(0,13) * W.e * W.V; //Transmembrane potential

    return work;
}

double solver::Energy_3Na_2K_18states(const W_matrix &W, const Eigen::VectorXd &v) const{
    double G = 0;
    G += J(1,0) * std::log(W.c_ADP/(W.c_ATP*W.K_h));
    G += J(8,7) * std::log(W.c_P);
    G *= W.kB*W.T;

    return G;
}

double solver::Qdot_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double result=1;
    double J_14_4 = get_current(W, P, 14, 4);
    double J_15_12 = get_current(W, P, 15, 12);
    double J_16_6 = get_current(W, P, 16, 6);
    double J_17_10 = get_current(W, P, 17, 10);

    for (int i=0; i < 13; i++){
        result *= W(i+1,i)/W(i,i+1);
    }

    result *= W(0,13)/W(13,0);
    result = J(1,0)*std::log(result);
    result += J_14_4*std::log(W(14,4)/W(4,14)) + J_15_12*std::log(W(15,12)/W(12,15)) + J_16_6*std::log(W(16,6)/W(6,16)) + J_17_10*std::log(W(17,10)/W(10,17));

    result = -W.kB*W.T*result;

    return result;
}

double solver::Qdot_X_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Qdot_x=0;


    Qdot_x = J(0,1)*std::log(W(0,1)/W(1,0)) + J(2,1)*std::log(W(2,1)/W(1,2)) + J(7,8)*std::log(W(7,8)/W(8,7)) + J(9,8)*std::log(W(9,8)/W(8,9));
    Qdot_x *= -W.kB*W.T;
    return Qdot_x;
}
double solver::Qdot_Y_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Qdot_y=0;
    double J_14_4 = get_current(W, P, 14, 4);
    double J_15_12 = get_current(W, P, 15, 12);
    double J_16_6 = get_current(W, P, 16, 6);
    double J_17_10 = get_current(W, P, 17, 10);

    Qdot_y += J(0,13)*std::log(W(0,13)/W(13,0)) + J(13,12)*std::log(W(13,12)/W(12,13)) + J(12,11)*std::log(W(12,11)/W(11,12)) \
           + J(11,10)*std::log(W(11,10)/W(10,11)) + J(10,9)*std::log(W(10,9)/W(9,10));
    Qdot_y += J(2,3)*std::log(W(2,3)/W(3,2)) + J(3,4)*std::log(W(3,4)/W(4,3)) + J(4,5)*std::log(W(4,5)/W(5,4)) \
           + J(5,6)*std::log(W(5,6)/W(6,5)) + J(6,7)*std::log(W(6,7)/W(7,6));
    Qdot_y += J_14_4*std::log(W(14,4)/W(4,14)) + J_15_12*std::log(W(15,12)/W(12,15)) + J_16_6*std::log(W(16,6)/W(6,16)) \
           + J_17_10*std::log(W(17,10)/W(10,17));
    Qdot_y *= -W.kB*W.T;
    return Qdot_y;
}

double solver::Wdot_X_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double U_ATP = W.kB*W.T*std::log(W.c_ADP*W.c_P/(W.c_ATP*W.K_h))/2.;

    double Wdot_x = -J(1,0)*U_ATP - J(8,7)*U_ATP;
    return Wdot_x;
}
double solver::Wdot_Y_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double J_14_4 = get_current(W, P, 14, 4);
    double J_15_12 = get_current(W, P, 15, 12);
    double J_16_6 = get_current(W, P, 16, 6);
    double J_17_10 = get_current(W, P, 17, 10);
    
    double U_Na1 = W.kB*W.T*std::log(W.c_Na_out) - W.e*W.V/2.;
    double U_Na2 = -W.kB*W.T*std::log(W.c_Na_in) - W.e*W.V/2.;
    double U_K1 = -W.kB*W.T*std::log(W.c_K_out) + W.e*W.V/2.;
    double U_K2 = W.kB*W.T*std::log(W.c_K_in) + W.e*W.V/2.;

    double Wdot_y = -U_Na1*(J(3,2)+J(4,3)+J(5,4)+J_16_6) - U_Na2*(J(12,11)+J(13,12)+J(0,13)+J_17_10) \
             -U_K1*(J(6,5)+J(7,6)+J_14_4) - U_K2*(J(10,9)+J(11,10)+J_15_12);

    return Wdot_y;
}

// Information flow in bipartite system
double solver::Idot_X_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Idot_x = 0;
    // The marginal probabilities of being on E1 or E2P
    double P_E1=0,P_E2P=0;

    P_E1 += P(0)+P(15)+P(17);
    for(int i=9; i <= 13; i++){P_E1 += P(i);}
    P_E2P += P(14)+P(16);
    for(int i=2; i <= 7; i++){P_E2P += P(i);}

    Idot_x += J(0,1)*std::log2(P(0)/P_E1) \
           + J(2,1)*std::log2(P(2)/P_E2P) \
           + J(7,8)*std::log2(P(7)/P_E2P) \
           + J(9,8)*std::log2(P(9)/P_E1);

    return Idot_x;
}
double solver::Idot_Y_18states(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Idot_y = 0;

    // Marginal probabilities
    double P_Na3, P_Na2, P_Na, P_K, P_K2, P_0, P_NaK, P_KNa;
    P_Na3 = P(0) + P(1) + P(2);
    P_Na2 = P(13) + P(3);
    P_Na = P(12) + P(4);
    P_0 = P(11) + P(5);
    P_K = P(10) + P(6);
    P_K2 = P(9) + P(8) + P(7);
    P_NaK = P(14) + P(15);
    P_KNa = P(16) + P(17);

    Idot_y += J(0,13) * std::log2(P(0)*P_Na2/(P(13)*P_Na3)) /*Sum of X=E1*/\
           + J(13,12) * std::log2(P(13)*P_Na/(P(12)*P_Na2)) \
           + J(12,11) * std::log2(P(12)*P_0/(P(11)*P_Na)) \
           + J(11,10) * std::log2(P(11)*P_K/(P(10)*P_0)) \
           + J(10,9) * std::log2(P(10)*P_K2/(P(9)*P_K)) \
           + J(2,3) * std::log2(P(2)*P_Na2/(P(3)*P_Na3)) /*Sum of X=E2P*/\
           + J(3,4) * std::log2(P(3)*P_Na/(P(4)*P_Na2)) \
           + J(4,5) * std::log2(P(4)*P_0/(P(5)*P_Na)) \
           + J(5,6) * std::log2(P(5)*P_K/(P(6)*P_0)) \
           + J(6,7) * std::log2(P(6)*P_K2/(P(7)*P_K)) \
           + J(15,12) * std::log2(P(15)*P_Na/(P(12)*P_NaK)) /*Dead end states*/\
           + J(14,4) * std::log2(P(14)*P_Na/(P(4)*P_NaK)) \
           + J(16,6) * std::log2(P(16)*P_K/(P(6)*P_KNa)) \
           + J(17,10) * std::log2(P(17)*P_K/(P(10)*P_KNa));


    return Idot_y;
}
