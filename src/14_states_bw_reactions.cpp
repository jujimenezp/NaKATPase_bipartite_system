#include "W_matrix.hpp"

int main(int argc, char **argv){

    std::cout.precision(7);

    std::ofstream output_file("results/14_states_bw_reactions.dat");
    output_file << "W(0,1)/W(7,8);Energy_flow_X;Information_flow_X;Probability_current" << std::endl;
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
    W.delete_state(18);
    W.delete_state(17);
    W.delete_state(16);
    W.delete_state(15);
    W.delete_state(14);
    
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
    std::cout << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
              << "Normalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    // Additional transition rates to ensure thermodynamic consistency
    double w01w78=1; double prop_w01_to_w78 = 1;

    for(int i=0; i<13; i++){
        w01w78 *= W(i+1,i);
        if(i!=0 && i!=7){
            w01w78 *= 1/W(i,i+1);
        }
    }
    w01w78 *= W(0,13)/W(13,0);
    w01w78 *= W.c_ADP*W.c_P/(W.c_ATP*W.K_h)*std::pow(W.c_Na_out/W.c_Na_in,3)*std::pow(W.c_K_in/W.c_K_out,2)*std::exp(-W.e*W.V/(W.kB*W.T));

    double Wdot_x, Qdot_x, Edot_x, Idot_x;
  
    for (double j=5e-5; j <= 1e6; j += 1.*pow(10, floor(log10(j))) ) {
        std::cout << "W(0,1)/W(7,8) = " << j << std::endl;

        //Update proportion of W01 to W78
        W(1,1) += prop_w01_to_w78*std::sqrt(w01w78);
        W(8,8) += std::sqrt(w01w78)/prop_w01_to_w78;
        prop_w01_to_w78 = j;
        W(0,1) = prop_w01_to_w78*std::sqrt(w01w78); W(1,1) -= prop_w01_to_w78*std::sqrt(w01w78);
        W(7,8) = std::sqrt(w01w78)/prop_w01_to_w78; W(8,8) -= std::sqrt(w01w78)/prop_w01_to_w78;
        std::cout << W << std::endl;

        // Initialize solver
        solv.initialize(W);
        eigenvalues = solv.get_eigenvalues(W);
        eigenvectors = solv.get_eigenvectors(W);
        eigenvectors = solv.normalize_columns(eigenvectors);
        i = solv.steady_state_index(eigenvalues, 1e-11);
        solv.get_main_cycle_currents(W, eigenvectors.col(i));

        Wdot_x = solv.Wdot_X(W, eigenvectors.col(i));
        Qdot_x = solv.Qdot_X(W, eigenvectors.col(i));
        Edot_x = (Wdot_x + Qdot_x)/(W.T*W.kB);
        Idot_x = solv.Idot_X(W, eigenvectors.col(i));

        output_file << j <<";"<< Edot_x/solv.J(1,0) <<";"<< Idot_x/std::log2(M_E) << ";"<< solv.J(1,0) <<std::endl;
    }
}
