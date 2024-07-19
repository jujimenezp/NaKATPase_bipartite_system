#include "W_matrix.hpp"

// Free energy available by ATP hydrolysis
double G_ATP(const solver &solv, const W_matrix &W, const Eigen::VectorXd &v);

int main(int argc, char **argv){

    std::ofstream output_file("results/19_states.dat");
    W_matrix W(std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]), std::stod(argv[4]),\
               std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]), std::stod(argv[8]), \
               std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11]), std::stod(argv[12]), \
               std::stod(argv[13]),std::stod(argv[14]),std::stod(argv[15]), std::stod(argv[16]), \
               std::stod(argv[17]),std::stod(argv[18]),std::stod(argv[19]), std::stod(argv[20]), \
               std::stod(argv[21]),std::stod(argv[22]),std::stod(argv[23]), std::stod(argv[24]), \
               std::stod(argv[25]),std::stod(argv[26]),std::stod(argv[27]));
    output_file << W << std::endl;
    
    solver solv;
    solv.initialize(W);
    Eigen::VectorXd eigenvalues = solv.get_eigenvalues(W);
    Eigen::MatrixXd eigenvectors = solv.get_eigenvectors(W);
    eigenvectors = solv.normalize_columns(eigenvectors);
    std::cout << "\nEigenvalues: \n" << eigenvalues << std::endl;
    int i = solv.steady_state_index(eigenvalues, 1e-11);
    output_file << "Steady state eigenvalue: " << eigenvalues[i] << std::endl
                << "\nNormalized steady state eigenvector: \n" << eigenvectors.col(i) << std::endl;

    //Normalize to obtain fluxes according to Clarke et al.
    //eigenvectors = 683.4895730078262*eigenvectors;

    // Cycles currents
    double J_E2K2_in = solv.get_current(W, eigenvectors.col(i), 7, 8);
    double J_E2K2_out = solv.get_current(W, eigenvectors.col(i), 8, 9);
    double J_E1PNa3_in = solv.get_current(W, eigenvectors.col(i), 0, 1);
    double J_E2Na2_in = solv.get_current(W, eigenvectors.col(i), 3, 14);
    double J_E2PK_in = solv.get_current(W, eigenvectors.col(i), 5, 6);
    double J_E1Na2_in = solv.get_current(W, eigenvectors.col(i), 12, 13);

    // Dead-ends currents
    double J_E2PNaK_in = solv.get_current(W, eigenvectors.col(i), 4, 15);
    double J_E1NaK_in = solv.get_current(W, eigenvectors.col(i), 12, 16);
    double J_E2PKNa_in = solv.get_current(W, eigenvectors.col(i), 6,17);
    double J_E1KNa_in = solv.get_current(W, eigenvectors.col(i), 10, 18);

    output_file << "\nCurrent from [E1Na+3] to [E1P(Na+)_3] (main cycle): " << J_E1PNa3_in << std::endl
                << "Current from [E2PNa+2] to [E2(Na+)_2] (secondary cycle): " << J_E2Na2_in << std::endl
                << "Current from [E2P] to [E2PK+] (main cycle): " << J_E2PK_in << std::endl
                << "Current from [E1Na+] to [E1Na+2] (main cycle): " << J_E1Na2_in << std::endl
                << "Current from [E2PK+2] to [E2(K+)_2] (main cycle): " << J_E2K2_in << std::endl
                << "Current from [E2(K+)_2] to [E1K+2] (main cycle): " << J_E2K2_out << std::endl;
    
    output_file << "\nCurrent from [E2PNa+] to [E2PNa+K+] (dead-end)" << J_E2PNaK_in << std::endl
                << "Current from [E1Na+] to [E1Na+K+] (dead-end)" << J_E1NaK_in << std::endl
                << "Current from [E2PK+] to [E2PK+Na+] (dead-end)" << J_E2PKNa_in << std::endl
                << "Current from [E1K+] to [E1K+Na+] (dead-end)" << J_E1KNa_in << std::endl;

    // Work done in the 3Na-2K path
    double work_3Na_2K = solv.Work_3Na_2K(W, eigenvectors.col(i));
    output_file << "\nWork done in the 3Na_2K path: " << work_3Na_2K << std::endl;

    output_file.close();

    //Trying to find the eigenstate of Clarke, et al.
    Eigen::VectorXd P(15);
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

    std::cout << "Determining the eigenstate of Clarke et al.\n"
        << std::setw(8) << "P(i)"<<std::setw(15)<<"Clarke et al"<<std::setw(20)<<"This simulation"<< std::setw(15)<<"Discrepancy" << std::endl;

    //Normalization
    eigenvectors.col(i)= eigenvectors.col(i)*683.4895730078262;

    for(int j=0; j < 15; j++){
        std::cout << std::setw(8) << "P("+std::to_string(j)+") = " << std::setw(15) << P(j) << std::setw(20) << eigenvectors(j,i) <<std::setw(15)<< (P(j)-eigenvectors(j,i))/P(j) << std::endl;
    }

    Eigen::VectorXd v = eigenvectors.col(i);

    std::cout << std::setw(8) << "Sum = " << std::setw(15) << P.sum() << std::setw(20) << v.head(15).sum() << std::endl;

    return 0;
}


double G_ATP(const solver &solv, const W_matrix &W, const Eigen::VectorXd &v){
    double G=0;


}
