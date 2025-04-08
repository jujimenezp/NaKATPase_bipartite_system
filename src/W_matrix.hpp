// Store the transition matrix W and define functions to calculate thermodynamic quantitites
#pragma once

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#ifndef W_MATRIX_H_
#define W_MATRIX_H_

// Cass that stores the transition matrix
class W_matrix : public Eigen::MatrixXd{
    private:
    public:
        //Parameters
        double T,V,F,R; // Temperaure, voltage, faraday constant, and ideal gas constant
        double FVoRT;
        double kB, e; // Boltzmann constant and electron charge
        double prop_w01_to_w78; // W_01 / W_78
        std::string E_rate_units;

        // Initial transition rates
        double k_1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0, \
            ko_dNv0, ko_bNv0, ko_dKv0, ko_bKv0, k_31,   \
            k_32, k_4f, k_4r, ki_dN1v0, ki_bN1v0,       \
            ki_dN, ki_bN, ki_dK, ki_bk;

        // Transition rates taking electrical potential into account
        double k_2f,k_2r,ko_dN1,ko_bN1,ko_dN,ko_bN,ko_dK,ko_bK,ki_dN1,ki_bN1;

        // Ion concentrations (mMol)
        double c_Na_out, c_Na_in, c_K_out, c_K_in, c_ATP, c_ADP, c_P, K_h;


        W_matrix(double  k_1i, double  k_2fv0i, double  k_2rv0i, double  ko_dN1v0i, double  ko_bN1v0i, \
                 double ko_dNv0i, double  ko_bNv0i, double  ko_dKv0i, double  ko_bKv0i, \
                 double  k_31i, double k_32i, double  k_4fi, double  k_4ri, double  ki_dN1v0i, \
                 double  ki_bN1v0i, double ki_dNi, double  ki_bNi, double  ki_dKi, double  ki_bki, \
                 double Ti, double Vi, double Fi, double Ri,            \
                 double c_Na_outi, double c_Na_ini, double c_K_outi, double c_K_ini, \
                 double c_ATPi, double c_ADPi, double c_Pi, double K_hi, std::string J_or_eV, \
                 double c_Na_out_prop, double c_Na_in_prop, double c_K_out_prop, double c_K_in_prop, \
                 double c_ATP_prop, double c_ADP_prop, double c_P_prop, double prop_w01_to_w78i) : Eigen::MatrixXd(19,19)
        {
            T=Ti; V=Vi; F=Fi; R=Ri;
            k_1=k_1i; k_2fv0=k_2fv0i; k_2rv0=k_2rv0i; ko_dN1v0=ko_dN1v0i; ko_bN1v0=ko_bN1v0i;
            ko_dNv0=ko_dNv0i; ko_bNv0=ko_bNv0i; ko_dKv0=ko_dKv0i; ko_bKv0=ko_bKv0i;
            k_31=k_31i; k_32=k_32i; k_4f=k_4fi; k_4r=k_4ri; ki_dN1v0=ki_dN1v0i;
            ki_bN1v0=ki_bN1v0i; ki_dN=ki_dNi; ki_bN=ki_bNi; ki_dK=ki_dKi; ki_bk=ki_bki;

            FVoRT = F*V/(R*T);

            k_2f = k_2fv0*std::exp(0.1*FVoRT);
            k_2r = k_2rv0*std::exp(-0.1*FVoRT);
            ko_dN1 = ko_dN1v0*std::exp(0.65*FVoRT);
            ko_bN1 = ko_bN1v0*std::exp(-0.65*FVoRT);
            ko_dN = ko_dNv0*std::exp(0.185*FVoRT);
            ko_bN = ko_bNv0*std::exp(-0.185*FVoRT);
            ko_dK = ko_dKv0*std::exp(0.185*FVoRT);
            ko_bK = ko_bKv0*std::exp(-0.185*FVoRT);
            ki_dN1 = ki_dN1v0*std::exp(-0.25*FVoRT);
            ki_bN1 = ki_bN1v0*std::exp(0.25*FVoRT);

            c_Na_out=c_Na_outi*c_Na_out_prop; c_Na_in=c_Na_ini*c_Na_in_prop; c_K_out=c_K_outi*c_K_out_prop; c_K_in=c_K_ini*c_K_in_prop;
            c_ATP=c_ATPi*c_ATP_prop; c_ADP=c_ADPi*c_ADP_prop, c_P=c_Pi*c_P_prop, K_h=K_hi;

            // Define if the results are in J or eV
             try{
                if(J_or_eV == "J"){
                    kB = 1.380649e-23; // J/K
                    e = 1.60217663e-19; // C
                    E_rate_units = " J/s";
                }
                else if (J_or_eV == "eV"){
                    kB = 8.617333262e-5; // eV/K
                    e = 1; // eV/V
                    E_rate_units = " eV/s";
                }
                else {
                    throw(9999);
                }
            }
            catch(int error){
                std::cerr << "Invalid energy units provided to solver constructor. Please use \"J\" for Joules or \"eV\" for electronvolts." << std::endl;
                exit(1);
            }

             // construccion of the W matrix
            (*this).setZero();

            (*this)(0,0) = -(k_1*c_ATP_prop+ki_dN1); (*this)(0,13) = ki_bN1*c_Na_in;
            (*this)(1,0) = k_1*c_ATP_prop; (*this)(1,1) = -k_2f; (*this)(1,2) = k_2r;
            (*this)(2,1) = k_2f; (*this)(2,2) = -(k_2r+ko_dN1); (*this)(2,3)=ko_bN1*c_Na_out;
            (*this)(3,2) = ko_dN1; (*this)(3,3) = -(ko_bN1*c_Na_out+k_31+2*ko_dN); (*this)(3,4) = ko_bN*c_Na_out;
            (*this)(4,3) = 2*ko_dN; (*this)(4,4) = -(ko_bN*c_Na_out+ko_bK*c_K_out+ko_dN); (*this)(4,5) = 2*ko_bN*c_Na_out; (*this)(4,15) = ko_dK;
            (*this)(5,4) = ko_dN; (*this)(5,5) = -(2*ko_bN*c_Na_out+2*ko_bK*c_K_out); (*this)(5,6) = ko_dK;
            (*this)(6,5) = 2*ko_bK*c_K_out; (*this)(6,6) = -(ko_dK+ko_bN*c_Na_out+ko_bK*c_K_out); (*this)(6,7) = 2*ko_dK; (*this)(6,17) = ko_dN;
            (*this)(7,6) = ko_bK*c_K_out; (*this)(7,7) = -(2*ko_dK+k_32);
            (*this)(8,7) = k_32; (*this)(8,8) = -k_4f; (*this)(8,9) = k_4r;
            (*this)(9,8) = k_4f; (*this)(9,9) = -(k_4r+2*ki_dK); (*this)(9,10) = ki_bk*c_K_in;
            (*this)(10,9) = 2*ki_dK; (*this)(10,10) = -(ki_bk*c_K_in+ki_bN*c_Na_in+ki_dK); (*this)(10,11) = 2*ki_bk*c_K_in; (*this)(10,18) = ki_dN;
            (*this)(11,10) = ki_dK; (*this)(11,11) = -(2*ki_bk*c_K_in+2*ki_bN*c_Na_in); (*this)(11,12) = ki_dN;
            (*this)(12,11) = 2*ki_bN*c_Na_in; (*this)(12,12) = -(ki_dN+ki_bk*c_K_in+ki_bN*c_Na_in); (*this)(12,13) = 2*ki_dN; (*this)(12,16) = ki_dK;
            (*this)(13,0) = ki_dN1; (*this)(13,12) = ki_bN*c_Na_in; (*this)(13,13) = -(2*ki_dN+ki_bN1*c_Na_in); (*this)(13,14) = k_4f;
            (*this)(14,3) = k_31; (*this)(14,14) = -k_4f;
            (*this)(15,4) = ko_bK*c_K_out; (*this)(15,15) = -ko_dK;
            (*this)(16,12) = ki_bk*c_K_in; (*this)(16,16) = -ki_dK;
            (*this)(17,6) = ko_bN*c_Na_out; (*this)(17,17) = -ko_dN;
            (*this)(18,10) = ki_bN*c_Na_in; (*this)(18,18) = -ki_dN;

            // Additional transition rates to ensure thermodynamic consistency
            double w01w78=1;  prop_w01_to_w78=prop_w01_to_w78i;

            for(int i=0; i<13; i++){
                w01w78 *= (*this)(i+1,i);
                if(i!=0 && i!=7){
                    w01w78 *= 1/(*this)(i,i+1);
                }
            }
            w01w78 *= (*this)(0,13)/(*this)(13,0);
            w01w78 *= c_ADP*c_P/(c_ATP*K_h)*std::pow(c_Na_out/c_Na_in,3)*std::pow(c_K_in/c_K_out,2)*std::exp(-e*V/(kB*T));

            (*this)(0,1) = std::sqrt(prop_w01_to_w78*w01w78)*c_ADP_prop; (*this)(1,1) -= std::sqrt(prop_w01_to_w78*w01w78)*c_ADP_prop;
            (*this)(7,8) = std::sqrt(w01w78/prop_w01_to_w78)*c_P_prop; (*this)(8,8) -= std::sqrt(w01w78/prop_w01_to_w78)*c_P_prop;
            // (*this)(0,1) = 0; (*this)(1,1) -= 0;
            // (*this)(7,8) = 0; (*this)(8,8) -= 0;
        }

        // Delete row and column corresponding to state index and also
        // substract its corresponding contribution to diagonal elements
        void delete_state(int index){
            Eigen::MatrixXd A, B, C, D;

            for(int i=0; i < (*this).cols(); i++){
                (*this)(i,i) += (*this)(index,i);
            }

            A = (*this)(Eigen::seq(0,index-1), Eigen::seq(0,index-1));
            B = (*this)(Eigen::seq(0,index-1), Eigen::seq(index+1, Eigen::last));
            C = (*this)(Eigen::seq(index+1, Eigen::last), Eigen::seq(0,index-1));
            D = (*this)(Eigen::seq(index+1, Eigen::last), Eigen::seq(index+1, Eigen::last));
            (*this).resize((*this).rows()-1,(*this).cols()-1);
            (*this) << A, B,
                       C, D;
        }
};


// Find nonequilibrium steady state and calculate thermodynamic quantitites
class solver{
    private:
    public:
        // Currents
        Eigen::MatrixXd J;
        Eigen::MatrixXd J_X;
        Eigen::MatrixXd J_Y;

        //Bipartite Subsystems
        std::unordered_map<int, std::vector<int>> X;
        std::unordered_map<int, std::vector<int>> Y;

        solver(int J_size){
            J = Eigen::MatrixXd::Zero(J_size,J_size);
        }
        Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
        void initialize(const Eigen::MatrixXd &W){eigen_solver = Eigen::EigenSolver<Eigen::MatrixXd>(W);}

        // Find eigenvalues
        Eigen::VectorXd get_eigenvalues(const Eigen::MatrixXd &W){return eigen_solver.eigenvalues().real();}

        // Find eigenvectors
        Eigen::MatrixXd get_eigenvectors(const Eigen::MatrixXd &W){return eigen_solver.eigenvectors().real();}

        // Find index of eigenvalue for the steady state.
        // Does not work correctly if there are many values under the given threshold.
        int steady_state_index(Eigen::VectorXd &eigenvalues, double threshold=1e-12);

        // Normalize a vector
        Eigen::VectorXd normalize_vector(const Eigen::VectorXd v){return v/v.sum();}

        // Normalize every column of a matrix
        Eigen::MatrixXd normalize_columns(const Eigen::MatrixXd W){
            Eigen::MatrixXd result = Eigen::MatrixXd::Zero(W.rows(),W.cols());

            for(int i=0; i < W.cols(); i++){
                result.col(i) = W.col(i) / (W.col(i)).sum();
            }
            return result;
        }

        // Calculate current J_ji going from state i to state j.
        double get_current(const Eigen::MatrixXd &, const Eigen::VectorXd &, int , int ) const;
        void get_main_cycle_currents(const Eigen::MatrixXd &, const Eigen::VectorXd &);

        // Define X and Y subsystems
        void bipartite_subsystems(const W_matrix &, const Eigen::VectorXd &, const std::unordered_map<int ,std::vector<int>>, const std::unordered_map<int ,std::vector<int>> );

        // Work rate done in the 3Na_2K cycle
        double Work_3Na_2K(const W_matrix &, const Eigen::VectorXd &) const;

        // Work done through ATP hydrolysis
        double Energy_3Na_2K(const W_matrix &, const Eigen::VectorXd &) const;

        // Heat rate in the main cycle
        double Qdot(const W_matrix&, const Eigen::VectorXd&) const;

        // Entropy of the system and the environment
        double System_entropy(const W_matrix&, const Eigen::VectorXd &) const;

        // Eficiency of the transport through the 3Na_2K path
        double Efficiency_3Na_2K(const W_matrix &, const Eigen::VectorXd &) const;

        void subsystems_currents(const W_matrix &, const Eigen::VectorXd &);

        // Heat rate in bipartite system according to Ehrich and Sivak (2023)
        double Qdot_X(const W_matrix &, const Eigen::VectorXd &) const;
        double Qdot_Y(const W_matrix &, const Eigen::VectorXd &) const;

        // Work rate in bipartite system according to Ehrich and Sivak (2023)
        double Wdot_X(const W_matrix &, const Eigen::VectorXd &) const;
        double Wdot_Y(const W_matrix &, const Eigen::VectorXd &) const;

        // Information flow in bipartite system
        double Idot_X(const W_matrix &, const Eigen::VectorXd &) const;
        double Idot_Y(const W_matrix &, const Eigen::VectorXd &) const;

        // Entropy flow for subsystems
        double Sdot_X(const W_matrix &, const Eigen::VectorXd &) const;
        double Sdot_Y(const W_matrix &, const Eigen::VectorXd &) const;



        // 18 states
        // Work rate done in the 3Na_2K cycle
        double Work_3Na_2K_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Work done through ATP hydrolysis
        double Energy_3Na_2K_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Heat rate in the main cycle
        double Qdot_18states(const W_matrix&, const Eigen::VectorXd&) const;

        // Entropy of the system and the environment
        double System_entropy_18states(const W_matrix&, const Eigen::VectorXd &) const;

        // Eficiency of the transport through the 3Na_2K path
        double Efficiency_3Na_2K_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Heat rate in bipartite system according to Ehrich and Sivak (2023)
        double Qdot_X_18states(const W_matrix &, const Eigen::VectorXd &) const;
        double Qdot_Y_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Work rate in bipartite system according to Ehrich and Sivak (2023)
        double Wdot_X_18states(const W_matrix &, const Eigen::VectorXd &) const;
        double Wdot_Y_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Information flow in bipartite system
        double Idot_X_18states(const W_matrix &, const Eigen::VectorXd &) const;
        double Idot_Y_18states(const W_matrix &, const Eigen::VectorXd &) const;

        // Entropy flow for subsystems
        double Sdot_X_18states(const W_matrix &, const Eigen::VectorXd &) const;
        double Sdot_Y_18states(const W_matrix &, const Eigen::VectorXd &) const;

};

int solver::steady_state_index(Eigen::VectorXd &eigenvalues, double threshold){
    try{
        for(int i=0;i < eigenvalues.size(); i++){
            if(fabs(eigenvalues[i]) < threshold){return i;}
        }
        throw (9999);
    }
    catch(int error){
        std::cerr << "No eigenvalue found with absolute value under " << threshold << std::endl
                  << "Eigenvalues:\n" << eigenvalues << std::endl;
        exit(1);
    }
}

double solver::get_current(const Eigen::MatrixXd &W, const Eigen::VectorXd &v, int i, int j)
const{
    return W(i,j)*v[j]-W(j,i)*v[i];
}

void solver::get_main_cycle_currents(const Eigen::MatrixXd &W, const Eigen::VectorXd &v){
    for (int i=0; i < 13; i++) {
        J(i+1,i) = W(i+1,i)*v(i) - W(i,i+1)*v(i+1);
        J(i,i+1) = -J(i+1,i);
    }
    J(13,0) = W(13,0)*v(0) - W(0,13)*v(13);
    J(0,13) = -J(13,0);
}

double solver::Work_3Na_2K(const W_matrix &W, const Eigen::VectorXd &v)
const{
    double work=0;
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

    work *= W.kB*W.T;

    // Effect of the transmembrane potential
    work -= J(0,13) * W.e * W.V; //Transmembrane potential

    return work;
}

double solver::Energy_3Na_2K(const W_matrix &W, const Eigen::VectorXd &v) const{
    double G = 0;
    G += J(1,0) * std::log(W.c_ADP/(W.c_ATP*W.K_h));
    G += J(8,7) * std::log(W.c_P);
    G *= W.kB*W.T;
    
    return G;
}

double solver::Qdot(const W_matrix &W, const Eigen::VectorXd &P) const{
    double result=1;

    for (int i=0; i < 13; i++){
        result *= W(i+1,i)/W(i,i+1);
    }

    result *= W(0,13)/W(13,0);
    result = -W.kB*W.T*J(1,0)*std::log(result);

    return result;
}

double solver::System_entropy(const W_matrix &W, const Eigen::VectorXd &v) const{
    double sum=0;
    // for(auto& p: v){
    //     sum -= p*std::log(p);
    // }

    for(int i = 0; i<13; i++){
        sum -= std::log(v(i+1)/v(i));
    }
    sum -= std::log(v(0)/v(13));
    sum *= J(1,0);

    // Global entropy production
    // for(int i = 0; i<13; i++){
    //     sum += std::log(W(i+1,i)*v(i)/(W(i,i+1)*v(i+1)));
    // }
    // sum += std::log(W(0,13)*v(13)/(W(13,0)*v(0)));
    // sum *= J(1,0);
    return sum;
}


double solver::Efficiency_3Na_2K(const W_matrix &W, const Eigen::VectorXd &v) const{

    double W_Na=0, W_K=0, W_ATP=0;

    W_Na += J(3,2) * std::log(W.c_Na_out);
    W_Na += J(4,3) * std::log(W.c_Na_out);
    W_Na += J(5,4) * std::log(W.c_Na_out);
    W_Na -= J(12,11) * std::log(W.c_Na_in);
    W_Na -= J(13,12) * std::log(W.c_Na_in);
    W_Na -= J(0,13) * std::log(W.c_Na_in);
    W_Na *= W.kB*W.T;
    W_Na -= 3*J(0,13) * W.e * W.V;

    W_K -= J(6,5) * std::log(W.c_K_out);
    W_K -= J(7,6) * std::log(W.c_K_out);
    W_K += J(10,9) * std::log(W.c_K_in);
    W_K += J(11,10) * std::log(W.c_K_in);
    W_K *= W.kB*W.T;
    W_K += 2*J(0,13) * W.e * W.V;

    W_ATP += J(1,0) * std::log(W.c_ADP/(W.c_ATP*W.K_h));
    W_ATP += J(8,7) * std::log(W.c_P);
    W_ATP *= W.kB*W.T;

    std::cout << "\nW_ATP: " << W_ATP << std::endl
              << "W_Na: " << W_Na << std::endl
              << "W_K: " << W_K << std::endl;

    // Effect of the transmembrane potential
    return (W_Na+W_K)/-W_ATP;

}

void solver::subsystems_currents(const W_matrix &W, const Eigen::VectorXd &P){

    for (const auto &[i, m]: X){
        for (const auto & [j,n]: X){
            if (i==j) continue;
            for (const auto & p: m){
                for (const auto & q:n){
                    if(W(p,q) == 0 || W(q,p) == 0) continue;
                    J_X(i,j) += J(p,q);
                }
            }
        }
    }


    for (const auto &[i, m]: Y){
        for (const auto & [j,n]: Y){
            if (i==j) continue;
            for (const auto & p: m){
                for (const auto & q:n){
                    if(W(p,q) == 0 || W(q,p) == 0) continue;
                    J_Y(i,j) += J(p,q);
                }
            }
        }
    }

    return;
}

double solver::Idot_X(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Idot_x = 0, P1, P2;

    for (const auto & [i,m]: X){
        P1 = 0;
        for (const auto & p : m) P1 += P(p);
        //P1 = std::reduce(m.begin(),m.end());
        for (const auto & [j,n]: X){
            if(i<=j) continue;
            P2 = 0;
            for (const auto & q : n) P2 += P(q);
            //P2 = std::reduce(n.begin(),n.end());
            for (const auto & p: m){
                for (const auto & q: n){
                    if(J(p,q)==0) continue;
                    Idot_x += J(p,q)*std::log2((P(p)/P1)*(P2/P(q)));
                }
            }
        }
    }

    return Idot_x/J(1,0);
}

double solver::Idot_Y(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Idot_y = 0, P1, P2;

    for (const auto & [i,m]: Y){
        P1 = 0;
        for (const auto & p : m) P1 += P(p);
        //P1 = std::reduce(m.begin(),m.end());
        for (const auto & [j,n]: Y){
            if(i<=j) continue;
            P2 = 0;
            for (const auto & q : n) P2 += P(q);
            //P2 = std::reduce(n.begin(),n.end());
            for (const auto & p: m){
                for (const auto & q: n){
                    if(J(p,q)==0) continue;
                    Idot_y += J(p,q)*std::log2((P(p)/P1)*(P2/P(q)));
                }
            }
        }
    }

    return Idot_y/J(1,0);
}

double solver::Qdot_X(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Qdot_x=0;

    for (const auto &[i, m]: X){
        for (const auto & [j,n]: X){
            if (i<=j) continue;
            for (const auto & p: m){
                for (const auto & q:n){
                    if(W(p,q) == 0 || W(q,p) == 0) continue;
                    Qdot_x += J(q,p)*std::log(W(q,p)/W(p,q));
                }
            }
            //Qdot_x += J_X(j,i)*std::log(num/den);
        }
    }
    Qdot_x *= -W.kB*W.T;
    return Qdot_x;
}

double solver::Qdot_Y(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Qdot_y=0;

    for (const auto &[i, m]: Y){
        for (const auto & [j,n]: Y){
            if (i<=j) continue;
            for (const auto & p: m){
                for (const auto & q:n){
                    if(W(p,q) == 0 || W(q,p) == 0) continue;
                    Qdot_y += J(q,p)*std::log(W(q,p)/W(p,q));
                }
            }
            // Qdot_y += J_Y(j,i)*std::log(num/den);
        }
    }
    Qdot_y *= -W.kB*W.T;
    return Qdot_y;
}

double solver::Wdot_X(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Wdot_x = 0;
    double U_ATP = W.kB*W.T*std::log(W.c_ADP*W.c_P/(W.c_ATP*W.K_h))/2.;
    double U_Na1 = W.kB*W.T*std::log(W.c_Na_out) - W.e*W.V/2.;
    double U_Na2 = -W.kB*W.T*std::log(W.c_Na_in) - W.e*W.V/2.;
    double U_K1 = -W.kB*W.T*std::log(W.c_K_out) + W.e*W.V/2.;
    double U_K2 = W.kB*W.T*std::log(W.c_K_in) + W.e*W.V/2.;

    for (const auto & [i,m]: X){
        if(std::find(m.begin(),m.end(),0) != m.end() && std::find(m.begin(),m.end(),1) == m.end()) Wdot_x -= J(1,0)*U_ATP;
        if(std::find(m.begin(),m.end(),7) != m.end() && std::find(m.begin(),m.end(),8) == m.end()) Wdot_x -= J(8,7)*U_ATP;
        if(std::find(m.begin(),m.end(),2) != m.end() && std::find(m.begin(),m.end(),3) == m.end()) Wdot_x -= J(3,2)*U_Na1;
        if(std::find(m.begin(),m.end(),3) != m.end() && std::find(m.begin(),m.end(),4) == m.end()) Wdot_x -= J(4,3)*U_Na1;
        if(std::find(m.begin(),m.end(),4) != m.end() && std::find(m.begin(),m.end(),5) == m.end()) Wdot_x -= J(5,4)*U_Na1;
        if(std::find(m.begin(),m.end(),11) != m.end() && std::find(m.begin(),m.end(),12) == m.end()) Wdot_x -= J(12,11)*U_Na2;
        if(std::find(m.begin(),m.end(),12) != m.end() && std::find(m.begin(),m.end(),13) == m.end()) Wdot_x -= J(13,12)*U_Na2;
        if(std::find(m.begin(),m.end(),0) != m.end() && std::find(m.begin(),m.end(),13) == m.end()) Wdot_x -= J(0,13)*U_Na2;
        if(std::find(m.begin(),m.end(),5) != m.end() && std::find(m.begin(),m.end(),6) == m.end()) Wdot_x -= J(6,5)*U_K1;
        if(std::find(m.begin(),m.end(),6) != m.end() && std::find(m.begin(),m.end(),7) == m.end()) Wdot_x -= J(7,6)*U_K1;
        if(std::find(m.begin(),m.end(),9) != m.end() && std::find(m.begin(),m.end(),10) == m.end()) Wdot_x -= J(10,9)*U_K2;
        if(std::find(m.begin(),m.end(),10) != m.end() && std::find(m.begin(),m.end(),11) == m.end()) Wdot_x -= J(11,10)*U_K2;
    }

    return Wdot_x;
}
double solver::Wdot_Y(const W_matrix &W, const Eigen::VectorXd &P) const{
    double Wdot_y = 0;
    double U_ATP = W.kB*W.T*std::log(W.c_ADP*W.c_P/(W.c_ATP*W.K_h))/2.;
    double U_Na1 = W.kB*W.T*std::log(W.c_Na_out) - W.e*W.V/2.;
    double U_Na2 = -W.kB*W.T*std::log(W.c_Na_in) - W.e*W.V/2.;
    double U_K1 = -W.kB*W.T*std::log(W.c_K_out) + W.e*W.V/2.;
    double U_K2 = W.kB*W.T*std::log(W.c_K_in) + W.e*W.V/2.;

    for (const auto & [i,m]: Y){
        if(std::find(m.begin(),m.end(),0) != m.end() && std::find(m.begin(),m.end(),1) == m.end()) Wdot_y -= J(1,0)*U_ATP;
        if(std::find(m.begin(),m.end(),7) != m.end() && std::find(m.begin(),m.end(),8) == m.end()) Wdot_y -= J(8,7)*U_ATP;
        if(std::find(m.begin(),m.end(),2) != m.end() && std::find(m.begin(),m.end(),3) == m.end()) Wdot_y -= J(3,2)*U_Na1;
        if(std::find(m.begin(),m.end(),3) != m.end() && std::find(m.begin(),m.end(),4) == m.end()) Wdot_y -= J(4,3)*U_Na1;
        if(std::find(m.begin(),m.end(),4) != m.end() && std::find(m.begin(),m.end(),5) == m.end()) Wdot_y -= J(5,4)*U_Na1;
        if(std::find(m.begin(),m.end(),11) != m.end() && std::find(m.begin(),m.end(),12) == m.end()) Wdot_y -= J(12,11)*U_Na2;
        if(std::find(m.begin(),m.end(),12) != m.end() && std::find(m.begin(),m.end(),13) == m.end()) Wdot_y -= J(13,12)*U_Na2;
        if(std::find(m.begin(),m.end(),0) != m.end() && std::find(m.begin(),m.end(),13) == m.end()) Wdot_y -= J(0,13)*U_Na2;
        if(std::find(m.begin(),m.end(),5) != m.end() && std::find(m.begin(),m.end(),6) == m.end()) Wdot_y -= J(6,5)*U_K1;
        if(std::find(m.begin(),m.end(),6) != m.end() && std::find(m.begin(),m.end(),7) == m.end()) Wdot_y -= J(7,6)*U_K1;
        if(std::find(m.begin(),m.end(),9) != m.end() && std::find(m.begin(),m.end(),10) == m.end()) Wdot_y -= J(10,9)*U_K2;
        if(std::find(m.begin(),m.end(),10) != m.end() && std::find(m.begin(),m.end(),11) == m.end()) Wdot_y -= J(11,10)*U_K2;
    }

    return Wdot_y;
}

double solver::Sdot_X(const W_matrix &W, const Eigen::VectorXd &P) const{
    double sum = 0;

     // The marginal probabilities of being on E1 or E2P
    double P_E1=0,P_E2P=0;

    P_E1 += P(0);
    for(int i=9; i <= 13; i++){P_E1 += P(i);}
    for(int i=2; i <= 7; i++){P_E2P += P(i);}
    
    sum -= std::log(P(1)/P_E1)+std::log(P_E2P/P(1))+std::log(P(8)/P_E2P)+std::log(P_E1/P(8));
    return sum;
}

void solver::bipartite_subsystems(const W_matrix &W, const Eigen::VectorXd &P, const std::unordered_map<int ,std::vector<int>> Xi, const std::unordered_map<int ,std::vector<int>> Yi){
    X=Xi;
    Y=Yi;
    J_X = Eigen::MatrixXd::Zero(X.size(), X.size());
    J_Y = Eigen::MatrixXd::Zero(Y.size(), Y.size());
    return;
}



#endif // W_MATRIX_H
