// Store the W matrix coefficients for multiple use
#pragma once

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#ifndef W_MATRIX_H_
#define W_MATRIX_H_

class W_matrix : public Eigen::MatrixXd{
    private:
    public:
        //Parameteres
        double T,V,F,R;
        double FVoRT;

        // Initial transition rates
        double k_1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0, \
            ko_dNv0, ko_bNv0, ko_dKv0, ko_bKv0, k_31,   \
            k_32, k_4f, k_4r, ki_dN1v0, ki_bN1v0,       \
            ki_dN, ki_bN, ki_dK, ki_bk;

        // Transition rates taking electrical potential into account
        double k_2f,k_2r,ko_dN1,ko_bN1,ko_dN,ko_bN,ko_dK,ko_bK,ki_dN1,ki_bN1;

        // Ion concentrations (mMol)
        double c_Na_out, c_Na_in, c_K_out, c_K_in;


        W_matrix(double T, double V, double F, double R, \
                 double c_Na_out, double c_Na_in, double c_K_out, double c_K_in, \
                 double  k_1, double  k_2fv0, double  k_2rv0, double  ko_dN1v0, double  ko_bN1v0, \
                 double ko_dNv0, double  ko_bNv0, double  ko_dKv0, double  ko_bKv0, \
                 double  k_31, double k_32, double  k_4f, double  k_4r, double  ki_dN1v0, \
                 double  ki_bN1v0, double ki_dN, double  ki_bN, double  ki_dK, double  ki_bk) : Eigen::MatrixXd(19,19)
        {
            T=T; V=V; F=F; R=R;
            k_1=k_1; k_2fv0=k_2fv0; k_2rv0=k_2rv0; ko_dN1v0=ko_dN1v0; ko_bN1v0=ko_bN1v0;
            ko_dNv0=ko_dNv0; ko_bNv0=ko_bNv0; ko_dKv0=ko_dKv0; ko_bKv0=ko_bKv0;
            k_31=k_31; k_32=k_32; k_4f=k_4f; k_4r=k_4r; ki_dN1v0=ki_dN1v0;
            ki_bN1v0=ki_bN1v0; ki_dN=ki_dN; ki_bN=ki_bN; ki_dK=ki_dK; ki_bk=ki_bk;

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

            c_Na_out=c_Na_out; c_Na_in=c_Na_in; c_K_out=c_K_out; c_K_in=c_K_in;

            //W = Eigen::MatrixXd::Zero(19,19);
            (*this).setZero();

            (*this)(0,0) = -(k_1+ki_dN1); (*this)(0,13) = ki_bN1*c_Na_in;
            (*this)(1,0) = k_1; (*this)(1,1) = -k_2f; (*this)(1,2) = k_2r;
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
        }
};

class solver{
    private:
    public:
        Eigen::EigenSolver<Eigen::MatrixXd> solver;
        void initialize(const Eigen::MatrixXd &W){solver = Eigen::EigenSolver<Eigen::MatrixXd>(W);}
        Eigen::VectorXd get_eigenvalues(const Eigen::MatrixXd &W){return solver.eigenvalues().real();}
        Eigen::MatrixXd get_eigenvectors(const Eigen::MatrixXd &W){return solver.eigenvectors().real();}

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
        double get_current(const Eigen::MatrixXd &W, const Eigen::VectorXd &v, int i, int j);
};

int solver::steady_state_index(Eigen::VectorXd &eigenvalues, double threshold){
    try{
        for(int i=0;i < eigenvalues.size(); i++){
            if(fabs(eigenvalues[i]) < threshold){return i;}
        }
        throw (9999);
    }
    catch(int error){
        std::cout << "No eigenvalue found with absolute value under " << threshold << std::endl;
        exit(1);
    }
}

double solver::get_current(const Eigen::MatrixXd &W, const Eigen::VectorXd &v, int i, int j){
    return W(j,i)*v[i]-W(i,j)*v[j];
}
#endif // W_MATRIX_H_
