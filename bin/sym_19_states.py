#!/usr/bin/python3
# Program to symbolically calculate the steady state probabilities for the
# 19 state Albers-Post model for the Sodium-Potassium pump

from pandas import read_csv
import sympy as sp

#transition_rates = read_csv("data/transition_rates.csv", skiprows=0, sep=',')
#parameters = read_csv("data/parameters.csv", skiprows=0, sep=',')

#T, V, F, R, c_Na_out, c_Na_in, c_K_out, c_K_in = parameters.iloc[:,1]
#k_1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0, ko_dNv0, ko_bNv0, Ko_dKv0, ko_bKv0, k_31, k_32, k_4f, k_4r, ki_dN1v0, ki_bN1v0, ki_dN, ki_bN, ki_dK, ki_bk = transition_rates.iloc[:,1]

T, V, F, R, c_Na_out, c_Na_in, c_K_out, c_K_in = sp.symbols('T V F R c_Na_out c_Na_in c_K_out c_K_in', real=True)
k_1, k_2fv0, k_2rv0, ko_dN1v0, ko_bN1v0, ko_dNv0 = sp.symbols('k1 k_2fv0 k_2rv0 ko_dN1v0 ko_bN1v0 ko_dNv0', real=True)
ko_bNv0, Ko_dKv0, ko_bKv0, k_31, k_32, k_4f = sp.symbols('ko_bNv0 Ko_dKv0 ko_bKv0 k_31 k_32 k_4f', real=True)
k_4r, ki_dN1v0, ki_bN1v0, ki_dN, ki_bN, ki_dK, ki_bk = sp.symbols('k_4r ki_dN1v0 ki_bN1v0 ki_dN ki_bN ki_dK ki_bk', real=True)
k_2f,k_2r,ko_dN1,ko_bN1,ko_dN,ko_bN,ko_dK,ko_bK,ki_dN1,ki_bN1 = sp.symbols('k_2f k_2r ko_dN1 ko_bN1 ko_dN ko_bN ko_dK ko_bK ki_dN1 ki_bN1', real=True)

W = sp.zeros(19,19)

W[0,0] = -(k_1+ki_dN1); W[0,13] = ki_bN1*c_Na_in;
W[1,0] = k_1; W[1,1] = -k_2f; W[1,2] = k_2r;
W[2,1] = k_2f; W[2,2] = -(k_2r+ko_dN1); W[2,3]=ko_bN1*c_Na_out;
W[3,2] = ko_dN1; W[3,3] = -(ko_bN1*c_Na_out+k_31+2*ko_dN); W[3,4] = ko_bN*c_Na_out;
W[4,3] = 2*ko_dN; W[4,4] = -(ko_bN*c_Na_out+ko_bK*c_K_out+ko_dN); W[4,5] = 2*ko_bN*c_Na_out; W[4,15] = ko_dK;
W[5,4] = ko_dN; W[5,5] = -(2*ko_bN*c_Na_out+2*ko_bK*c_K_out); W[5,6] = ko_dK;
W[6,5] = 2*ko_bK*c_K_out; W[6,6] = -(ko_dK+ko_bN*c_Na_out+ko_bK*c_K_out); W[6,7] = 2*ko_dK; W[6,17] = ko_dN;
W[7,6] = ko_bK*c_K_out; W[7,7] = -(2*ko_dK+k_32);
W[8,7] = k_32; W[8,8] = -k_4f; W[8,9] = k_4r;
W[9,8] = k_4f; W[9,9] = -(k_4r+2*ki_dK); W[9,10] = ki_bk*c_K_in;
W[10,9] = 2*ki_dK; W[10,10] = -(ki_bk*c_K_in+ki_bN*c_Na_in+ki_dK); W[10,11] = 2*ki_bk*c_K_in; W[10,18] = ki_dN;
W[11,10] = ki_dK; W[11,11] = -(2*ki_bk*c_K_in+2*ki_bN*c_Na_in); W[11,12] = ki_dN;
W[12,11] = 2*ki_bN*c_Na_in; W[12,12] = -(ki_dN+ki_bk*c_K_in+ki_bN*c_Na_in); W[12,13] = 2*ki_dN; W[12,16] = ki_dK;
W[13,0] = ki_dN1; W[13,12] = ki_bN*c_Na_in; W[13,13] = -(2*ki_dN+ki_bN1*c_Na_in); W[13,14] = k_4f;
W[14,3] = k_31; W[14,14] = -k_4f;
W[15,4] = ko_bK*c_K_out; W[15,15] = -ko_dK;
W[16,12] = ki_bk*c_K_in; W[16,16] = -ki_dK;
W[17,6] = ko_bN*c_Na_out; W[17,17] = -ko_dN;
W[18,10] = ki_bN*c_Na_in; W[18,18] = -ki_dN;

b = sp.zeros(19,1)
#sp.pprint(W)

p = W.solve(b)
sp.pprint(p)

#eigenvalues = W.eigenvals()
#eigenvectors = W.eigenvects()
#print(f'Eigenvalues:\n{eigenvalues}')
#print(f'Eigenvectors:\n{eigenvectors}')
