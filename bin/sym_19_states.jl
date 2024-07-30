# Program to symbolically calculate the steady state probabilities for the
# 19 state Albers-Post model for the Sodium-Potassium pump

using Symbolics

@variables T V F R c_Na_out c_Na_in c_K_out c_K_in
@variables k_1 k_31  k_32  k_4f k_4r ki_dN  ki_bN  ki_dK  ki_bk
@variables k_2f k_2r ko_dN1 ko_bN1 ko_dN ko_bN ko_dK ko_bK ki_dN1 ki_bN1

#W = Matrix{Real}(undef, 19,19)
W = zeros(19,19)
W = convert(Matrix{Real},W)

W[1,1] = -(k_1+ki_dN1); W[1,14] = ki_bN1*c_Na_in;
W[2,1] = k_1; W[2,2] = -k_2f; W[2,3] = k_2r;
W[3,2] = k_2f; W[3,3] = -(k_2r+ko_dN1); W[3,4]=ko_bN1*c_Na_out;
W[4,3] = ko_dN1; W[4,4] = -(ko_bN1*c_Na_out+k_31+2*ko_dN); W[4,5] = ko_bN*c_Na_out;
W[5,4] = 2*ko_dN; W[5,5] = -(ko_bN*c_Na_out+ko_bK*c_K_out+ko_dN); W[5,6] = 2*ko_bN*c_Na_out; W[5,16] = ko_dK;
W[6,5] = ko_dN; W[6,6] = -(2*ko_bN*c_Na_out+2*ko_bK*c_K_out); W[6,7] = ko_dK;
W[7,6] = 2*ko_bK*c_K_out; W[7,7] = -(ko_dK+ko_bN*c_Na_out+ko_bK*c_K_out); W[7,8] = 2*ko_dK; W[7,18] = ko_dN;
W[8,7] = ko_bK*c_K_out; W[8,8] = -(2*ko_dK+k_32);
W[9,8] = k_32; W[9,9] = -k_4f; W[9,10] = k_4r;
W[10,9] = k_4f; W[10,10] = -(k_4r+2*ki_dK); W[10,11] = ki_bk*c_K_in;
W[11,10] = 2*ki_dK; W[11,11] = -(ki_bk*c_K_in+ki_bN*c_Na_in+ki_dK); W[11,12] = 2*ki_bk*c_K_in; W[11,19] = ki_dN;
W[12,11] = ki_dK; W[12,12] = -(2*ki_bk*c_K_in+2*ki_bN*c_Na_in); W[12,13] = ki_dN;
W[13,12] = 2*ki_bN*c_Na_in; W[13,13] = -(ki_dN+ki_bk*c_K_in+ki_bN*c_Na_in); W[13,14] = 3*ki_dN; W[13,17] = ki_dK;
W[14,1] = ki_dN1; W[14,13] = ki_bN*c_Na_in; W[14,14] = -(2*ki_dN+ki_bN1*c_Na_in); W[14,15] = k_4f;
W[15,4] = k_31; W[15,15] = -k_4f;
W[16,5] = ko_bK*c_K_out; W[16,16] = -ko_dK;
W[17,13] = ki_bk*c_K_in; W[17,17] = -ki_dK;
W[18,7] = ko_bN*c_Na_out; W[18,18] = -ko_dN;
W[19,11] = ki_bN*c_Na_in; W[19,19] = -ki_dN;

display(W)
