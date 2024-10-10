clear all
clc
close all

syms psi_11 psi_12 psi_13 real
syms psi_21 psi_22 psi_23 real
syms psi_31 psi_32 psi_33 real
syms psi_41 psi_42 psi_43 real

syms K_11 K_12 K_13 real 
syms K_21 K_22 K_23 real 
syms K_31 K_32 K_33 real 


Psi = [psi_11, psi_12,psi_13;
psi_21, psi_22,psi_23;
psi_31, psi_32,psi_33;
psi_41, psi_42,psi_43]';

K =[K_11 K_12 K_13 
K_21 K_22 K_23 
K_31 K_32 K_33];

Coeff = Psi' * K * Psi;

i = 3;
j = 2;

A_ij = 0;
for k = 1:3
    for l = 1:3
        A_ij = A_ij + Psi(k, i) * K(k, l) * Psi(l, j);
    end
end

disp(A_ij)

disp(Coeff(i, j))