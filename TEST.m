%%TEST TEST TEST
%% Testing section (comment out when function is working as intended)
clc ; clear ; close ;
N = 100;
m = (1:1:N)+10;
L = 10;
d = ones(1,N);
C = 22;

M = Mass_matrix(N,m,L,d,C);