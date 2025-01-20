%%TEST TEST TEST
%% mass matrix test - bogus variable generation
clc ; clear ; close ;
N = 1000;
m = (1:1:N)+10;
L = 10;
d = ones(1,N);
C = ones(1,9*N)*24;

M = Mass_matrix(N,m,L,d,C);