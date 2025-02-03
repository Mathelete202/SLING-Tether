
N = 100000;
offset = 0;
m = 100*ones(N,1);
d = 1*ones(N,1);
L = 10;
C = 10*ones(9*N,1);


[a,b,c] = B1_matrix(N,offset,m,d,L,C);
sparse(a,b,c)