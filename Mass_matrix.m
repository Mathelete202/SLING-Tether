%% Mass Matrix generator
% The purpose of this function is to generate the mass matrix for a N-link 
% MET simulation
function [M] = Mass_matrix(N,m,L,d,C)

%Mass matrix allocation
M = zeros(12*N,12*N);

%% Filling rows corresponding to constraint m1x (center of mass constraint)
%Filling cells for acceleration P_mx
M(1:12:end,1:3:(3*N)) = repmat(m,N,1);
M(2:12:end,2:3:(3*N)) = repmat(m,N,1);
M(3:12:end,3:3:(3*N)) = repmat(m,N,1);

%Filling cells for acceleration C_m1x
M(1:12:end,(3*N+1):9:end) = repmat(m .* d,N,1);
M(2:12:end,(3*N+2):9:end) = repmat(m .* d,N,1);
M(3:12:end,(3*N+3):9:end) = repmat(m .* d,N,1);

%% Filling rows corresponding to constraints m2x (link to link pivot constraint)
%Filling cells for acceleration P_mx
for i = 4:12:(12*N)
    M(i,(3*floor(i/12)+1)) = 1;
    M((i+1),(3*floor(i/12)+2)) = 1;
    M((i+2),(3*floor(i/12)+3)) = 1;

    if (ceil(i/12) ~= N)
        M(i,(3*(1+floor(i/12))+1)) = -1;
        M((i+1),(3*(1+floor(i/12))+2)) = -1;
        M((i+2),(3*(1+floor(i/12))+3)) = -1;
    end
end

%Filling cells for acceleration C_m1x
for i = 4:12:(12*N)
    M(i,(3*N+9*floor(i/12)+1)) = L/2;
    M((i+1),(3*N+9*floor(i/12)+2)) = L/2;
    M((i+2),(3*N+9*floor(i/12)+3)) = L/2;

    if (ceil(i/12) ~= N)
        M(i,(3*N+9*(floor(i/12)+1)+1)) = -L/2;
        M((i+1),(3*N+9*(floor(i/12)+1)+2)) = -L/2;
        M((i+2),(3*N+9*(floor(i/12)+1)+3)) = -L/2;
    end
end

%% Filling rows corresponding to constraints m3x (DCM elements constraints)
%Constraint m31,m32,m33
%Since these constraints have a clear repeated pattern, we use a nested for loops
for j = 0:2 
    for i = (6+j):12:(12*N)
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+2)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C;
    end
end

%Constraint m34
for i = 9:12:(12*N)
    for j = 0:1
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+2)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C;
    end
end
%Constraint m35
for i = 10:12:(12*N)
    for j = 0:2:2
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+2)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C;
    end
end

%Constraint m36
for i = 11:12:(12*N)
    for j = 1:2
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+2)) = C;
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C;
    end
end
end
