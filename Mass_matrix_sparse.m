%% Mass Matrix generator V2
%% THIS VERSION USES SPARSE MATRICES, MUCH MORE OPTIMIZED 
% NOT USABLE FOR N > 1000
% The purpose of this function is to generate the mass matrix for a N-link 
% MET simulation
%% Inputs list
% N - number of links
% m - array containing the mass of each link (in kg)
% L - array containing the length of each link (in m)
% d - array containing the CG offset for each link (in m)
% C - array containing DCM elements from previous timestep

%% Important notes: 
% to optimize memory usage, utilize this function as a nested function,
% this will prevent unecessary data duplication

% Note that the indexing system used in this function is derived from the
% SLING structures indexing system (see google drive structures > 
% Tether structure simulation research (sept 23 - ongoing) > READ ME)

%% Function body
function [M] = Mass_matrix_sparse(N,m,L,d,C)

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

clear m d;

%% Filling rows corresponding to constraints m2x (link to link pivot constraint)
%Filling cells for acceleration P_mx
for i = 4:12:(12*N) %adding left pivot constraint
    M(i,(3*floor(i/12)+1)) = 1;
    M((i+1),(3*floor(i/12)+2)) = 1;
    M((i+2),(3*floor(i/12)+3)) = 1;

    if (ceil(i/12) ~= N) %adding right pivot constraint unless @ last link
        M(i,(3*(1+floor(i/12))+1)) = -1;
        M((i+1),(3*(1+floor(i/12))+2)) = -1;
        M((i+2),(3*(1+floor(i/12))+3)) = -1;
    end
end

%Filling cells for acceleration C_m1x
for i = 4:12:(12*N) %adding left pivot constraint
    M(i,(3*N+9*floor(i/12)+1)) = L/2;
    M((i+1),(3*N+9*floor(i/12)+2)) = L/2;
    M((i+2),(3*N+9*floor(i/12)+3)) = L/2;

    if (ceil(i/12) ~= N) %adding right pivot constraint unless @ last link
        M(i,(3*N+9*(floor(i/12)+1)+1)) = -L/2;
        M((i+1),(3*N+9*(floor(i/12)+1)+2)) = -L/2;
        M((i+2),(3*N+9*(floor(i/12)+1)+3)) = -L/2;
    end
end

clear L;

%% Filling rows corresponding to constraints m3x (DCM elements constraints)
%Constraint m31,m32,m33
%Since these 3 constraints have the same pattern, we use a nested for loops
for j = 0:2 
    for i = (6+j):12:(12*N)
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C(9*(ceil(i/12)-1)+1+j);
        M(i,(3*N+9*floor(i/12)+3*j+2)) =  C(9*(ceil(i/12)-1)+4+j);
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C(9*(ceil(i/12)-1)+7+j);
    end
end

%Constraint m34, m35, m36, nested for similar reasons
for i = 9:12:(12*N)
    %m34
    for j = 0:1
        M(i,(3*N+9*floor(i/12)+3*j+1)) = C(9*(ceil(i/12)-1)+j+1);
        M(i,(3*N+9*floor(i/12)+3*j+2)) = C(9*(ceil(i/12)-1)+j+4);
        M(i,(3*N+9*floor(i/12)+3*j+3)) = C(9*(ceil(i/12)-1)+7+j);
    end

    %m35
    for j = 0:2:2
        M(i+1,(3*N+9*floor(((i+1)+1)/12)+3*j+1)) = C(9*(ceil((i+1)/12)-1)+j+1);
        M(i+1,(3*N+9*floor((i+1)/12)+3*j+2)) = C(9*(ceil((i+1)/12)-1)+j+4);
        M(i+1,(3*N+9*floor((i+1)/12)+3*j+3)) = C(9*(ceil((i+1)/12)-1)+7+j);
    end

    %m36
    for j = 1:2
        M(i+2,(3*N+9*floor((i+2)/12)+3*j+1)) = C(9*(ceil((i+2)/12)-1)+j+1);
        M(i+2,(3*N+9*floor((i+2)/12)+3*j+2)) = C(9*(ceil((i+2)/12)-1)+j+4);
        M(i+2,(3*N+9*floor((i+2)/12)+3*j+3)) = C(9*(ceil((i+2)/12)-1)+7+j);
    end

end

end
