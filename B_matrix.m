%% Mass Matrix B generator
%% WHILE THIS IS TECHNICALLY FUNCTIONAL, IT CANNOT HANDLE SPARSE MATRICES, 
% NOT USABLE FOR N > 1000 (memory constraints)
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
function [M] = B_matrix(N,m,L,d,C)

%Mass matrix allocation
M = zeros(12*N,12*N);

%% Filling rows corresponding to constraint m1x (center of mass constraint)

% compute sum of mi*di
midi = sum(m.*d);

for i = 1:N
    % fill in the -m_m d_m terms multiplied by position accelerations for A
    % matrix
    M(3*i-2,12*i-11) = midi;
    M(3*i-2,12*i-8) = 1;
    if i < N
        M(3*i-2,12*(i+1)-8) = 1;
    end

    M(3*i-1,12*i-10) = midi;
    M(3*i-1,12*i-7) = 1;
    if i < N
        M(3*i-1,12*(i+1)-7) = 1;
    end
    
    M(3*i,12*i-9) = midi;
    M(3*i,12*i-6) = 1;
    if i < N
        M(3*i,12*(i+1)-6) = 1;
    end

    M(3*N+9*i-8,12*i-11) = midi;
    M(3*N+9*i-8,12*i-8) = L/2;
    if i < N
        M(3*N+9*i-8,12*(i+1)-8) = L/2;
    end
    M(3*N+9*i-8,12*i-5) = 2*C(9*i-8);
    M(3*N+9*i-8,12*i-2) = C(9*i-5);
    M(3*N+9*i-8,12*i) = C(9*i-2);
    
    
    M(3*N+9*i-7,12*i-10) = midi;
    M(3*N+9*i-7,12*i-7) = L/2;
    if i < N
        M(3*N+9*i-7,12*(i+1)-7) = L/2;
    end
    M(3*N+9*i-7,12*i-5) = 2*C(9*i-6);
    M(3*N+9*i-7,12*i-2) = C(9*i-4);
    M(3*N+9*i-7,12*i-1) = C(9*i-1);
    
    
    M(3*N+9*i-6,12*i-9) = midi;
    M(3*N+9*i-6,12*i-6) = L/2;
    if i < N
        M(3*N+9*i-6,12*(i+1)-6) = L/2;
    end
    M(3*N+9*i-6,12*i-5) = 2*C(9*i-6);
    M(3*N+9*i-6,12*i-2) = C(9*i-3);
    M(3*N+9*i-6,12*i-1) = C(9*i);


    M(3*N+9*i-5,12*i-4) = 2*C(9*i-5);
    M(3*N+9*i-5,12*i-2) = C(9*i-8);
    M(3*N+9*i-5,12*i) = C(9*i-2);

    M(3*N+9*i-4,12*i-4) = 2*C(9*i-4);
    M(3*N+9*i-4,12*i-2) = C(9*i-7);
    M(3*N+9*i-4,12*i) = C(9*i-1);

    M(3*N+9*i-3,12*i-4) = 2*C(9*i-3);
    M(3*N+9*i-3,12*i-2) = C(9*i-6);
    M(3*N+9*i-3,12*i) = C(9*i);

    M(3*N+9*i-2,12*i-3) = 2*C(9*i-2);
    M(3*N+9*i-2,12*i-1) = C(9*i-8);
    M(3*N+9*i-2,12*i) = C(9*i-5);

    M(3*N+9*i-1,12*i-3) = 2*C(9*i-1);
    M(3*N+9*i-1,12*i-1) = C(9*i-7);
    M(3*N+9*i-1,12*i) = C(9*i-4);

    M(3*N+9*i,12*i-3) = 2*C(9*i);
    M(3*N+9*i,12*i-1) = C(9*i-6);
    M(3*N+9*i,12*i) = C(9*i-3);


clear m d;



end
