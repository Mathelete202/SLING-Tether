%% Mass Matrix generator for A Matrix
%% WHILE THIS IS TECHNICALLY FUNCTIONAL, IT CANNOT HANDLE SPARSE MATRICES, 
% NOT USABLE FOR N > 1000 (memory constraints)
% The purpose of this function is to generate the mass matrix for a N-link 
% MET simulation
%% Inputs list
% N - number of links
% m - array containing the mass of each link (in kg)
% d - array containing the CG offset for each link (in m)

%% Important notes: 
% to optimize memory usage, utilize this function as a nested function,
% this will prevent unecessary data duplication

% Note that the indexing system used in this function is derived from the
% SLING structures indexing system (see google drive structures > 
% Tether structure simulation research (sept 23 - ongoing) > READ ME)

%% Function body
function [M] = A_matrix(N,m,d)

%Mass matrix allocation
M = zeros(12*N,12*N);

%% Filling rows corresponding to constraint m1x (center of mass constraint)
%Filling cells for acceleration P_mx
for i = 1:N
    % fill in the -m_m d_m terms multiplied by position accelerations for A
    % matrix
    M(3*i-2,3*i-2) = -m(i) * d(i);
    M(3*i-1,3*i-1) = -m(i) * d(i);
    M(3*i,3*i) = -m(i) * d(i);
    
    M(3*N+9*i-8,3*i-2) = -m(i) * d(i);
    M(3*N+9*i-7,3*i-2) = -m(i) * d(i);
    M(3*N+9*i-6,3*i-2) = -m(i) * d(i);
    M(3*N+9*i-5,3*i-1) = -m(i) * d(i);
    M(3*N+9*i-4,3*i-1) = -m(i) * d(i);
    M(3*N+9*i-3,3*i-1) = -m(i) * d(i);
    M(3*N+9*i-2,3*i) = -m(i) * d(i);
    M(3*N+9*i-1,3*i) = -m(i) * d(i);
    M(3*N+9*i,3*i) = -m(i) * d(i);
    
    % fill in the -m_m d_m terms multiplied by DCM accelerations for A
    % matrix
    M(3*i-2,3*N+9*i-8) = -m(i) * d(i) * d(i);
    M(3*i-1,3*N+9*i-7) = -m(i) * d(i) * d(i);
    M(3*i,3*N+9*i-6) = -m(i) * d(i) * d(i);
    
    M(3*N+9*i-8,3*N+9*i-8) = -m(i) * d(i) * d(i);
    M(3*N+9*i-7,3*N+9*i-7) = -m(i) * d(i) * d(i);
    M(3*N+9*i-6,3*N+9*i-6) = -m(i) * d(i) * d(i);
    M(3*N+9*i-5,3*N+9*i-5) = -m(i) * d(i) * d(i);
    M(3*N+9*i-4,3*N+9*i-4) = -m(i) * d(i) * d(i);
    M(3*N+9*i-3,3*N+9*i-3) = -m(i) * d(i) * d(i);
    M(3*N+9*i-2,3*N+9*i-2) = -m(i) * d(i) * d(i);
    M(3*N+9*i-1,3*N+9*i-1) = -m(i) * d(i) * d(i);
    M(3*N+9*i,3*N+9*i) = -m(i) * d(i) * d(i);
    
end

clear m d;

end
