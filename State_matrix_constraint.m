%% State Matrix generator Constraint portion V1
%% WHILE THIS IS TECHNICALLY FUNCTIONAL, IT CANNOT HANDLE SPARSE MATRICES, 
% NOT USABLE FOR N > 1000 (memory constraints)
% The purpose of this function is to generate the mass matrix for a N-link 
% MET simulation
%% Inputs list
% x - state at previous timestep
% N - number of links
% m - array containing the mass of each link (in kg)
% L - array containing the length of each link (in m)
% d - array containing the CG offset for each link (in m)

%% Important notes: 
% to optimize memory usage, utilize this function as a nested function,
% this will prevent unecessary data duplication

% Note that the indexing system used in this function is derived from the
% SLING structures indexing system (see google drive structures > 
% Tether structure simulation research (sept 23 - ongoing) > READ ME)

%% Function body
function [M] = State_matrix_constraint(x, N,m,L,d)

%State matrix allocation
M = zeros(12*N,1);

for i = 1:N
    M(i) = x()
end

end
