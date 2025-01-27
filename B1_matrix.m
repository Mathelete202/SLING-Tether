%% Mass Matrix B1 generator
%% Inputs list
% N - number of links
% offset - number of rows before this matrix starts
% m - array containing the mass of each link (in kg)
% d - array containing the CG offset for each link (in m)
% L - array containing the length of each link (in m)
% C - array containing DCM elements from previous timestep


%% Function body
function [is,js,vals] = B1_matrix(N,offset,m,d,L,C)

% compute number of entries
Rmn = 1;
Pmnv = 2; % variables + lagrange multipliers
Pmnl = 3;
Pmn = Pmnv+Pmnl;
Cm1nv = 2;
Cm1nl = 6;
Cm1n = Cm1nv + Cm1nl;
Cm2nv = 2;
Cm2nl = 3;
Cm2n = Cm2nv + Cm2nl;
Cm3nv = 2;
Cm3nl = 3;
Cm3n = Cm3nv + Cm3nl;

rmLag = 6; % number of lagrange multipliers referenced from next link (inapplicable for last link)

totper = 3*Rmn + 3*Pmn + 3*Cm1n + 3*Cm2n + 3*Cm3n;
tot = totper * N-rmLag;
lo = 12*N; % lagrange multiplier offset (how many P and C variables)

% allocate 3 zeros vectors to store i, j, value
is = zeros(tot,1);
js = zeros(tot,1);
vals = zeros(tot,1);

%% Filling rows corresponding to constraint m1x (center of mass constraint)

% compute sum of mi*di
tot_md = sum(m.*d);

for i = 1:N
    midi = m(i)*d(i);
    mididi = m(i)*d(i)*d(i);

    index_offset = (i-1)*totper;
    
    % fill in the -m_m d_m terms multiplied by position accelerations for A
    % matrix

    % positions
    curr = index_offset + 1;
    is(curr) = offset+3*i-2; js(curr) = 3*i-2;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*i-1; js(curr) = 3*i-1;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*i; js(curr) = 3*i;
    vals(curr) = -midi;

    % DCMs
    curr = curr+1;
    is(curr) = offset+3*N+9*i-8; js(curr) = 3*i-2;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-7; js(curr) = 3*i-2;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-6; js(curr) = 3*i-2;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-5; js(curr) = 3*i-1;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-4; js(curr) = 3*i-1;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-3; js(curr) = 3*i-1;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-2; js(curr) = 3*i;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-1; js(curr) = 3*i;
    vals(curr) = -midi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i; js(curr) = 3*i;
    vals(curr) = -midi;
    
    % fill in the -m_m d_m^2 terms multiplied by DCM accelerations for A
    % matrix
    
    % positions
    curr = curr+1;
    is(curr) = offset+3*i-2; js(curr) = 3*N+9*i-8;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*i-1; js(curr) = 3*N+9*i-7;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*i; js(curr) = 3*N+9*i-6;
    vals(curr) = -mididi;

    % DCMs
    % col 1
    curr = curr+1;
    is(curr) = offset+3*N+9*i-8; js(curr) = 3*N+9*i-8;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-7; js(curr) = 3*N+9*i-5;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-6; js(curr) = 3*N+9*i-2;
    vals(curr) = -mididi;

    % col 2
    curr = curr+1;
    is(curr) = offset+3*N+9*i-5; js(curr) = 3*N+9*i-7;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-4; js(curr) = 3*N+9*i-4;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-3; js(curr) = 3*N+9*i-1;
    vals(curr) = -mididi;

    % col 3
    curr = curr+1;
    is(curr) = offset+3*N+9*i-2; js(curr) = 3*N+9*i-6;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i-1; js(curr) = 3*N+9*i-3;
    vals(curr) = -mididi;

    curr = curr+1;
    is(curr) = offset+3*N+9*i; js(curr) = 3*N+9*i;
    vals(curr) = -mididi;


    % Lagrange multiplier mass matrix stuff now

    % fill in the -sum of m_i d_i terms multiplied by lagrange multipliers
    
    curr = curr+1;
    is(curr)=offset+3*i-2; js(curr)=lo+12*i-11;
    vals(curr) = -tot_md;

    curr = curr+1;
    is(curr)=offset+3*i-1; js(curr)=lo+12*i-10;
    vals(curr) = -tot_md;

    curr = curr+1;
    is(curr)=offset+3*i; js(curr)=lo+12*i-9;
    vals(curr) = -tot_md;

    curr = curr+1;
    is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*i-11;
    vals(curr) = -tot_md;

    curr = curr+1;
    is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*i-10;
    vals(curr) = -tot_md;

    curr = curr+1;
    is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*i-9;
    vals(curr) = -tot_md;

    % fill in where there is a coefficient of -1 for lagrange multipliers
    curr = curr+1;
    is(curr)=offset+3*i-2; js(curr)=lo+12*i-8;
    vals(curr) = -1;

    curr = curr+1;
    is(curr)=offset+3*i-1; js(curr)=lo+12*i-7;
    vals(curr) = -1;

    curr = curr+1;
    is(curr)=offset+3*i; js(curr)=lo+12*i-6;
    vals(curr) = -1;

    if i < N
        curr = curr+1;
        is(curr)=offset+3*i-2; js(curr)=lo+12*(i+1)-8;
        vals(curr) = -1;

        curr = curr+1;
        is(curr)=offset+3*i-1; js(curr)=lo+12*(i+1)-7;
        vals(curr) = -1;
    
        curr = curr+1;
        is(curr)=offset+3*i; js(curr)=lo+12*(i+1)-6;
        vals(curr) = -1;
    end
    
    % fill in where there is a coefficient of -L/2 for lagrange multipliers
    curr = curr+1;
    is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*i-8;
    vals(curr) = -L/2;

    curr = curr+1;
    is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*i-7;
    vals(curr) = -L/2;

    curr = curr+1;
    is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*i-6;
    vals(curr) = -L/2;

    if i < N
        curr = curr+1;
        is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*(i+1)-8;
        vals(curr) = -L/2;

        curr = curr+1;
        is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*(i+1)-7;
        vals(curr) = -L/2;
    
        curr = curr+1;
        is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*(i+1)-6;
        vals(curr) = -L/2;
    end

    % fill in where there is a coefficient of -2C for lagrange multipliers
    curr = curr+1;
    is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*i-5;
    vals(curr) = -2*C(9*i-8);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*i-5;
    vals(curr) = -2*C(9*i-6);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*i-5;
    vals(curr) = -2*C(9*i-6);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-5; js(curr)=lo+12*i-4;
    vals(curr) = -2*C(9*i-5);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-4; js(curr)=lo+12*i-4;
    vals(curr) = -2*C(9*i-4);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-3; js(curr)=lo+12*i-4;
    vals(curr) = -2*C(9*i-3);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-2; js(curr)=lo+12*i-3;
    vals(curr) = -2*C(9*i-2);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-1; js(curr)=lo+12*i-3;
    vals(curr) = -2*C(9*i-1);

    curr = curr+1;
    is(curr)=offset+3*N+9*i; js(curr)=lo+12*i-3;
    vals(curr) = -2*C(9*i);

    % fill in where there is a coefficient of -C for lagrange multipliers
    curr = curr+1;
    is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-5);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-8; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-2);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-4);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-7; js(curr)=lo+12*i-1;
    vals(curr) = -C(9*i-1);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-3);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-6; js(curr)=lo+12*i-1;
    vals(curr) = -C(9*i);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-5; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-8);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-5; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-2);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-4; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-7);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-4; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-1);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-3; js(curr)=lo+12*i-2;
    vals(curr) = -C(9*i-6);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-3; js(curr)=lo+12*i;
    vals(curr) = -C(9*i);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-2; js(curr)=lo+12*i-1;
    vals(curr) = -C(9*i-8);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-2; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-5);

    curr = curr+1;
    is(curr)=offset+3*N+9*i-1; js(curr)=lo+12*i-1;
    vals(curr) = -C(9*i-7);
    curr = curr+1;
    is(curr)=offset+3*N+9*i-1; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-4);

    curr = curr+1;
    is(curr)=offset+3*N+9*i; js(curr)=lo+12*i-1;
    vals(curr) = -C(9*i-6);
    curr = curr+1;
    is(curr)=offset+3*N+9*i; js(curr)=lo+12*i;
    vals(curr) = -C(9*i-3);

clear m d;



end
