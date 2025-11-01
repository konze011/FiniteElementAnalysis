% Frames and Grids HW 4 Q1
clc;
clearvars;

% Constants
E = 200*10^9;
G = 80 * 10^9; % Shear Modulus [Pa]
P = -80*1000; %N
T = -120*1000; %Nm
% Element specific constants
% Element 1
I1 = 1.57 * 10^(-4);
J1 = 7.85 * 10^(-5); % Polar Moment of Inertia [m^4]
L1 = 1; % Element Length [m]
% Element 2
I2=4.91*10^-6; 
J2=2.45*10^-6;
L2=1;

% Get stiffness matrix
[K_global, F_global, u_global] = beamElementMatrix(E, I1, G, J1, L1);
% Display result
disp('elemet 1 stiffness matrix:')
disp(K_global)

[K_global2, F_global2, u_global2] = beamElementMatrix(E, I2, G, J2, L2);
disp('element 2 stiffness matrix:')
disp(K_global2)

% compute reduced matrices
fixed_DOFs1 = [1, 2, 5];   % indices of fixed DOFs
free_DOFs1 = setdiff(1:6, fixed_DOFs1);  % remaining DOFs
K_reduced1 = K_global(free_DOFs1, free_DOFs1);
disp('element 1 reduced matrix:')
disp(K_reduced1)

fixed_DOFs2 = [3, 4, 6];   % indices of fixed DOFs
free_DOFs2 = setdiff(1:6, fixed_DOFs2);  % remaining DOFs
K_reduced2 = K_global2(free_DOFs2, free_DOFs2);
disp('element 2 reduced matrix:')
disp(K_reduced2)

%combine reduced matricies
combined_matrix=K_reduced1+K_reduced2;
disp('combined reduced matricies:')
disp(combined_matrix)

%set f2, m2, and t2 matrix
knowns=[P;0;T];
unknowns=combined_matrix\knowns;
disp('solved unknowns at node 2:')
disp(unknowns)

% initialize global stiffness for 3 nodes × 3 DOFs/node = 9 DOFs
K_total = zeros(9);

% mappings that match your element DOF ordering [w1,th1,w2,th2,phi1,phi2]
map1 = [1, 2, 4, 5, 3, 6];   % element 1: node1-node2
map2 = [4, 5, 7, 8, 6, 9];   % element 2: node2-node3

% assemble element 1 (K_global is 6x6 from your beamElementMatrix)
for i = 1:6
    for j = 1:6
        K_total(map1(i), map1(j)) = K_total(map1(i), map1(j)) + K_global(i,j);
    end
end

% assemble element 2 (K_global2 is 6x6)
for i = 1:6
    for j = 1:6
        K_total(map2(i), map2(j)) = K_total(map2(i), map2(j)) + K_global2(i,j);
    end
end

disp('Assembled global stiffness matrix (9x9):')
disp(K_total)

solved_column=[0;0;0;unknowns(1,1);unknowns(2,1);unknowns(3,1);0;0;0];
disp('solved column')
disp(solved_column)

solved_forces=K_total*solved_column;
disp('solved forces')
disp(solved_forces)

function [K_global, F_global, u_global] = beamElementMatrix(E, I, G, J, L)
    %% Bending stiffness (4x4)
    K_bending = (E*I/L^3) * [
        12      6*L     -12     6*L;
        6*L     4*L^2   -6*L    2*L^2;
       -12     -6*L      12    -6*L;
        6*L     2*L^2   -6*L    4*L^2
    ];

    %% Torsional stiffness (2x2)
    K_torsion = (G*J/L) * [
        1  -1;
       -1   1
    ];

    %% Assemble into 6x6 global matrix
    K_global = zeros(6);
    K_global(1:4,1:4) = K_bending;
    K_global(5:6,5:6) = K_torsion;

    %% Define symbolic force and displacement vectors
    u_global = [w1; th1; w2; th2; phi1; phi2];
    F_global = [f1y; m1; f2y; m2; t1x; t2x];

end