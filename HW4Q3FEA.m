% Frames and Grids HW 4 Q3
E = 210*10^9; % [Pa]
G= 84*10^9;
A = 1e-2; %[m^2]
I = 2e-4; %[m^4]
J=1*10^-4;
L = 5; % length of all elements
f3x=0;
f3y=-20*1000;
m3=0;
knowns=[f3x;f3y;m3];

alpha_element1=0; %vertical element [deg]
alpha_element2=-60; %horizontal element [deg]
alpha_element3=-120;
alpha1=deg2rad(alpha_element1);
alpha2=deg2rad(alpha_element2);
alpha3=deg2rad(alpha_element3);


%calcuate K Hat (local K)
Khat = calcKhat(E, A, I, L); % Calculate Khat for all elements
disp('K hat output same for elements 1, 2, & 3:')
disp(Khat)

%calculate transform matrix for both elements
T1 = transform(alpha1); % Transform matrix for element 1
T2 = transform(alpha2); % Transform matrix for element 2
T3 = transform(alpha3);
disp('Transform matrices for elements 1, 2, and 3')
tol = 1e-10;
T1(abs(T1) < tol) = 0; %accounts for numerical round off from pi approx
format short g
disp(T1)
disp(T2)
disp(T3)

%calculate global K
K1g = T1' * Khat * T1;
K2g = T2' * Khat * T2;
K3g = T3' *Khat * T3;
disp('element 1 global K:')
disp(K1g)
disp('element 2 global K:')
disp(K2g)
disp('element 3 global K:')
disp(K3g)

% Initialize global stiffness matrix
K_total = zeros(9,9);

% Define DOF mappings for each element
map1 = [1 2 3 4 5 6]; % Element 1: node 1–2
map2 = [1 2 3 7 8 9]; % Element 2: node 1–3
map3 = [4 5 6 7 8 9]; % Element 3: node 2–3

% Assembly of element stiffness matrices into global K
for i = 1:6
    for j = 1:6
        K_total(map1(i), map1(j)) = K_total(map1(i), map1(j)) + K1g(i,j);
        K_total(map2(i), map2(j)) = K_total(map2(i), map2(j)) + K2g(i,j);
        K_total(map3(i), map3(j)) = K_total(map3(i), map3(j)) + K3g(i,j);
    end
end

disp('Assembled Global Stiffness Matrix K_total:')
disp(K_total)


fixed_DOFs1 = [1, 2, 3, 4, 5, 6];   % indices of fixed DOFs
free_DOFs1 = setdiff(1:9, fixed_DOFs1);  % remaining DOFs
K_reduced = K_total(free_DOFs1, free_DOFs1);
disp('Reduced Matrix:')
disp(K_reduced)

% put known force and moment values into a matrix:
unknowns=K_reduced\knowns;
disp('Solved displacements')
disp(unknowns)
% 
% index (row,column)
 u3=unknowns(1,1);
 v3=unknowns(2,1);
 theta3=unknowns(3,1);
% 
% known displacement matrix
 displacement=[0;0;0;0;0;0;u3;v3;theta3];
 
% output forces matrix
forces=K_total*displacement;
disp('Solved Forces')
disp(forces)

%% Functions
function [Khat]=calcKhat(E,A,I,L)
% build local Khat function (returns 6x6)
Khat= [
  (E*A/L) 0 0 -(E*A/L) 0 0;
  0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
  0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
  -(E*A/L) 0 0 (E*A/L) 0 0;
  0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
  0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L
];
end

function [T]=transform(alpha)
T = [
  cos(alpha)  sin(alpha)  0   0           0           0;
 -sin(alpha)  cos(alpha)  0   0           0           0;
   0          0           1   0           0           0;
   0          0           0   cos(alpha)  sin(alpha)  0;
   0          0           0  -sin(alpha)  cos(alpha)  0;
   0          0           0   0           0           1
];
end