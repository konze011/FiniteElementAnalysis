% Frames and Grids HW 4 Q2
E = 70e9; % [Pa]
A = 4e-2; %[m^2]
I = 2e-4; %[m^4]
L1 = 4; % vertical element 1 length
L2 = 4; % horizontal element 2 length
alpha_element1=90; %vertical element [deg]
alpha_element2=0; %horizontal element [deg]
alpha1=deg2rad(alpha_element1);
alpha2=deg2rad(alpha_element2);
f2x=30*1000; %[N]
f2y=0; %[N]
m2=30*1000; %[kN*m]
f3x=0; 
m3=0; %[kN*m]
r = 0.113;

%calcuate K Hat (local K)
Khat = calcKhat(E, A, I, L1); % Calculate Khat for element 1
disp('K hat output same for elements 1 and 2:')
disp(Khat)

%calculate transform matrix for both elements
T1 = transform(alpha1); % Transform matrix for element 1
T2 = transform(alpha2); % Transform matrix for element 2
disp('Transform matrices for elements 1 and 2:')
tol = 1e-10;
T1(abs(T1) < tol) = 0; %accounts for numerical round off from pi approx
format short g
disp(T1)
disp(T2)

%calculate global K
K1g = T1' * Khat * T1;
K2g = T2' * Khat * T2;
disp('element 1 global K:')
disp(K1g)
disp('element 2 global K:')
disp(K2g)

% assemble into global K (9x9)
K = zeros(9);
% mapping for element1 (node1-node2): global dofs [1:6]
idx1 = [1:6];
K(idx1,idx1) = K(idx1,idx1) + K1g;

% element2 (node2-node3): global dofs [4:9]
idx2 = [4:9];
K(idx2,idx2) = K(idx2,idx2) + K2g;
disp('Global combined K')
disp(K)

fixed_DOFs1 = [1, 2, 3, 8];   % indices of fixed DOFs
free_DOFs1 = setdiff(1:9, fixed_DOFs1);  % remaining DOFs
K_reduced = K(free_DOFs1, free_DOFs1);
disp('Reduced Matrix:')
disp(K_reduced)

%put known force and moment values into a matrix:
knowns=[f2x;f2y;m2;f3x;m3];
unknowns=K_reduced\knowns;
disp('Solved displacements')
disp(unknowns)

%index (row,column)
u2=unknowns(1,1);
v2=unknowns(2,1);
theta2=unknowns(3,1);
u3=unknowns(4,1);
theta3=unknowns(5,1);

% known displacement matrix
displacement=[0;0;0;u2;v2;theta2;u3;0;theta3];

% output forces matrix
forces=K*displacement;
disp('Solved Forces')
disp(forces)

sum_forces=forces(3,1)+forces(6,1)+forces(9,1);
disp('Summed Torques')
disp(sum_forces)

%% Stress Calculation
J = pi*r^4/2;

u_elem1_global=displacement(1,1);
u_elem1_local=T1*u_elem1_global;
F_elem1_local=Khat*u_elem1_local;

% ---------- inputs (from your model) ----------
% E, G, L1 (use the element 1 length), r (half depth)
% u (global 9x1 displacement vector), T1 (6x6 transform for element 1)
% idx1 = 1:6 for element 1


L = L1;           % element 1 length
idx1 = 1:6;

% 1) local element displacements 
u_elem1_global = displacement(idx1);            
u_hat = T1 * u_elem1_global;        
% 2) Axial stress (uniform)
axial_op = [-1/L, 1/L];              
sigma_axial = E * ( axial_op * [ u_hat(1); u_hat(4) ] );  
% 3) Bending curvature operator at x=0 (left node) using Hermite shape fns
B2 = [-6/L^2, -4/L, 6/L^2, -2/L];   
vb_theta = [ u_hat(2); u_hat(3); u_hat(5); u_hat(6) ];  
curvature = B2 * vb_theta;          
% 4) Bending normal stress at outer fiber points
sigma_bend_top    = - ( + r ) * E * curvature;   % at P (y = +r)
sigma_bend_bottom = - ( - r ) * E * curvature;   % at R (y = -r)
% 5) final stresses at P,Q,R,S
sigma_P = sigma_axial + sigma_bend_top;
sigma_R = sigma_axial + sigma_bend_bottom;
sigma_Q = sigma_axial;   
sigma_S = sigma_axial;
% 6) von Mises
sigma_vm_P = sqrt(sigma_P.^2 + 3 * tau_P.^2);
sigma_vm_Q = sqrt(sigma_Q.^2 + 3 * tau_Q.^2);
sigma_vm_R = sqrt(sigma_R.^2 + 3 * tau_R.^2);
sigma_vm_S = sqrt(sigma_S.^2 + 3 * tau_S.^2);
% Note: generative AI was used to assist in coding and debugging stress
% functions



% 8) Print nicely (convert to MPa)
fprintf('sigma_P = %.4f MPa, sigma_vm_P = %.4f MPa\n', sigma_P/1e6, sigma_vm_P/1e6);
fprintf('sigma_Q = %.4f MPa, sigma_vm_Q = %.4f MPa\n', sigma_Q/1e6, sigma_vm_Q/1e6);
fprintf('sigma_R = %.4f MPa, sigma_vm_R = %.4f MPa\n', sigma_R/1e6, sigma_vm_R/1e6);
fprintf('sigma_S = %.4f MPa, sigma_vm_S = %.4f MPa\n', sigma_S/1e6, sigma_vm_S/1e6);

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