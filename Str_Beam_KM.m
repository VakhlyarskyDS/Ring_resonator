%% Function of calculating stiffness and mass matrices of beam Bernoulli
%% element
% K_e - Stiffness matrix of the elemenet
% M_e - Mass matrix of the element
% E_t - Structure containing element parameters
% Structure of element data E_t
% E_t.x1 - x coordinate of the 1-nd node of the element
% E_t.y1 - y coordinate of the 1-st node of the element
% E_t.x2 - x coordinate of the 2-nd node of the element
% E_t.y2 - y coordinate of the 2-nd node of the element
% E_t.b - height of the element section
% E_t.h - thickness of the element section
% E_t.E - Young modulus
% E_t.p - density of material

function [K_e, M_e] = Str_Beam_KM(E_t)
% Length of the element
L = sqrt((E_t.x2-E_t.x1)^2+(E_t.y2-E_t.y1)^2);
% Area of cross section of the element
A = E_t.b*E_t.h;
% Moment of inertia of cross section of the element
J = E_t.b*E_t.h^3/12;


% Stiffness matrix
K_e = E_t.E*[A/L,          0,        0, -A/L,         0,        0;
               0,   12*J/L^3,  6*J/L^2,    0, -12*J/L^3,  6*J/L^2;
               0,    6*J/L^2,    4*J/L,    0,  -6*J/L^2,    2*J/L;
            -A/L,          0,        0,  A/L,         0,        0;
               0,  -12*J/L^3, -6*J/L^2,    0,  12*J/L^3, -6*J/L^2;
               0,    6*J/L^2,    2*J/L,    0,  -6*J/L^2,    4*J/L];

% Mass matrix
M_e = E_t.p*A*L*[1/3,         0,        0, 1/6,         0,         0;
                   0,     13/35, 11*L/210,   0,      9/70, -13*L/420;
                   0,  11*L/210,  L^2/105,   0,  13*L/420,  -L^2/140;
                 1/6,         0,        0, 1/3,         0,         0;
                   0,      9/70, 13*L/420,   0,     13/35, -11*L/210;
                   0, -13*L/420, -L^2/140,   0, -11*L/210,   L^2/105];

%% Rotation of matrices to global coordinate system
% Cos ans sin of angle of element
C_th = (E_t.x2-E_t.x1)/L;
S_th = (E_t.y2-E_t.y1)/L;
% Rotation matrix 3x3
P_m = [C_th, S_th, 0;
      -S_th, C_th, 0;
          0,    0, 1];
% Zeromatrix 3x3
O = zeros(3);
% Matrix of transformation
P_m = [P_m O; O P_m];
% Rotating matrices
K_e = P_m'*K_e*P_m;
M_e = P_m'*M_e*P_m;
% Getting symmetry part of matrices
K_e = (K_e+K_e')/2;
M_e = (M_e+M_e')/2;
