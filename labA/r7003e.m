% -*-Octave-*-

% Known system constants from table 4, page 64
g = 9.8;
b_f = 0;
m_b = 0.381;
l_b = 0.112;
I_b = 0.00616;
m_w = 0.036;
l_w = 0.021;
I_w = 0.00000746;
R_m = 4.4;
L_m = 0;
b_m = 0;
K_e = 0.444;
K_t = 0.470;

%
% Task 3.3.1: Derive linearized EOM in state-space representation
%

g11 = l_w*(m_w+m_b)+(I_w/l_w);
g12 = l_w*l_b*m_b;
g21 = m_b*l_b;
g22 = I_b+m_b*(l_b)^2;
var_gamma = [g11, g12; g21, g22];

a11 = 0;
a12 = -((K_t*K_e)/(R_m*l_w));
a13 = 0;
a14 = ((K_t*K_e)/R_m);
a21 = 0;
a22 = (K_t*K_e)/(R_m*l_w);
a23 = m_b*l_b*g;
a24 = -((K_t*K_e)/R_m);
var_alfa = [a11 a12 a13 a14; a21 a22 a23 a24];

var_beta = [K_t/R_m; -K_t/R_m]; %K_t/R_m; -K_e/R_m

A = var_gamma\var_alfa;
B = var_gamma\var_beta;

% Inject remaining known states
a_1 = [0, 1, 0, 0];
a_2 = A(1:1,1:4);
a_3 = [0, 0, 0, 1];
a_4 = A(2:2,1:4);
A = [a_1; a_2; a_3; a_4]

% Ditto
b_1 = 0;
b_2 = B(1:1,1:1);
b_3 = 0;
b_4 = B(2:2,1:1);
B = [b_1; b_2; b_3; b_4]

C = [0, 0, 1, 0]

D = 0

%
% Task 3.4.1: Determine the transfer function to the system
%

sys = ss(A, B, C, D);
poles = pole(sys);
[zeros, gain] = zero(sys);

%
% Task 3.5.1 can we symbolically solve it with MATLAB?
% Has already been solved on paper, however.
%
