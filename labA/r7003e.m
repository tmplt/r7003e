% -*-Octave-*-

%% Known system constants from table 4, page 64
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
%% Task 3.3: Derive linearized EOM in state-space representation
g11 = l_w*(m_w+m_b)+(I_w/l_w);
g12 = l_w*l_b*m_b;
g21 = m_b*l_b;
g22 = I_b+m_b*(l_b)^2;
gamma = [g11, g12; g21, g22];

a11 = 0;
a12 = -((K_t*K_e)/(R_m*l_w));
a13 = 0;
a14 = ((K_t*K_e)/R_m);
a21 = 0;
a22 = (K_t*K_e)/(R_m*l_w);
a23 = m_b*l_b*g;
a24 = -((K_t*K_e)/R_m);
alpha = [a11 a12 a13 a14; a21 a22 a23 a24];

beta = [K_t/R_m; -K_t/R_m];

A = gamma\alpha;
B = gamma\beta;

% Inject remaining known states
a_1 = [0, 1, 0, 0];
a_2 = A(1:1,1:4);
a_3 = [0, 0, 0, 1];
a_4 = A(2:2,1:4);
A = [a_1; a_2; a_3; a_4];

% Ditto
b_1 = 0;
b_2 = B(1:1,1:1);
b_3 = 0;
b_4 = B(2:2,1:1);
B = [b_1; b_2; b_3; b_4];

C = [0, 0, 1, 0];

D = 0;
%% Task 3.4: Determine the transfer function to the system
[z, p, k] = ss2zp(A, B, C, D);
% We have a pole and zero in s = 0; that pole is thus cancelled.
% However, remaining zero is close enough to zero to assume numerical zero.
z = [ 0 ];
p = p(2:end, 1);
% Plant TF is thus
[num, den] = zp2tf(z, p, k);
plant = tf(num, den);
%% Task 3.5 Design a PID controller stabilizing the TF
% Of our poles, two are in LHP and are thus stable.
% Our single real root on the RHP is unstable and must be moved.
% Our desired poles are thus (where the unstable pole has been moved):
dp = abs([p(p < 0, :); -70]);

den = den(1, 2:end); % we are not interested in the coeff. of s^3

% We find our PID parameters after equating them on paper
Kp = (dp(1)*dp(2) + dp(2)*dp(3) + dp(1)*dp(3) - den(2)) / k;
Ki = (dp(1)*dp(2)*dp(3) - den(3)) / k;
Kd = (sum(dp) - den(1)) / k;

% We verify
controller = pid(Kp, Ki, Kd);
system = feedback(plant, controller); % closed-loop system
[zc, pc, kc] = tf2zp(system.Numerator{1,1:11}, system.Denominator{1,1:11});
% Compare the above with [ z, p, k ]
% Visually inspect system response via:
%impulse(cs)
%% Task 3.7 Check if everything is working as it should be
Cf = eye(4);
Df = zeros(4, 2);

% model our novel input, d
push = inv(gamma) * [ l_w; l_b ];
Bf = [ 0 0; B(2) push(1); 0 0; B(4) push(2) ];

%% Task 3.8 Convert the controller to the discrete domain
% Find the bandwidth of the system.
sys_bw = bandwidth((controller*plant) / (1 + controller * plant));

% NOTE: the below assumes the while feedback system as the controller,
% which is probably incorrect.
% system = zpk(zc, pc, kc);
% system = minreal(system);
% sys_bw = bandwidth(system);
% In the above, we redeclare our system in zpk-form so that minreal works.
% minreal cancels the zero/pole at s = 0 that comes from adding system
% feedback.

% We find our sampling frequency as according to 8.3.6 in FPE7.
% See also example 8.1, page 623.
sampling_freq = sys_bw * 25;
% We convert to Hz
sampling_freq = sampling_freq / (2 * pi);
T = 1 / sampling_freq;

% We compute the discrete controller
controllerd = c2d(pid(Kp, Ki, Kd, T), T, 'zoh');

%Get the obervabilty matrix
ob = obsv(A, C)

%Check the range of ob matrix  
range(ob)

%Get controlability matrix
co = ctrb(A, B)

%Check the range of co matrix
range(co)
