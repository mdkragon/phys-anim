% test script for generating sound for cube collision
%

% sample duration
duration = 0.5;

% sampling frequency
fq = 44100;
dt = 1/fq;

% spring constant
%   k = (Young's Modulus) * (thickness);
% Young's Modulus for steel
Y = 200; 
% TODO: what thickness to use?
thickness = 0.01;

k = Y * thickness;


% mass of particles
%   m = density * thickness * area
% for now just assume homogenous object
% density for steel is 7.85 g/cm^3
mass = 0.785 * 3 * thickness;

% fluid damping constant (must be negative?)
%gamma = .1;
gamma = 0.00001;

% viscoelastic damping constant (must be negative?)
%n_eta = .1;
n_eta = 0.1;

% vertices
vert = [ -1 -1  1;  ... p111
         -1 -1 -1;  ... p211
          1 -1 -1;  ... p221
          1 -1  1;  ... p121
          -1 -1  1;  ... p111
         -1  1  1;  ... p112
         -1  1 -1;  ... p212
          1  1 -1; ... p222
          1  1  1;  ... p122
         -1  1  1];  ... p112
vert2 = [ 1  1 -1; ... p222
          1 -1 -1;  ... p221
          1 -1 1; ...
          1  1 1; ...
          1  1 -1; ...
          -1  1 -1; ...
          -1  -1 -1]; ...

 
% edges
%{
edges = [ 2 3 5;  ...
          4 6;    ...
          4 7;    ...
          8;      ...
          6 7;    ...
          8;      ...
          8]; 
%}

% number of vertices
%n = size(vert,1);
n = 8;
%n = 3;
%n = 4;

% K matrix (elastic force matrix)
%   TODO: this is only the 1D case (should be 3nx3n)
%K = [ 0 1 1 0 1 0 0 0 ; ...
%      0 0 0 1 0 1 0 0 ; ...
%      0 0 0 1 0 0 1 0 ; ...
%      0 0 0 0 0 0 0 1 ; ...
%      0 0 0 0 0 1 1 0 ; ...
%      0 0 0 0 0 0 0 1 ; ...
%      0 0 0 0 0 0 0 1 ; ...
%      0 0 0 0 0 0 0 0 ];
  
% reconstructed K to have += -k on diagonal for every math on col and row
K = [ -3 1 1 0 1 0 0 0 ; ...
      0 -3 0 1 0 1 0 0 ; ...
      0 0 -3 1 0 0 1 0 ; ...
      0 0 0 -3 0 0 0 1 ; ...
      0 0 0 0 -3 1 1 0 ; ...
      0 0 0 0 0 -3 0 1 ; ...
      0 0 0 0 0 0 -3 1 ; ...
      0 0 0 0 0 0 0 -3 ];
  
K = K + K'; % making symmetric matrix
K = [ K, zeros(n), zeros(n); ...
      zeros(n), K, zeros(n); ...
      zeros(n), zeros(n), K];
%K = k*K*K';
K = -1 * k * K;

%test tetrahedral
% K = [ 0 1 1 1; ...
%      1 0 1 1; ...
%      1 1 0 1; ...
%      1 1 1 0];

% following diagonal += -k
% K = [-6 1 1 1; ...
%      1 -6 1 1; ...
%      1 1 -6 1; ...
%      1 1 1 -6];
% K = [K, zeros(n), zeros(n); ...
%     zeros(n), K, zeros(n); ... 
%     zeros(n), zeros(n), K];
% K = -1 * k*K;

% test one triangle mesh
%K = [ 0 1 1; ...
%      1 0 1; ...
%      1 1 0];
%K = [ K, zeros(n), zeros(n); ...
%      zeros(n), K, zeros(n); ...
%      zeros(n), zeros(n), K];
%K = k*K;

% diagonalize K
[G D] = eig(K);
lambda = diag(D);
Ginv = inv(G);
% force positive eigen values // nope C:
%lambda = abs(lambda);

% mass matrix
M = mass .* eye(3*n);
m = diag(M);

% construct w vector
%   w_i = (-(gamma*lambda_i + n) +/- sqrt((gamma*lamba_i + n)^2 - 4*lambda_i))/2
t1 = -(gamma .* lambda + n_eta);
t2 = sqrt((gamma.*lambda + n_eta).^2 - 4.*lambda);
w_plus = (t1 + t2)./2.0;
w_minus = (t1 - t2)./2.0;

% intialize constant terms to 0;
c = zeros(3*n,1);

% force/impulse vector
%   This is given from collisions
f = zeros(3*n,1);
f(1) = 1;

% compute transformed impulse
g = Ginv * f;

% time of impact/collision
t0 = 0;

% update constants
c = c + g./((m .* (w_plus - w_minus)) .* exp(w_plus * t0));
c_bar = conj(c);


% derive mode amplitude from the mode velocity
%   amplitude is related to the particle vecolicty which 
%   is linearly related to the mode velocity
% dz/dt = c * w_plus * exp(w_plus * t) + c_bar * w_minus * exp(w_minus * t)
mode_vel = @(t) c .* w_plus .* exp(w_plus .* t) + c_bar .* w_minus .* exp(w_minus .* t);


% initialize current time
t = 0;

% TODO: should there be some amplitude scaling factor? 
%         the relation between mode vel and particle
a = mode_vel(0);


nmode = length(lambda);
nsample = fq*duration;
mode_resp = zeros([nmode, nsample]);
for ind = 1:fq*duration
  % actual time
  %t = dt*ind
  t = pi/100*ind;
  % compute amplitude
  v = mode_vel(t);
  % compute sample for each mode
  %mode_resp(:,ind) = v .* (c .* exp(w_plus*t) + c_bar .* exp(w_minus * t));
  mode_resp(:,ind) = v .* (c .* exp(w_plus .* t) + c_bar .* exp(w_minus .* t));
end

% aggregate mode responses
sample = sum(mode_resp,1);


%% for test purposes
c1 = c(1);
c1_bar = c_bar(1);
g1 = g(1);
m1 = m(1);
wp1 = w_plus(1);
wm1 = w_minus(1);
v1 = c1 * wp1 * exp(wp1 * t) + c1_bar * wm1 * exp(wm1 * t);



nmodes = size(mode_resp,1);
c = jet(nmodes);
figure;
set(gcf, 'Color', 'w');
modenames = {};
for ind = 1:nmodes
  modenames{ind} = sprintf('mode %d', ind);
  plot(mode_resp(ind,:), 'Color', c(ind,:));
  hold on;
end

%legend(modenames);

