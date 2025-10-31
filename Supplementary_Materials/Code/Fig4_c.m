clear all;close all;clc;
%% Set parameters 
N = 41;             % the number of antennas
lambda = 1;         % wave length

% Non-uniform distribution
rx  = [0.3749 0.6299 1.5302 1.8494 2.3497 2.8973 3.2995 3.8098 4.6065 ...
       5.0000 5.3749 5.6299 6.5302 6.8494 7.3497 7.8973 8.2995 8.8098 ...
       9.6065 10.000]'; 
rx  = [-flipud(rx); 0; rx].*lambda; 

num = 1801;         % the number of discrete angles
L_0 = 21.0;         % the initial mainlobe beamwidth
Gain_min = 6.0;     % the dBi of G_{min}
thetal = 90;        % the central angle
rho = 10^(-12/10);  % SLL
iterMax = 8;        % the maximum iteration times of Loop A (I_{m})
eta = 0.955;        % the degrading factor of G
Delta_max = 0.03;   % the upper bound of ||x_{Delta}||

%% Implement the proposed algorithm 
[gp_propose, tht, w, L_ML] = proposed_algorithm(N,lambda,rx,num,L_0,Gain_min,thetal,rho,iterMax,eta,Delta_max);