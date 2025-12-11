clear all;close all;clc;
%% Set parameters 
N       = 41;             % the number of antennas
lambda  = 1;              % wave length

% Non-uniform distribution 
rx  = [0.62, 1.04, 1.47, 1.95, 2.50, 2.98, 3.44, 4.01, 4.46, 4.95, 5.39,...
        5.88, 6.31, 6.89, 7.45, 7.99, 8.41, 8.92, 9.42, 10.00]';
rx  = [-flipud(rx); 0; rx].*lambda;

num         = 1801;         % the number of discrete angles
L_0         = 12.0;         % the initial mainlobe beamwidth
Gain_min    = 8.0;          % the dBi of G_{min}
thetal      = 90;           % the central angle
rho         = 10^(-12/10);  % SLL 
iterMax     = 8;            % the maximum iteration times of Loop A (I_{m})
eta         = 0.955;        % the degrading factor of G
Delta_max   = 0.03;         % the upper bound of ||x_{Delta}||

% Load AEP
load aEphi_01degree.mat; 

%% Implement the proposed algorithm
[gp_propose, tht, w, L_ML] = proposed_algorithm(N,lambda,rx,num,L_0,Gain_min,thetal,rho,iterMax,eta,Delta_max,aep);