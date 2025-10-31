% 20251030 by zhuoruijiang@outlook.com and swlei@uestc.edu.cn
function [gp_propose, tht, w, L_ML] = proposed_algorithm(N, lambda, rx, num, L_0, Gain_min, thetal, rho, iterMax, eta, Delta_max, aep)
%   This function implements Algorithm 1 for mainlobe beamwidth maximization.
%
%   Inputs:
%       N         - the number of antennas
%       lambda    - the wave length
%       rx        - the array distribution
%       num       - the number of discrete angles
%       L_0       - the initial AMB
%       Gain_min  - the dBi of G_{min}
%       thetal    - the central angle
%       rho       - SLL 
%       iterMax   - the maximum iteration times of Loop A (I_{m})
%       eta       - the degrading factor of G
%       Delta_max - the upper bound of ||x_{Delta}||
%       aep       - (Optional) The data matrix loaded from aEphi_01degree.mat.
%                   If not provided, an ideal array model is used for calculation.
%
%   Outputs:
%       gp_propose - The final calculated power gain pattern (in dBi)
%       tht        - The corresponding angle vector (in degrees)
%       w          - The final array excitation/weight vector
%       L_ML       - The final mainlobe beamwidth  

%% Internal calculations based on input parameters
k = 2*pi/lambda;            % the wave number
tht = linspace(0,180,num)'; % the angular distribution

L_old = 0;
L_new = L_0;
cout_p = 0;                 % the count of the whole loop
cout_end = 0; 
res = 0;                    % indicator of the first convex problem

% Calculate Array Factor a
% Check if the 'aep' parameter was provided
if nargin < 12
     % If 'aep' is not provided, calculate 'a' under the ideal array model
    disp('*** "aep" parameter not provided. Calculating steering vector "a" under the ideal array model. ***');
    a = exp(1j*k*rx*cosd(tht')); % the array factor
else
     % If 'aep' is provided, process the data accordingly
    disp('*** "aep" parameter detected. Calculating steering vector "a" from provided data. ***');
    E = aep((num-1)/2+1:(num-1)/2+num,:);
    Eabs = abs(E);
    E = E./repmat(max(Eabs),num,1);
    a = E';                      % the array factor
end

% Calculate A and P
A = zeros(N,N);
for i=1:num
    A = A+a(:,i)*a(:,i)'*sind(tht(i))*pi/num;
end
C = sqrtm(A);
P = [];
for i=1:length(tht)
    P = [P,C^(-1)'*a(:,i)*a(:,i)'*C^(-1)];
end

% Create a figure window for the initial wide beam
figure(1);  

%% Algorithm 1 implementation
while(1) %% Loop C begin 
    half_ML = L_new/2;  % half of the current AMB
    half_TB = 5;        % half of the transition band
    idxML = find(tht<=thetal+half_ML&tht>=thetal-half_ML); % the indicies of the mainlobe region
    idxTB = find((tht<=thetal+half_ML+half_TB+0.001&tht>thetal+half_ML)|(tht>=thetal-half_ML-half_TB-0.001&tht<thetal-half_ML)); % the indicies of the transition band region
    idxSL = setdiff((1:length(tht))',[idxML;idxTB]); % the indicies of the sidelobe region
    
    if(cout_p==0||res==1) % determine whether to enter Loop C
    % Wide beam initialization / minimum mainlobe power gain (MMPG) maximization
        if(cout_p==0)
            xx = ones(N,1)/norm(ones(N,1)); 
        end
        gap_p1 = 1;
        delta = 0.01;
        alpha = 0.20;
        iter_p1 = 0;
        % The first convex optimization problem (inspired by the exsiting PGPS-based algorithm)
        while(alpha>0&&gap_p1>=1e-6)
            cvx_begin quiet
                variable xr(N,1) complex;
                variable G_0;
                minimize -G_0;
                subject to
                    for i=1:length(idxML)
                        real(2*(xx'/C'*a(:,idxML(i))*a(:,idxML(i))'/C*xr)) >= G_0; 
                    end
                    for i=1:length(idxSL)
                        real(2*xx'/C'*a(:,idxSL(i))*a(:,idxSL(i))'/C*xr) <= rho*G_0; 
                    end
                    norm(xr) <= 1;                
            cvx_end
            xx = (0.5+alpha)*xr+(0.5-alpha)*xx;
            alpha = alpha-delta;
            gap_p1 = norm(xx-xr)
            iter_p1 = iter_p1+1
            res = 0;
            L_old = L_new;
        end
        
        % Ensure G_0 > G_min
        if (10*log10(G_0) <= Gain_min)
            L_new = L_new / 2;
            continue;
        end

        % Compute the power gain pattern obtained from the existing PGPS-based algorithm
        w = C\xx; 
        gp = zeros(num,1);
        for i=1:num
            gp(i) = 2*w'*a(:,i)*a(:,i)'*w/(w'*A*w);
        end
        % Plot the power gain pattern of the initial wide beam
        if(cout_p==0) 
            plot(tht, 10*log10(abs(gp)));grid on;hold on;
            xlabel('\theta (degree)');
            ylabel('Power gain (dBi)');    
            xlim([0, 180]);
            ylim([-25, 15]);
            title('Initial Power Gain Pattern');
        end
    end
    
    %% Loop B begin 
    % In actual programming, we integrate Loop B and Loop C under 'while(1)' 
    % in order to share some common codes and reduce redundancy.
    % That's why no 'while' statement for Loop B here.
    L = length(idxTB);
    if(mod(L,2)~=0)
        if(tht(idxTB(end))-thetal > thetal-tht(idxTB(1)))
            idxTB = [idxTB(1);idxTB(end)-1];
        else
            idxTB = [idxTB(1)+1;idxTB(end)];
        end
        L = L+1;
    end
    % Initialization for matrix 'U' and weighting vector 'q'
    R = zeros(round(L/2-1), round(L/2));
    for n = 1:round(L/2-1)
        R(n,n)    = -1;
        R(n, n+1) = 1;
    end
    U = [R, zeros(size(R)); zeros(size(R)), -R];
    q = zeros(L, 1);
    q(1:round(L/2-1)) = 1:round(L/2-1);
    q(round(L/2):L)   = (L-round(L/2)+1):-1:1;
    
    iter_p2 = 0;
    G_0 = eta*G_0;
    if(10*log10(G_0)<=Gain_min) % if G_0 is below G_{min}, then set G_0 = G_{min}
        G_0 = 10^(Gain_min/10);
    end
    
    % Reset iterMax to allow it to be incremented inside Loop A
    currentIterMax = iterMax; 
    
    %% Loop A begin: the second convex optimization problem
    while(iter_p2<currentIterMax)
        cvx_begin quiet
            variable xi(N,1) complex;
            variable t(L,1);
            minimize (q'*t);
            subject to
                t >= 0;
                U*t <= 0;
                for i=1:length(idxML)
                    real(4*xx'*P(:,N*(idxML(i)-1)+1:N*idxML(i))*xi+2*xx'*P(:,N*(idxML(i)-1)+1:N*idxML(i))*xx) >= G_0;
                end
                for i=1:length(idxTB)
                    real(4*xx'*P(:,N*idxTB(i)-N+1:N*idxTB(i))*xi+2*xx'*P(:,N*idxTB(i)-N+1:N*idxTB(i))*xx) + t(i) >= G_0;
                end
                for i=1:length(idxSL)
                    real(4*xx'*P(:,N*idxSL(i)-N+1:N*idxSL(i))*xi+2*xx'*P(:,N*idxSL(i)-N+1:N*idxSL(i))*xx) <= rho*G_0;
                end
                norm(xx+xi) <= 1;
                norm(xi) <= Delta_max;
        cvx_end
        
        xx = xx+xi;
        gap_p2 = norm(xi)
        w = C\xx;
        for i=1:num
            gp(i) = 2*w'*a(:,i)*a(:,i)'*w/(w'*A*w);
        end
        
        % Adjust mainlobe and transition band regions
        idxK = idxTB(find(t>=1e-4));
        Lk = length(idxK);
        if mod(Lk,2) > 0.2
            idxK = [idxK; idxK(ceil(Lk/2))];
        end
        idxML = sort([idxML;setdiff(idxTB, idxK)]);
        idxTB = sort(idxK);
        L = length(idxTB);
        if(gap_p2<=1e-4 && min(10*log10(abs(gp(idxML)))) < 10*log10(G_0)+0.001) 
            iter_p2 = iter_p2 +1;
            break;
        end

        % update matrix U and weighting vector q
        R = zeros(round(L/2-1), round(L/2));
        for n = 1:round(L/2-1)
            R(n,n)    = -1;
            R(n, n+1) = 1;
        end
        U = [R, zeros(size(R)); zeros(size(R)), -R];
        tt = find(t>=1e-4);
        Lt = length(tt);
        if mod(length(tt),2) > 0.2 
            tt = [tt;tt(ceil(Lt/2))];
        end
        tt = sort(tt);
        q = 1./t(tt);
        iter_p2 = iter_p2 + 1
        if(iter_p2==currentIterMax) 
            if(min(10*log10(abs(gp(idxML)))) > 10*log10(G_0)+0.001)
                currentIterMax = currentIterMax + 1; % if not converged
            end
        end
    end
    cout_p = cout_p + 1;
    L_ML_l = thetal-tht(idxML(1));   % left part of mainlobe
    L_ML_r = tht(idxML(end))-thetal; % right part of mainlobe
    L_new = L_ML_l + L_ML_r;
    if(min(10*log10(abs(gp(idxML))))<=Gain_min+0.001) % Loop B ends when G reaches G_{min}
        if(L_new<=L_old && L_ML_l == L_ML_r || cout_end > 3) % to ensure convergence of results
            L_ML = 2*min(L_ML_l,L_ML_r);
            break; % Loop C ends
        else
            res = 1; % re-run Loop C
            cout_end = cout_end + 1;  
        end
    end
end

%% Plot the final power gain pattern
w = C\xx; % the final array excitaion
for i=1:num
    gp(i) = 2*w'*a(:,i)*a(:,i)'*w/(w'*A*w);
end
gp_propose = 10*log10(abs(gp));

figure(2); 
plot(tht, gp_propose);
grid on;
xlim([0, 180]);
ylim([-25, 15]);
xlabel('\theta (degree)');
ylabel('Final Power Gain (dBi)');
title('Final Optimized Power Gain Pattern');

end