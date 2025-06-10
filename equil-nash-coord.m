%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Finds equilibria of example in "Global Price Shocks and International Monetary
% Coordination" by Guerrieri, Lorenzoni, Werning (2005)
% 
% Compares uncoordinated and coordinated equilibria
% Plots comparative statics for oil shock
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Parameter definitions                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% production function
eta = 1/6;
ax  = .2;
al  = 1- ax;

% labor supply 
phi = 2;
Psi = al^(-phi);

% fraction of firms that can't reset prices
lam = .5;

% CES aggregator for non-traded goods varieties
veps = 3;

% oil endowment 
x = ax;

% structure containing all parameters, useful to pass into functions
v = struct('eta', eta, ...
    'ax', ax, ...
    'al', al, ...
    'phi', phi, ...
    'Psi', Psi, ...
    'lam', lam, ...
    'veps', veps, ...
    'x', x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Compare private and planner optimization                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% oil shock
dx = .05;
v.X = exp(-dx)*ax;    

% vector of possible output levels Y
n_y = 100;
y_vals = linspace(-.032,-.02,n_y);

% find Nash equilibrium using function FOC_Nash
fun = @(Y)FOC_Nash(Y, v);
Y_Nash = fzero(fun, .9);
[~, ~, ~, P_N_Nash, ~, Q_Nash] = equilibrium(Y_Nash, v);    
u_Nash = U_i(Y_Nash,Q_Nash,v);

% find planner optimum using function FOC_Plan
fun = @(Y)FOC_Plan(Y, v);
Y_Plan = fzero(fun, .9);
[u_Plan, P_N_Plan, Q_Plan] = U_p(Y_Plan,v);

u_i = zeros(1,n_y);
u_p = zeros(1,n_y);
P_N_i = zeros(1,n_y);
P_N_2 = zeros(1,n_y);
P_N_p = zeros(1,n_y);
Q_p = zeros(1,n_y);

for i = 1:n_y
    [u_i(i), P_N_i(i)] = U_i(exp(y_vals(i)),Q_Nash,v);
    [u_p(i), P_N_p(i), Q_p(i)] = U_p(exp(y_vals(i)),v);
    [~, P_N_2(i)] = U_i(exp(y_vals(i)),Q_Plan,v);
end

% plot Phillips curve for individual country and world
figure(1)
plot(100*y_vals,100*log(P_N_i),100*y_vals,100*log(P_N_p),100*y_vals,100*log(P_N_2),...
    100*log(Y_Nash),100*log(P_N_Nash),'ob',100*log(Y_Plan),100*log(P_N_Plan),'*r','LineWidth',2)
grid on
hold on
axis("tight")
xlabel('Output $$y$$','Interpreter','LaTex','FontSize',16)
ylabel('Inflation $$\pi_N$$','Interpreter','LaTex','FontSize',16)
grid on
legend('Local Phillips curve, high oil price','World Phillips curve',...
        'Local Phillips curve, low oil price','Location','southeast')

%%%%% lines to create tikz figure

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 14, 14]); % or desired size
% set(gca, 'FontSize', 11);
% matlab2tikz('fig-PC.tex')

figure(2)
plot(100*y_vals,100*(exp(u_i-u_Plan)-1),100*y_vals,100*(exp(u_p-u_Plan)-1),...
    100*log(Y_Nash),100*(exp(u_Nash-u_Plan)-1),'o',100*log(Y_Plan),0,'o','LineWidth',2)
axis("tight")
xlabel('Output $$y$$','Interpreter','LaTex','FontSize',16)
ylabel('Utility','Interpreter','LaTex','FontSize',16)

%%%%% lines to create tikz figure

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 10]); % or desired size
% matlab2tikz('fig-U.tex')

% compute decomposition of FOC

f = zeros(1,n_y); f1 = zeros(1,n_y); f2 = zeros(1,n_y); f3 = zeros(1,n_y)

for i = 1:n_y
    f(i) = FOC_Nash(exp(y_vals(i)), v);
    f1(i) = FOC_Nash1(exp(y_vals(i)), v);
    f2(i) = FOC_Nash2(exp(y_vals(i)), v);
    f3(i) = FOC_Nash3(exp(y_vals(i)), v);
    f4(i) = FOC_Nash4(exp(y_vals(i)), v);
end

figure(3)
subplot(1,2,1)
plot(100*y_vals,f,...
    100*y_vals,f1,100*y_vals,f2,100*y_vals,f3,...
    100*log(Y_Nash),0,'ob','LineWidth',2)
axis([-3.2 -2 -0.08 0.06])
grid on
title('Country','Interpreter','LaTex','FontSize',16)
xlabel('Output $$y$$','Interpreter','LaTex','FontSize',16)
legend('Marginal benefit','Labor wedge',...
    'Inflation distortion','Oil wedge',...
    'Interpreter','LaTex','FontSize',12,'Location','SouthEast')

subplot(1,2,2)
plot(100*y_vals,f1+f2+f4,...
    100*y_vals,f1,100*y_vals,f2+f4,...,
    100*log(Y_Plan),0,'ob','LineWidth',2)
axis([-3.2 -2 -0.08 0.06])
grid on
title('World','Interpreter','LaTex','FontSize',16)
xlabel('Output $$y$$','Interpreter','LaTex','FontSize',16)
legend('Marginal benefit','Labor wedge',...
    'Inflation distortion',...
    'Interpreter','LaTex','FontSize',12,'Location','SouthEast')
legend('Location', 'northeast')

%%%%% lines to create tikz figure

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 26, 12]); % or desired size
% matlab2tikz('fig-welfare.tex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Comparative statics with oil shock 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid for different sizes of oil shocks dx
num_vals = 100;
dx_max = .05;
dx_vals = linspace(-0.01, dx_max, num_vals);

% vectors of optimal responses for the different oil shocks
Y_Nash = zeros(1, num_vals);
Y_Plan = zeros(1, num_vals);

% calculate best response for each value of oil shock
for i = 1:num_vals
    v.X = exp(-dx_vals(i))*ax;    
    % calculate Nash equilibrium output
    fun = @(y)FOC_Nash(y, v);
    Y_Nash(i) = fzero(fun, .9);
    % calculate planner output
    fun = @(y)FOC_Plan(y, v);
    Y_Plan(i) = fzero(fun, .9);
end

figure(4)
plot(-100*dx_vals,100*log(Y_Nash),-100*dx_vals,100*log(Y_Plan),'LineWidth',2)
legend('No coordination','Coordination',...
    'Interpreter','LaTex','Location','SouthEast','FontSize',14)
grid on
hold on
xlabel('Oil shock $$\bar{x}$$','Interpreter','LaTex','FontSize',16)
ylabel('Output $$y$$','Interpreter','LaTex','FontSize',16)

%%%%% lines to create tikz figure

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 10]); % or desired size
% matlab2tikz('fig-comp-stat.tex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Function definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% private utility
function [U, pN] = U_i(Y,Q,v)
    fun = @(p)f_PN(p, Y, Q, v);
    pN = fzero(fun, 1);
    pT = Y*pN;
    px = Q*pT;
    pI = (v.al + v.ax*px^(1-v.eta))^(1/(1-v.eta));      
    xl = v.ax/v.al*px^(-v.eta);
    Del = (pN.^v.veps).*(v.lam+(1-v.lam)*pI.^(-v.veps));
    factor = (v.al^(1/v.eta) + v.ax^(1/v.eta)*xl^(1-1/v.eta))^(v.eta/(v.eta-1));        
    l = Del*Y/factor;
    x = xl*l;
    U = log(Y) - v.Psi/(1+v.phi)*l^(1+v.phi) - Q*(x-v.X);
end

% function to find pN
function F = f_PN(P_N, Y, Q, v)
    P_T = Y*P_N;
    P_X = Q*P_T;
    P_Nt = (v.al + v.ax*P_X^(1-v.eta))^(1/(1-v.eta));      
    P_Ni = (v.lam + (1-v.lam)*P_Nt^(1-v.veps))^(1/(1-v.veps));     
    F = P_Ni-P_N;
end

% planner utility
function [U, P_N, Q] = U_p(Y,v)
    [L, ~, ~, P_N, ~, Q] = equilibrium(Y, v);    
    U = log(Y) - v.Psi/(1+v.phi)*L^(1+v.phi);
end

% function for FOC of Nash equilibrium
function DU = FOC_Nash(Y, v)
    [L, ~, ~, ~, ~, Q, ~, sL, sX, h, del] = equilibrium(Y, v);    
    DU = 1/Y - v.al^(-v.phi)*L^v.phi*(1+(v.eta*sX+del*h)/(1-h))*(L/Y)...
                - Q*(1-(v.eta*sL-del*h)/(1-h))*v.X/Y;
end

% function for FOC of planner
function DU = FOC_Plan(Y, v)
    [L, ~, ~, ~, ~, ~, ~, sl, sx, h, del] = equilibrium(Y, v);    
    % DU is expression in FOC, best response where DU == 0
    DU = 1/Y - v.al^(-v.phi)*L^v.phi*(1+(v.eta*sx+del*h)/(1-h))*(L/Y)...
        - v.al^(-v.phi)*L^(1+v.phi)*(v.eta*sx+del*h)/(v.eta*sl-del*h)/v.X...
        *(1-(v.eta*sl-del*h)/(1-h))*v.X/Y;
end

% functions for Nash FOC decomposition
function DU = FOC_Nash1(Y, v)
    [L, ~, ~, ~, ~, ~, Del, sL, ~, ~, ~] = equilibrium(Y, v);    
    DU = 1/Y - 1/sL*v.al^(-v.phi)*L^v.phi*(L/Y);
end
function DU = FOC_Nash2(Y, v)
    [L, ~, ~, ~, ~, ~, Del, sL, ~, h, del] = equilibrium(Y, v);    
    DU = -(L/Y)*(1/sL)*v.al^(-v.phi)*L^v.phi*del*h/(1-h);
end
function DU = FOC_Nash3(Y, v)
    [L, ~, ~, ~, ~, Q, ~, sL, sX, h, del] = equilibrium(Y, v);    
    DU = -1/Y*(-sX/sL*v.al^(-v.phi)*L^(1+v.phi)+Q*v.X)*(1-(v.eta*sL-del*h)/(1-h));
end
function DU = FOC_Nash4(Y, v)
    [L, ~, ~, ~, ~, ~, ~, sL, ~, h, del] = equilibrium(Y, v);  
    DU = -(L/Y)*(1/sL)*v.al^(-v.phi)*L^v.phi*del*h/(1-h)*((1-h)/(v.eta*sL-del*h)-1);
end

% symmetric equilibrium prices and quantities given output Y
% and given oil shock X 
function [L, P_X, P_Nt, P_N, P_T, Q, Del, sL, sX, h, del] = equilibrium(Y, v)
        
    % iterative algorithm based on guessing the distortion Delta,
    % solving all other variables, updating Delta, checking convergence

    % start with guess Delta = 1
    Del = 1;
    for i = 1:100    
        % labor needed to produce output y, given distortion delta    
        L = (((Del*Y)^(1-1/v.eta)-v.ax^(1/v.eta)*v.X^(1-1/v.eta))/(v.al^(1/v.eta)))^(v.eta/(v.eta-1));        
        % prices
        P_X = (v.ax*L/(v.al*v.X))^(1/v.eta);        
        P_Nt = (v.al+v.ax*P_X^(1-v.eta))^(1/(1-v.eta));      
        P_N = (v.lam + (1-v.lam)*P_Nt^(1-v.veps))^(1/(1-v.veps));     
        P_T = Y*P_N;        
        Q = P_X./P_T;
        % update Delta and check convergence
        Deli = (P_N.^v.veps).*(v.lam+(1-v.lam)*P_Nt.^(-v.veps));
        if max(abs(Deli-Del))<1e-10
            break
        else
            Del = Deli;
        end
        if i == 100
            disp("ERROR: fixed point did not converge");
        end
    end
    % compute shares    
    temp = (v.al/v.ax)^(1/v.eta)*(L/v.X).^(1-1/v.eta);
    sL = temp/(1+temp);
    sX = 1-sL;

    % pass-through of oil price into non-tradable price p_N
    h = (1-v.lam*P_N^(v.veps-1))*sX;
    
    % Delta'(p_N): derivative of the distortion function w.r.t. pN
    DDel = v.veps*Del/P_N - v.veps*((1-v.lam)^(1/(1-v.veps)))*(P_N^(1-v.veps)-v.lam)^(-1/(1-v.veps));
    
    % delta: elasticity of distortion w.r.t. pN
    del = DDel/Del*P_N;

end
