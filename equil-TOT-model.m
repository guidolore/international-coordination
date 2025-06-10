%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% "Global Price Shocks and International Monetary Coordination" 
% by Guerrieri, Lorenzoni, Werning (2005)
% 
% Solves for model with Home and Foreign goods a la Gali-Monacelli in Section 6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Parameter definitions                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% production function
eta = 1/6;
aX  = .2;
aL  = 1-aX;

% home bias in consumption
ome = .8;

% labor supply
phi = 2;
Psi = aL^(-phi);

% fraction of firms that can't reset prices
lam = 0.5;

% CES aggregator for non-traded goods varieties
veps = 3;

% oil endowment at initial steady state
X = aX;

% structure containing all parameters, useful to pass into functions
v = struct('eta', eta, ...
    'aX', aX, ...
    'aL', aL, ...
    'ome',ome,...
    'phi', phi, ...
    'Psi', Psi, ...
    'lam', lam, ...
    'veps', veps, ...
    'X', X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Find no coordination steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fun = @(L) ome*(PF(L,v))^(1/eta-1)*(L/aL)^(-1/eta) - Psi*L^phi;
L0 = fzero(fun,1);
C0 = PF(L0,v);
W0 = v.Psi*L0^phi*C0;
PX0 = W0*(L0/aL)^(1/eta);

% store steady state nominal wage in structure v
v.W = W0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Find no coordination equilibrium with oil shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% oil shock
dx = 0.1;
v.X = (1-dx)*aX;    

% find Nash equilibrium using function FOC_Nash
v.mu = 1/ome; % set markup to Nash level
fun = @(L)FOC_Nash(L, v);
L_Nash = fzero(fun,L0);
[PX_Nash, ~, P_Nash, C_Nash, Del_Nash, ~, ~, ~, ~, U_Nash] = equilibrium(L_Nash, v);
%[PX, P_Ht, P_H, C, Del, del, sL, sX, h]

PX = PX_Nash; Ps = P_Nash; CFs = (1-v.ome)*C_Nash;
[P_Ht, P_H, CH, Del, L, X, U] = equilibrium_local(1, v, PX, Ps, CFs);

% vector of possible exchange rates
n_E = 100;
E_vals = linspace(.5,2,n_E);

U_i = zeros(1,n_E);
P_H_i = zeros(1,n_E);

for i = 1:n_E
    PX = PX_Nash; Ps = P_Nash; CFs = (1-v.ome)*C_Nash;
    [~, P_H_i(i), ~, ~, L_i(i), ~, U_i(i)] = equilibrium_local(E_vals(i), v, PX, Ps, CFs);
end

% plot to check optimality 
plot(L_i,U_i,L_Nash,U_Nash,'*')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Comparative statics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_shocks = 100;
shocks = linspace(0,.2,n_shocks);

for i = 1:n_shocks

    v.X = exp(-shocks(i))*aX;    


    % no coordination
    v.mu = 1/ome; v.W = W0;
    fun = @(L)FOC_Nash(L, v);
    L_NC(i) = fzero(fun,L0);
    [PX_NC(i), ~, P_NC(i)] = equilibrium(L_NC(i), v);

    % coordination
    v.mu = 1; v.W = 1; 
    fun = @(L)FOC_Plan(L, v);
    L_C(i) = fzero(fun,1);
    [PX_C(i), ~, P_C(i)] = equilibrium(L_C(i), v);

end
subplot(3,1,1)
plot(-100*shocks,100*log(PX_NC/PX0),-100*shocks,100*log(PX_C),'LineWidth',2)
grid on
xlabel('Oil shock','Interpreter','LaTex','FontSize',16)
ylabel('Oil price','Interpreter','LaTex','FontSize',16)

subplot(3,1,2)
plot(-100*shocks,100*log(P_NC),-100*shocks,100*log(P_C),'LineWidth',2)
grid on
xlabel('Oil shock','Interpreter','LaTex','FontSize',16)
ylabel('Inflation','Interpreter','LaTex','FontSize',16)

subplot(3,1,3)
plot(-100*shocks,100*log(L_NC./L0),-100*shocks,100*log(L_C./aL),'LineWidth',2)
grid on
xlabel('Oil shock','Interpreter','LaTex','FontSize',16)
ylabel('Employment','Interpreter','LaTex','FontSize',16)
legend('No coordination equilibrium','Coordination equilibrium',...
    'Interpreter','LaTex','Location','SouthEast','FontSize',14)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 16, 20]); % or desired size
matlab2tikz('fig-comp-stat-HF.tex')

%fig = gcf;
%set(fig, 'Units', 'normalized');
%set(fig, 'PaperPositionMode', 'auto');
%set(fig, 'InvertHardcopy', 'off'); 
%exportgraphics(fig,'fig6-CS-HF.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Functions definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Production function 
function Y = PF(L,v)
    Y = (v.aL^(1/v.eta)*L.^(1-1/v.eta)+v.aX^(1/v.eta)*v.X^(1-1/v.eta)).^(v.eta/(v.eta-1));
end

% Symmetric equilibrium prices and quantities given L
function [PX, P_Ht, P_H, C, Del, del, sL, sX, h, U] = equilibrium(L, v)
        
    % oil price and good prices
    PX = ((v.aX/v.aL)*(L/v.X))^(1/v.eta)*v.W;
    P_Ht = v.mu*(v.aL*v.W^(1-v.eta) + v.aX*PX^(1-v.eta))^(1/(1-v.eta));
    P_H = (v.lam + (1-v.lam)*P_Ht^(1-v.veps))^(1/(1-v.veps));
    
    % distortion Delta and elasticity of Delta to good price
    Del = P_H^v.veps*(v.lam + (1-v.lam)*P_Ht^(-v.veps));
    del = v.veps*(1-(P_H/P_Ht)*(1/Del));
    
    % consumption
    C = PF(L,v)/Del;

    % labor and oil shares
    sL = v.W*L/(v.W*L+PX*v.X);
    sX = 1-sL;

    % pass-through of exchange rate to home good price
    h = (1-v.lam)*(P_Ht/P_H)^(1-v.veps)*sX;

    % utility
    U = log(C) - v.Psi/(1+v.phi)*L^(1+v.phi);

end

function [P_Ht, P_H, CH, Del, L, X, U] = equilibrium_local(E, v, PX, Ps, CFs)
        
    % good prices
    P_Ht = v.mu*(v.aL*v.W^(1-v.eta) + v.aX*(E*PX)^(1-v.eta))^(1/(1-v.eta));
    P_H = (v.lam + (1-v.lam)*P_Ht^(1-v.veps))^(1/(1-v.veps));
    
    % distortion Delta and elasticity of Delta to good price
    Del = P_H^v.veps*(v.lam + (1-v.lam)*P_Ht^(-v.veps));
    
    % consumption
    CH = v.ome/(1-v.ome)*(P_H/(E*Ps))^(-1)*CFs;
    XLratio = v.aX/v.aL*(E*PX/v.W)^(-v.eta);
    YLratio = (v.aL^(1/v.eta) + v.aX^(1/v.eta)*XLratio^(1-1/v.eta)).^(v.eta/(v.eta-1));
    L = Del*CH/(v.ome*YLratio);
    X = XLratio*L;

    % utility
    U = v.ome*(log(CH)-log(v.ome)) + (1-v.ome)*(log(CFs)-log(1-v.ome)) - v.Psi/(1+v.phi)*L^(1+v.phi) + (1-v.ome)/CFs*PX/Ps*(v.X-X);

end

% Expression in FOC of domestic planner
function DU = FOC_Nash(L, v)
    % get candidate symmetric equilibrium prices and quantities
    [PX, ~, PH, C, ~, del, ~, sX, h] = equilibrium(L, v);  

    % derivatives of domestic equilibrium relation w.r.t. the exchange rate
    pH_prime = h;
    cH_prime = 1-pH_prime;
    l_prime = del*pH_prime + cH_prime + sX*v.eta;
    x_prime = l_prime - v.eta;

    % compute expression in FOC
    DU = v.ome*cH_prime - v.Psi*L^(1+v.phi)*l_prime...
                - 1/C*PX/PH*v.X*x_prime;
end

% Expression in FOC of global planner
function DU = FOC_Plan(L, v)    
    [~, ~, ~, ~, ~, del, sL, ~, h] = equilibrium(L, v);
    % derivatives of global equilibrium relation w.r.t. employment
    pH_prime = h*1/v.eta;
    c_prime = sL - del*pH_prime;

    % compute expression in FOC
    DU = c_prime - v.Psi*L^(1+v.phi);
end

function U = Ut(C,L,v)
    % utility function
    U = log(C) - v.Psi/(1+v.phi)*L^(1+v.phi);
end



