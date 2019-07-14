% soluzione numerica del sistema non lineare che descrive la soluzione stazionaria
function [u_ss,v_ss]=steady_state(R, H, nu, epsilon_0, ...
    k_hyd, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex, ...
    u_tent, v_tent, tol_stat)

fprintf('\nSoluzione dello steady-state\n');

% superficie laterale
Sigma_rod=2*pi*R*H;

% dato iniziale
x0(1)=u_tent;
x0(2)=v_tent;

% risolve
x = fsolve(@fun, x0, optimset('TolFun', tol_stat), epsilon_0, nu, ...
    k_hyd, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex, Sigma_rod); 

% mette la soluzione in u_ss e v_ss
u_ss=x(1);
v_ss=x(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun(x, epsilon_0, nu, ...
    k_hyd, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex, Sigma_rod)
% estrae le incognite u,v da x
u=x(1);
v=x(2);

% output
y=zeros(1,2);

% prima equazione:
% bilancio dei flussi del cGMP sulle facce dei dischi
y(1)=k_hyd*PDE_s*u-(alpha_min+(alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v^m_cyc))*(1/2*nu*epsilon_0);



% seconda equazione:
% bilancio dei flussi del Ca2+ sulla membrana plasmatica
y(2)=j_ex_sat/(Sigma_rod*B_ca*F)*v/(K_ex+v)-j_cg_max*f_ca/(Sigma_rod*B_ca*F*2)*u^m_cg/(K_cg^m_cg+u^m_cg);
