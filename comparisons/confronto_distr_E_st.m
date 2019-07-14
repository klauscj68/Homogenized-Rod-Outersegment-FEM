% confronto fra la distribuzione teorica della E_st e quella Monte Carlo
function confronto_distr_E_st


% carica la distribuzione Monte Carlo
salamander=false;
if salamander
    load E_st_salamander
    nome='Salamander';
else
    load E_st_mouse
    nome='Mouse';
end
path('../common',path);

% legge i dati
[R, H, n_sez, flag_geom_sp, dz_0, ...
    n_sd, Z_sd, ...
    n_inc, inc, ...
    taglia, tol_R, tol_angle, ...
    n_ref_cyto, n_ref_id, metodo_cyto, metodo_id, ...
    solver, theta, alpha, tol_fix, norma_inf, peak_delta, ...
    t_fin, n_step_t, downsample, ...
    plot_mesh, plot_num, plot_niagara, plot_pool , inspect, ...
    u_tent, v_tent, tol_stat, ...
    flag_model, flag_model_disc, flag_Ca_clamp, ...
    mode_time, mode_space, n_step_R, equal_step, tau, ...
    nu_RE, nu_RG, cc_R_st, D_R_st, ...
    n_Phi, Phi, random_location, ...
    k_GE, cc_G_st, D_G_st, ...
    k_hyd, PDE_s, k_st, ...
    cc_E_st, D_E_st, k_E, ...
    MC_disk, ...
    n_sample, ...
    epsilon_0, nu, sigma, ...
    cc_u, kk_u, cc_v, kk_v, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex]=data;
% % % 'riduco i campioni a 10000'
% % % n_sample=10000;
% % % massa_E_st{1}=massa_E_st{1}(1:10000,:);
% % % pause

% % passo di integrazione temporale
% t_step=t_fin/n_step_t;
% % istanti di calcolo
% time=(0:n_step_t)'*t_step;
time=time_downsample;

% funzione teorica di densità di probabilità della massa di E*
% modello R*-->E* (flag_model=1)
if n_step_R>1
    error('deve essere n_step_R=1')
end
% tasso di produzione di E_st
if flag_model_disc~=1
    warning('la distribuzione teorica si riferisce alla cascata senza G*')
    nu=nu_RG(1);
else
    nu=nu_RE(1);
end
% costante di decadimento di R_st
tau_R=tau(1);

% tempi su cui la funzione densità è calcolata
n_step_dens=1000;
t_fin_dens=time(end);
% t_fin_dens=1.5;
time_dens=linspace(t_fin_dens/n_step_dens,t_fin_dens,n_step_dens);


% media di E_st, per tutti i tempi time_dens
mu_E_st=zeros(1,n_step_dens);
std_E_st=zeros(1,n_step_dens);
for j=1:n_step_dens
    t_tilda=time_dens(j);
    % valore discriminante della massa di diesterasi
    m_tilda=nu/k_E*(1-exp(-k_E*t_tilda));
    mu_E_st(j) = quadv(@(e)media(e,nu,k_E,t_tilda,tau_R),0,m_tilda) + m_tilda*exp(-t_tilda/tau_R);
    std_E_st(j) = sqrt(quadv(@(e)var(e,nu,k_E,t_tilda,tau_R,mu_E_st(j)),0,m_tilda) + (m_tilda-mu_E_st(j))^2*exp(-t_tilda/tau_R));
end


figure(1)
hold on
title(nome)
plot(time,mean(massa_E_st{1}),'r','Linewidth',2)
plot(time_dens,mu_E_st,'g','Linewidth',2)
xlabel('time [s]')
ylabel('mean(E*) [#]')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 time(end) 0 12])
print -depsc mu

figure(2)
hold on
title(nome)
plot(time,std(massa_E_st{1}),'r','Linewidth',2)
plot(time_dens,std_E_st,'g','Linewidth',2)
xlabel('time [s]')
ylabel('std(E*) [#]')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 time(end) 0 12])
print -depsc std

figure(3)
hold on
title(nome)
plot(time,std(massa_E_st{1})./mean(massa_E_st{1}),'r','Linewidth',2)
plot(time_dens,std_E_st./mu_E_st,'g','Linewidth',2)
xlabel('time [s]')
ylabel('CV(E*) [#]')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 time(end) 0 10])
print -depsc cv


return


function f=media(e,nu,k_E,t_tilda,tau_R)
    f=e*f_distr(e,nu,k_E,t_tilda,tau_R);
return
function f=var(e,nu,k_E,t_tilda,tau_R,mu_e)
    f=(e-mu_e)^2*f_distr(e,nu,k_E,t_tilda,tau_R);
return
function f=f_distr(e,nu,k_E,t_tilda,tau_R)
    f=1/(nu*tau_R)*exp(k_E*t_tilda)*(1+k_E/nu*e*exp(k_E*t_tilda))^(-1/(tau_R*k_E)-1);
return
