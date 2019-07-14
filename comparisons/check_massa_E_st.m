function check_massa_E_st
% verifica delle funzioni che danno la distribuzione statistica della massa di E*
% modello R*-->E*
% uno o più step di disattivazione

close all

% aggiunge alla path la directory id
path('../common',path);

% parametri su cui è fatto il check

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
    j_ex_sat, K_ex]=data_n_step2;

% tasso di produzione di E_st
if flag_model_disc~=1
    warning('la distribuzione teorica si riferisce alla cascata senza G*')
    nu=nu_RG;
else
    nu=nu_RE;
end
% aggiunge uno zero in fondo a nu, per uniformare le formule
% significa che lo stato spento di R* non produce E*
nu=[nu,0];

% tasso di spegnimento della E_st
k=k_E;

% sottocampiona
n_step_t=n_step_t/downsample;

% numero di campioni (metodo Monte Carlo)
% sovrascrive il valore in data.m
n_sample=50000;
if (n_step_t+1)*n_sample*8>200e6
    error('troppi campioni o troppi istanti temporali')
end

% passo di integrazione temporale
t_step=t_fin/n_step_t;
% istanti di calcolo
time=(0:n_step_t)'*t_step;

% DISTRIBUZIONE TEORICA
mu_teo=zeros(n_step_t+1,1);
sigma_teo=zeros(n_step_t+1,1);
if n_step_R==1
    % caso n_step_R=1
    tau_R=tau(1);
    % media di E_st, per tutti i tempi
    for j=1:n_step_t+1
        t_tilda=time(j);
        % evita il calcolo in t=0
        if t_tilda<1000*eps
            t_tilda=1000*eps;
        end
        % valore discriminante della massa di E*
        m_tilda=nu(1)/k*(1-exp(-k*t_tilda));
        % media e scarto quadratico medio, come integrali di e*f_E(e) ed (e-mu)^2*f_E(E)
        % la f_E contiene una delta: quindi c'è un pezzo aggiuntivo, oltre l'integrale
        mu_teo(j) = quadv(@(e)media(e,nu(1),k,t_tilda,tau_R),0,m_tilda) + m_tilda*exp(-t_tilda/tau_R);
        sigma_teo(j) = sqrt(quadv(@(e)var(e,nu(1),k,t_tilda,tau_R,mu_teo(j)),0,m_tilda) + (m_tilda-mu_teo(j))^2*exp(-t_tilda/tau_R));
    end
else
    warning('caso teorico n>1 non implementato')
end
cv_teo=sigma_teo./mu_teo;


% METODO MONTE CARLO
% durate effettive dei singoli stati: una colonna per ogni esperimento
t=zeros(n_step_R,n_sample);
if equal_step
    t=exprnd(tau,n_step_R,n_sample);
else
    for i=1:n_step_R
        t(i,:)=exprnd(tau(i),1,n_sample);
    end
end
% tempi di transizione
t_transition=cumsum(t,1);

% per ogni istante t_tilda entro time, determina la massa E(t_tilda)
E=zeros(n_step_t+1,n_sample);
for passo=1:n_step_t+1
    % istante
    t_tilda=time(passo);
    % calcola il vettore state
    state=sum(t_tilda>t_transition,1)+1;
    % calcola la massa, aggregando i campioni per stato in cui si trova la R*
    % campioni nello stato 1, trattati a parte
    set_sample= (state==1);
    E(passo,set_sample)=nu(1)/k*(1-exp(-k*t_tilda));
    for n=2:n_step_R+1
        % campioni nello stato n
        set_sample= (state==n);
        % primo termine dell'espressione di E(t)
        E(passo,set_sample)=nu(1)/k*(exp(k*(t_transition(1,set_sample)-t_tilda))-exp(-k*t_tilda));
        % termini da 2 a n-1 dell'espressione di E(t)
        for j=2:n-1
            E(passo,set_sample)=E(passo,set_sample)+ ...
                nu(j)/k*(exp(k*(t_transition(j,set_sample)-t_tilda))-exp(k*(t_transition(j-1,set_sample)-t_tilda)));
        end
        % ultimo termine dell'espressione di E(t)
        E(passo,set_sample)=E(passo,set_sample)+nu(n)/k*(1-exp(k*(t_transition(n-1,set_sample)-t_tilda)));
    end
end

% statistiche, istante per istamte
mu_MC=mean(E')';
sigma_MC=std(E')';
cv_MC=sigma_MC./mu_MC;


% PLOT
figure(1)
hold on
% title('Salamander')
plot(time,mu_MC,'r','Linewidth',2)
plot(time,mu_teo,'g','Linewidth',2)
xlabel('time [s]')
ylabel('mean(E*) [#]')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
print -depsc mu

figure(2)
hold on
% title('Salamander')
plot(time,sigma_MC,'r','Linewidth',2)
plot(time,sigma_teo,'g','Linewidth',2)
xlabel('time [s]')
ylabel('std(E*) [#]')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
print -depsc std

figure(3)
hold on
% title('Salamander')
plot(time,cv_MC,'r','Linewidth',2)
plot(time,cv_teo,'g','Linewidth',2)
xlabel('time [s]')
ylabel('CV(E*)')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
print -depsc cv

figure(4)
hold on
% title('Salamander')
e_treshold=nu(1)/k;
n_punti=500;
e_teo= e_treshold/n_punti:e_treshold/n_punti:1.2*e_treshold;
n_curve=5;
delta_step=fix((n_step_t+1)/n_curve);
for passo=delta_step:delta_step:n_step_t+1
    [F_MC,e_MC]=stairs(sort(E(passo,:)));
    F_MC=F_MC/n_sample;
    plot(e_MC,F_MC,'r','Linewidth',2)
    if n_step_R==1
        F_teo=F_distr(e_teo,nu(1),k,time(passo),tau_R);
        plot(e_teo,F_teo,'g','Linewidth',2)
    end
end
xlabel('mass E* [#]')
ylabel('cumulative probability, for different times')
axis([0 max(e_teo) 0 1])
% set(gca,'XScale','log')
% set(gca,'YScale','log')
print -depsc F_distr


return



function f=media(e,nu,k,t_tilda,tau_R)
    f=e*f_distr(e,nu,k,t_tilda,tau_R);
return

function f=var(e,nu,k,t_tilda,tau_R,mu_e)
    f=(e-mu_e)^2*f_distr(e,nu,k,t_tilda,tau_R);
return

function f=f_distr(e,nu,k,t_tilda,tau_R)
    f=1/(nu*tau_R)*exp(k*t_tilda)*(1+k/nu*e*exp(k*t_tilda))^(-1/(tau_R*k)-1);
return

function f=F_distr(e_val,nu,k,t_tilda,tau_R)
% nota: qui e_val è un vettore
n_val=length(e_val);
f=zeros(size(e_val));
e_treshold=nu/k*(1-exp(-k*t_tilda));
for j=1:n_val
    if e_val(j)<=e_treshold
        f(j) = quadv(@(e)f_distr(e,nu,k,t_tilda,tau_R),0,e_val(j));
    else
        f(j)=1;
    end
end
return
