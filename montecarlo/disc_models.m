% Modelli di disco, Holcman e nostro
function disc_models

close all
clear all

% legge i dati
[ ... % parametri geometrici
    R, H, n_sez, ...
    ... % parametri microstrutturali
    epsilon_0, nu, sigma, ...
    ... % dati sulle incisure
    n_inc, inc, ...
    ... % parametri di discretizzazione
    taglia, tol_R, tol_angle, ...
    n_ref_cyto, n_ref_id, ...
    ... % dati relativi all'integrazione nel tempo
    theta, alpha, tol_fix, norma_inf, ...
    t_fin, n_step_t, t_step, time, ...
    ... % flag di disegno e controllo
    plot_mesh, plot_num, inspect, ...
    ... % dati iniziali di tentativo per lo steady-state
    u_tent, v_tent, tol_stat, ...
    ... % modello 
    flag_model, flag_Ca_clamp, ...
    ... % dati sulla R*
    cc_R_st, D_R_st, ...
    n_step_R, ...
    RK_tot, Rec_tot, M_bind, K_bind, ...
    Ca_val, lambda, lambda_dark, mu, ...
    ... % dati sulla T*
    cc_T_st, D_T_st, nu_RT, k_TE_vol, T_tot, ...
    ... % dati sulla E*
    cc_E_st, D_E_st, PDE, ...
    tau_E_dark, tau_E_adapt, Ca_adapt, k_tau_E2dark, k_tau_E2adapt, ...
    ... % coefficienti di diffusione nel cytosol
    cc_u, kk_u, cc_v, kk_v, ...
    ... % deplition del cGMP: attività catalitica della E*
    k_hyd, k_st, ...
    ... % produzione di cGMP: dati sulla ciclasi
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    ... % dati comuni sui canali del Ca2+
    B_ca, F, ...
    ... % dati sul canale cGMP-operato
    j_cg_max, f_ca, m_cg, K_cg_min, K_cg_max, K_CaM, n_CaM, ...
    ... % dati sullo scambiatore
    j_ex_sat, K_ex, ...
    ... % dati sull'illuminazione di background e flash aggiuntivo
    I_bkgr_ss, I_bkgr_history, Phi_flash_history]=data_nuovo;

% tasso di produzione temporale di E* da parte di T*
k_TE_PDE=k_TE_vol*(2*PDE); % [#/s]

% NOTA BENE: utilizza il vettore lambda_dark
lambda=lambda_dark;

% NOTA BENE: k_E è definito come reciproco di tau_E_dark
k_E=1/tau_E_dark;

% per pulizia cancella i dati non necessari
clear R H n_sez
clear epsilon_0 nu sigma
clear n_inc inc
clear taglia tol_R tol_angle
clear n_ref_cyto n_ref_id
clear theta alpha tol_fix norma_inf
clear plot_mesh plot_num inspect
clear u_tent v_tent tol_stat
clear flag_model flag_Ca_clamp
clear cc_R_st D_R_st
clear RK_tot Rec_tot M_bind K_bind
clear Ca_val lambda_dark
clear cc_T_st D_T_st k_TE_vol T_tot, ...
clear cc_E_st D_E_st PDE
clear tau_E_dark tau_E_adapt Ca_adapt k_tau_E2dark k_tau_E2adapt
clear cc_u kk_u cc_v kk_v
clear k_hyd k_st
clear alpha_max alpha_min m_cyc k_cyc
clear B_ca F
clear j_cg_max f_ca m_cg K_cg_min K_cg_max K_CaM n_CaM
clear j_ex_sat K_ex
clear I_bkgr_ss I_bkgr_history Phi_flash_history

% numero di campioni
n_sample=100000;

% nu_RT=nu_RT*10

% sovracampionamento temporale
% necessario erché di fatto usiamo un'integrazione esplicita
oversample=100;

% passo temporale (sovracampionato)
t_step_over=t_step/oversample;

% numeri di molecole: una riga per ogni campione; una colonna per ogni
% istante temporale
% modello monte carlo (granulare) e continuo
% numeri di particelle in ogni istante
R=false(n_sample,n_step_R,n_step_t+1); % 1 molecola: presente/assente
T_gr=zeros(n_sample,n_step_t+1);
E_gr=zeros(n_sample,n_step_t+1);
T_cont=zeros(n_sample,n_step_t+1);
E_cont=zeros(n_sample,n_step_t+1);

% inizializza
R_fin=false(n_sample,n_step_R);
T_gr_fin=zeros(n_sample,1);
E_gr_fin=zeros(n_sample,1);
T_cont_fin=zeros(n_sample,1);
E_cont_fin=zeros(n_sample,1);

% cutoff
exp_mu_dt=exp(-mu*t_step_over);
exp_lambda_dt=exp(-lambda*t_step_over);
exp_k_TE_PDE_dt=exp(-k_TE_PDE*t_step_over);
exp_k_E_dt=exp(-k_E*t_step_over);

% dato iniziale: una R* per ogni campione, in stato 1
R(:,1,1)=true;

% ciclo sui tempi
for ind_t=1:n_step_t
    disp([num2str(ind_t),'/',num2str(n_step_t)])

    % ricopia dal passo precedente
    R_in=R(:,:,ind_t);
    T_gr_in=T_gr(:,ind_t);
    E_gr_in=E_gr(:,ind_t);
    T_cont_in=T_cont(:,ind_t);
    E_cont_in=E_cont(:,ind_t);

    % ciclo su tutti gli istanti di sovracampionamento
    for ov=1:oversample

        % ricopia dal passo oversample precedente
        R_fin=R_in;
        T_gr_fin=T_gr_in;
        E_gr_fin=E_gr_in;
        T_cont_fin=T_cont_in;
        E_cont_fin=E_cont_in;

        % ciclo su tutti gli stati di R*
        for s=1:n_step_R
            % campioni con R* che sono in stato s al tempo t
            stato_s=find(R_in(:,s));
            % loro numero
            n_stato_s=length(stato_s);
            % produzione di T*_gr
            media=nu_RT(s)*t_step_over;
            % NOTA: è la distribuzione poissoniana corretta per la generazione di T* ?
            T_gr_fin(stato_s)=T_gr_fin(stato_s)+poissrnd(media,n_stato_s,1);
            % produzione di T*_cont
            T_cont_fin(stato_s)=T_cont_fin(stato_s)+media;

            % verifica se la R* deve spegnersi
            shut_rand=rand(n_stato_s,1);
            ind_shutoff=shut_rand>exp_mu_dt(s);
            ind_alive=shut_rand<=exp_mu_dt(s);
            R_fin(stato_s(ind_shutoff),s)=false;

            % limitatamente a quelle che non si sono spente, verifica se devono transitare in stato s+1
            % è corretto considerare i transiti solo su quelle che non si sono spente?
            ind_transit=rand(n_stato_s,1)>exp_lambda_dt(s);
            R_fin(stato_s(ind_alive&ind_transit),s)=false;
            if s<n_step_R
                R_fin(stato_s(ind_alive&ind_transit),s+1)=true;
            end

        end
        
        % legame T*E
        % granulare
        % ciascuna T_gr_in ha una probabilità exp(-k_TE_PDE*t_step_over) di non aver reagito a fine step
% % %         % usa la distribuzione esponenziale
% % %         num=zeros(n_sample,1);
% % %         for samp=1:n_sample
% % %             if T_gr_in(samp)>0
% % %                 reag_rand=rand(T_gr_in(samp),1);
% % %                 num(samp)=length(find(reag_rand>exp_k_TE_PDE_dt));
% % %                 
% % %             end
% % %         end
        % usa la distribuzione binomiale
        % NOTA: così la media è T_gr_in * (1-exp(-k_TE_PDE*t_step_over)) * exp(-k_TE_PDE*t_step_over)
        % e non è esattamente uguale a media_cont definita sotto;
        % sono tanto più vicini quanto più t_step_over è piccolo (oversample grande): discutere!
        % altra alternativa sarebbe: num=binornd(T_gr_in,k_TE_PDE*t_step_over);
        num=binornd(T_gr_in,1-exp_k_TE_PDE_dt);
        T_gr_fin=T_gr_fin-num;
        E_gr_fin=E_gr_fin+num;
        % continuo
        media_cont=k_TE_PDE*T_cont_in*t_step_over;
        T_cont_fin=T_cont_fin-media_cont;
        E_cont_fin=E_cont_fin+media_cont;

        % disattivazione di E*
        % granulare
        % ciascuna E_gr_in ha una probabilità exp(-k_E*t_step_over) di essere disattivata a fine step
% % %         % usa la distribuzione esponenziale
% % %         num=zeros(n_sample,1);
% % %         for samp=1:n_sample
% % %             if E_gr_in(samp)>0
% % %                 reag_rand=rand(E_gr_in(samp),1);
% % %                 num(samp)=length(find(reag_rand>exp_k_E_dt));
% % %             end
% % %         end
        % usa la distribuzione binomiale
        % NOTA: così la media è E_gr_in * (1-exp(-k_E*t_step_over)) * exp(-k_E*t_step_over)
        % e non è esattamente uguale a media_cont definita sotto;
        % sono tanto più vicini quanto più t_step_over è piccolo (oversample grande): discutere!
        % altra alternativa sarebbe: num=binornd(E_gr_in,k_E*t_step_over);
        num=binornd(E_gr_in,1-exp_k_E_dt);
        E_gr_fin=E_gr_fin-num;
        % continuo
        media_cont=k_E*E_cont_in*t_step_over;
        E_cont_fin=E_cont_fin-media_cont;

        % prepara per il passo oversample successivo
        R_in=R_fin;
        T_gr_in=T_gr_fin;
        E_gr_in=E_gr_fin;
        T_cont_in=T_cont_fin;
        E_cont_in=E_cont_fin;

    end
    
    % salva
    R(:,:,ind_t+1)=R_fin;
    T_gr(:,ind_t+1)=T_gr_fin;
    E_gr(:,ind_t+1)=E_gr_fin;
    T_cont(:,ind_t+1)=T_cont_fin;
    E_cont(:,ind_t+1)=E_cont_fin;

end



save('risultati','R','T_gr','E_gr','T_cont','E_cont')

% disegna

figure()
hold on
box on
c='bgrcmykbgrcmykbgrcmyk';
for s=1:n_step_R
    Y=reshape(R(:,s,:),n_sample,n_step_t+1);
    plot(time,mean(Y),c(s),'Linewidth',2)
end
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel('mean(R*) [#]');
set(aa,'fontsize',18);
set(gca,'Xlim',[0 t_fin])
print('-depsc','mean_R')


figure()
hold on
box on
c='bgrcmykbgrcmykbgrcmyk';
for s=1:n_step_R
    Y=reshape(R(:,s,:),n_sample,n_step_t+1);
    plot(time,std(Y),c(s),'Linewidth',2)
end
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel('std(R*) [#]');
set(aa,'fontsize',18);
set(gca,'Xlim',[0 t_fin])
print('-depsc','std_R')


figure()
hold on
box on
c='bgrcmykbgrcmykbgrcmyk';
for s=1:n_step_R
    Y=reshape(R(:,s,:),n_sample,n_step_t+1);
    plot(time,std(Y)./mean(Y),c(s),'Linewidth',2)
end
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel('CV(R*)');
set(aa,'fontsize',18);
set(gca,'Xlim',[0 t_fin])
print('-depsc','CV_R')


grafici(time,T_cont,T_gr,'T')

grafici(time,E_cont,E_gr,'E')

return


function grafici(time,Y_cont,Y_gr,nome)
% disegna
figure()
hold on
box on
a=plot(time,mean(Y_cont),'b','Linewidth',2);
b=plot(time,mean(Y_gr),'r','Linewidth',2);
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel(['mean(',nome,'*) [#]']);
set(aa,'fontsize',18);
[legh,objh]=legend([a b],'Ours','Holcman''s','location','NorthEast');
set(legh,'fontsize',18);
set(gca,'Xlim',[min(time) max(time)])
print('-depsc',['mean_',nome])

figure()
hold on
box on
a=plot(time,std(Y_cont),'b','Linewidth',2);
b=plot(time,std(Y_gr),'r','Linewidth',2);
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel(['std(',nome,'*) [#]']);
set(aa,'fontsize',18);
[legh,objh]=legend([a b],'Ours','Holcman''s','location','NorthEast');
set(legh,'fontsize',18);
set(gca,'Xlim',[min(time) max(time)])
print('-depsc',['std_',nome])

figure()
hold on
box on
a=plot(time,std(Y_cont)./mean(Y_cont),'b','Linewidth',2);
b=plot(time,std(Y_gr)./mean(Y_gr),'r','Linewidth',2);
set(gca,'fontsize',18);
aa=xlabel('time [s]');
set(aa,'fontsize',18);
aa=ylabel(['CV(',nome,'*)']);
set(aa,'fontsize',18);
[legh,objh]=legend([a b],'Ours','Holcman''s','location','SouthEast');
set(legh,'fontsize',18);
set(gca,'Xlim',[min(time) max(time)])
print('-depsc',['CV_',nome])

return
