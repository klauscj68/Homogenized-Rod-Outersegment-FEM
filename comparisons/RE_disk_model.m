function RE_disk_model
% statistiche dell'attività di R_st e della massa di E_st al variare di n_step
% modello R*-->E*

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
    j_ex_sat, K_ex]=data;

salamander=false;
if salamander
    t_fin=4;
    nome='salamander';
else
    t_fin=1;
    nome='mouse';
end

% vita media della R*, globale
if equal_step
    tau_tot=n_step_R*tau;
else
    tau_tot=sum(tau);
end
warning(['vita media totale di R* pari a ',num2str(tau_tot)])

% colori
c=['b','g','r','c','m','y','k'];

% tasso di produzione di E_st
if flag_model_disc~=1
    warning('qui è implementata la cascata senza G*')
    nu=nu_RG;
else
    nu=nu_RE;
end

% utilizza nu_1 come nu_max
nu_max=nu(1);
warning('sovrascrivo nu')

% tasso di spegnimento della E_st
k=k_E;

% sottocampiona
% downsample=1;
n_step_t=n_step_t/downsample;

% numero di campioni (metodo Monte Carlo)
% sovrascrive il valore in data.m
n_sample=10000;
if (n_step_t+1)*n_sample*8>200e6
    error('troppi campioni o troppi istanti temporali')
end

% passo di integrazione temporale
t_step=t_fin/n_step_t;
% istanti di calcolo
time=(0:n_step_t)'*t_step;

if (n_step_t+1)*n_sample*8>500e6
    error('troppi campioni')
end

% valori desiderati di n_step_R
valori_n_step_R=[1,2,4];

% inizializza
n_casi=5;

for caso=1:n_casi

    switch caso
        case 1
            ratio_step=1;
            ratio_actv=1;
            nome_caso='\tau_{i+1}=\tau_i; \nu_{i+1}=\nu_i';
        case 2
            ratio_step=1;
            ratio_actv=1/2;
            nome_caso='\tau_{i+1}=\tau_i; \nu_{i+1}=\nu_i/2';
        case 3
            ratio_step=1/2;
            ratio_actv=1;
            nome_caso='\tau_{i+1}=\tau_i/2; \nu_{i+1}=\nu_i';
        case 4
            ratio_step=1/2;
            ratio_actv=1/2;
            nome_caso='\tau_{i+1}=\tau_i/2; \nu_{i+1}=\nu_i/2';
        case 5
            ratio_step=2;
            ratio_actv=1/2;
            nome_caso='\tau_{i+1}=\tau_i*2; \nu_{i+1}=\nu_i/2';
    end

    % METODO MONTE CARLO

    % per differenti valori del numero di passi di spegnimento
    for caso_n_step_R=1:length(valori_n_step_R)
        n_step_R=valori_n_step_R(caso_n_step_R);

        % calcola la durata tau_1 del primo step in modo che la somma delle durate sia tau_tot
        % quindi calcola le durate di tutti i passi in progressione
        % geometrica di passo ratio_step
        tau=(tau_tot/sum(ratio_step.^(0:n_step_R-1)))*(ratio_step.^(0:n_step_R-1));        
        
        % calcola l'attività nu_1 del primo step in modo che la somma delle durate*attività sia tau_tot*nu_max
        % quindi calcola l'attività di tutti i passi in progressione geometrica di passo ratio_step
        nu=(nu_max*sum(ratio_step.^(0:n_step_R-1))/sum((ratio_step*ratio_actv).^(0:n_step_R-1)))*(ratio_actv.^(0:n_step_R-1));
        
        % stampa i valori asintotici per tempi grandi di CV(E_int)
        disp(['n=',num2str(n_step_R)])
        disp(nome_caso)
        disp(['valore asintotico per tempi grandi di CV(E_int)): ',num2str(sqrt(sum((nu.^2).*(tau.^2)))/sum(nu.*tau))])
%         pause

        % aggiunge uno zero in fondo a nu
        % perché lo stato spento non produce R*
        nu=[nu,0];
        
        % durate effettive dei singoli stati: una colonna per ogni esperimento
        t=zeros(n_step_R,n_sample);
        if equal_step
            t=exprnd(tau(1),n_step_R,n_sample);
        else
            for i=1:n_step_R
                t(i,:)=exprnd(tau(i),1,n_sample);
            end
        end
        % tempi di transizione
        t_transition=cumsum(t,1);

        % per ogni istante t_tilda entro time, determina 
        % l'attività  di R e la massa di E
        R=zeros(n_step_t+1,n_sample);
        E=zeros(n_step_t+1,n_sample);
        for passo=1:n_step_t+1
            % istante
            t_tilda=time(passo);
            % calcola il vettore state
            state=sum(time(passo)>t_transition,1)+1;
            % calcola la massa, aggregando i campioni per stato in cui si trova la R*
            % campioni nello stato 1, trattati a parte
            set_sample= (state==1);
            R(passo,set_sample)=nu(1);
            E(passo,set_sample)=nu(1)/k*(1-exp(-k*t_tilda));
            for n=2:n_step_R+1
                % campioni nello stato n
                set_sample= (state==n);
                R(passo,set_sample)=nu(n);
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

        % statistiche di R-actv e E-massa, istante per istamte
        mu_R=mean(R,2);
        sigma_R=std(R,[],2);
        CV_R=sigma_R./mu_R;
        if n_step_R==1
            CV_R_teo=exp(time/tau/2).*sqrt(1-exp(-time/tau));
        end

        R_int=cumtrapz(R,1)*t_step;
        mu_R_int=mean(R_int,2);
        sigma_R_int=std(R_int,[],2);
        CV_R_int=sigma_R_int./mu_R_int;

        mu_E=mean(E,2);
        sigma_E=std(E,[],2);
        CV_E=sigma_E./mu_E;
        if n_step_R==1
% % %             Gamma_E_teo=nu(1)^2*(...
% % %                 +(exp(-time/tau)-exp(-2*k_E*time) )/k_E/(2*k_E-1/tau) ...
% % %                 -(exp(-k_E*time-time/tau)-exp(-2*k_E*time) )/k_E/(k_E-1/tau) ...
% % %                 +(exp(-time/tau)-exp(-k_E*time-time/tau) )/k_E/(k_E-1/tau) ...
% % %                 -(exp(-time/tau)-exp(-2*k_E*time) )/(2*k_E-1/tau)/(k_E-1/tau) ...
% % %                 );
            Gamma_E_teo=2*nu(1)^2*(...
                +exp(-2*k_E*time)/(k_E-1/tau)/(2*k_E-1/tau) ...
                +exp(-time/tau)/k_E/(2*k_E-1/tau) ...
                -exp(-time/tau-k_E*time)/k_E/(k_E-1/tau) ...
                );
            mu_E_teo=nu(1)/(k_E-1/tau)*(exp(-time/tau)-exp(-k_E*time));
            CV_E_teo=sqrt(Gamma_E_teo-mu_E_teo.^2)./mu_E_teo;
        end

        E_int=cumtrapz(E,1)*t_step;
        mu_E_int=mean(E_int,2);
        sigma_E_int=std(E_int,[],2);
        CV_E_int=sigma_E_int./mu_E_int;

% % %         max_R=max(R,[],1);
% % %         mu_R_max=mean(max_R);
% % %         sigma_R_max=std(max_R);
% % %         CV_R_max=sigma_R_max/mu_R_max;
        
% % %         max_E=max(E,[],1);
% % %         mu_E_max=mean(max_E);
% % %         sigma_E_max=std(max_E);
% % %         CV_E_max=sigma_E_max/mu_E_max;


        % PLOT
% % %         figure(1+3*(caso-1))
% % %         plot(time,mu_R/nu_max,[c(caso_n_step_R),'+:'],'Linewidth',1);
% % %         plot(time,mu_R_int/(nu_max*tau_tot),[c(caso_n_step_R),'o:'],'Linewidth',1);
% % %         plot(time,mu_E/(nu_max*tau_tot),[c(caso_n_step_R),'+-'],'Linewidth',1);
% % %         plot(time,mu_E_int/(nu_max*tau_tot^2),[c(caso_n_step_R),'o-'],'Linewidth',1);
% % % 
% % %         figure(2+3*(caso-1))
% % %         plot(time,sigma_R/nu_max,[c(caso_n_step_R),':'],'Linewidth',1);
% % %         plot(time,sigma_E/(nu_max*tau_tot),[c(caso_n_step_R),'-'],'Linewidth',1);

        figure(length(valori_n_step_R)*(caso-1)+caso_n_step_R)
        set(gca,'Fontsize',16)
        set(gca,'Xlim',[0,t_fin])
        box on
        hold on
        title([nome,'; n=',num2str(n_step_R),';  ',nome_caso])
        xlabel('time [s]')
        ylabel('CV')
        axis([0 t_fin 0 1])

        a=plot(time,CV_R,'r.','Linewidth',1.5);
        b=plot(time,CV_R_int,'r--','Linewidth',1.5);
        c=plot(time,CV_E,'g.','Linewidth',1.5);
        d=plot(time,CV_E_int,'g--','Linewidth',1.5);
        if n_step_R==1
            e=plot(time,CV_R_teo,'r-','Linewidth',1.5);
            f=plot(time,CV_E_teo,'g-','Linewidth',1.5);
            [legh,objh]=legend([a,b,c,d,e,f],'R(t)','\int_0^tR(\tau)','E(t)','\int_0^tE(\tau)','Th.R(t)','Th.E(t)','Location','SouthEast');
        else
            [legh,objh]=legend([a,b,c,d],'R(t)','\int_0^tR(\tau)','E(t)','\int_0^tE(\tau)','Location','SouthEast');
        end
        set(legh,'Fontsize',14)

% % %         figure(4)
% % %         for passo=delta_step:delta_step:n_step_t+1
% % %             [F,e]=stairs(sort(E(passo,:)));
% % %             F=F/n_sample;
% % %             plot(e,F,[c(caso_n_step_R),type(caso)],'Linewidth',1)
% % %         end

    end % ciclo su n_step_R

    
end % ciclo sui casi equal / not equal

% figure(1)
% [legh,objh]=legend([aa(1:length(valori_n_step_R)) aa(1:length(valori_n_step_R):end)],['n\_step\_R=' num2str(valori_n_step_R(1))]...
%     ,['n\_step\_R=' num2str(valori_n_step_R(2))],['n\_step\_R=' num2str(valori_n_step_R(3))],'eq\_step, eq\_actv',...
%     'eq\_step,NOT\_eq\_actv','NOT\_eq\_step, eq\_actv','NOT\_eq\_step, NOT\_eq\_actv');
% set(legh,'Fontsize',14)
% 
% figure(2)
% [legh,objh]=legend([bb(1:length(valori_n_step_R)) bb(1:length(valori_n_step_R):end)],['n\_step\_R=' num2str(valori_n_step_R(1))]...
%     ,['n\_step\_R=' num2str(valori_n_step_R(2))],['n\_step\_R=' num2str(valori_n_step_R(3))],'eq\_step, eq\_actv',...
%     'eq\_step,NOT\_eq\_actv','NOT\_eq\_step, eq\_actv','NOT\_eq\_step, NOT\_eq\_actv','Location','South');
% set(legh,'Fontsize',14)
% 
% figure(3)
% [legh,objh]=legend([cc(1:length(valori_n_step_R)) cc(1:length(valori_n_step_R):end)],['n\_step\_R=' num2str(valori_n_step_R(1))]...
%     ,['n\_step\_R=' num2str(valori_n_step_R(2))],['n\_step\_R=' num2str(valori_n_step_R(3))],'eq\_step, eq\_actv',...
%     'eq\_step,NOT\_eq\_actv','NOT\_eq\_step, eq\_actv','NOT\_eq\_step, NOT\_eq\_actv','Location','NorthWest');
% set(legh,'Fontsize',14)


% save('stat','mu_E_max','sigma_E_max','CV_E_max','mu_E_int','sigma_E_int','CV_E_int')

% % % mu_E_max
% % % sigma_E_max
% % % CV_E_max
% % % 
% % % max_teo=nu(1)/k*(1-exp(-k*t_transition));
% % % max=max(E,[],1);
% % % [max_teo;max]
% % % norm(max_teo-max,inf)/norm(max_teo,inf)
% % % 
% % % [mean(max),mean(max_teo),nu(1)*tau(1)/(1+k*tau(1))]
% % % 
% % % [std(max), std(max_teo), nu(1)*tau(1)/(1+k*tau(1))/sqrt(2*tau(1)*k+1)]
% % % 
% % % % 
% % % save xxx
% PLOT
% figure(1)
% print('-depsc',['mu_E_',nome])
%  
% figure(2)
% print('-depsc',['std_',nome])
% 
% figure(3)
% print('-depsc',['cv_',nome])

% % % figure(4)
% % % print('-depsc',['F_distr_',nome])

return


