% Modelli TWS e GWS
function cyto_model

% aggiunge alla path la directory id
path('../common',path);
path('../id',path);
path('../cyto',path);

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
    mode_time, mode_space, n_step_R, mu, lambda, ...
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

% metodo utilizzato
% flag_theta=metodo_cyto(1)==1;
% flag_theta=false
flag_theta=true;

% nome del file di ouput
nome=[];

% tipo di soluzione
% flag_gws=(flag_model==3);
% flag_gws=true;
flag_gws=false;
if flag_gws
    nome=[nome,'t'];
else
    nome=[nome,'f'];
end

% modello nonlineare (true) o linearizzato (false)
flag_nonlin=true;
% flag_nonlin=false;
if flag_nonlin
    nome=[nome,'t'];
else
    nome=[nome,'f'];
end

% termine saturativo assente/presente
% cioè u costante al valore dark/pari al suo valore effettivo
% nelò termine di attivazione
saturation_effect=true;
% saturation_effect=false;
if saturation_effect
    nome=[nome,'t'];
else
    nome=[nome,'f'];
end

% metodologia di calcolo della corrente (linearizzata/non linearizzata)
flag_curr_nonlin=true;
% flag_curr_nonlin=false;
if flag_curr_nonlin
    nome=[nome,'t'];
else
    nome=[nome,'f'];
end

% Calcium clamp
 flag_Ca_clamp=false;
%flag_Ca_clamp=true;
if flag_Ca_clamp
    nome=[nome,'t'];
else
    nome=[nome,'f'];
end

if flag_Ca_clamp && ~flag_theta
    error('il calcium clamp non è programmato col metodo ode')
end

% numero di intervalli in cui è suddiviso il segmento (0,H/2)
nz=50;
% intervallo di discretizzazione (1/2 dell'altezza)
dz=(H/2)/nz;
if flag_gws
    % zero per il GWS
    nz=0;
end

% annullamento dei flussi al bordo per il calcolo della soluzione stazionaria
% steady-state solution computed by imposing vanishing boundary fluxes
[u_ss,v_ss]=steady_state(R, H, nu, epsilon_0, ...
    k_hyd, PDE_s, alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, j_ex_sat, K_ex, ...
    u_tent, v_tent, tol_stat);

% calcolo dei coefficienti della matrice linearizzata
theta_0=1/(1+nu);
Sigma_rod=2*pi*R*H;
V_cyto=(1-theta_0)*pi*R^2*H;
a=zeros(2,2);
% derivata del termine a destra dell'equazione per cG rispetto a cG
a(1,1) = -k_hyd*(PDE_s/(1/2*nu*epsilon_0));
% derivata del termine a destra dell'equazione per cG rispetto a Ca
a(1,2) = -(alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v_ss^m_cyc)^2*m_cyc*v_ss^(m_cyc-1);
% derivata del termine a destra dell'equazione per Ca rispetto a cG
a(2,1) = j_cg_max*f_ca/(2*B_ca*F*V_cyto)*m_cg*K_cg^m_cg*u_ss^(m_cg-1)/(K_cg^m_cg+u_ss^m_cg)^2;
% derivata del termine a destra dell'equazione per Ca rispetto a Ca
a(2,2) = -j_ex_sat/(B_ca*F*V_cyto)*K_ex/(K_ex+v_ss)^2;

% autovalori e autovettori
disp('matrice dei coefficienti del sistema linearizzato')
disp(a)
[V,D]=eig(a);
disp('autovettori della matrice')
disp(diag(D))
disp('autovalori della matrice')
disp(V)

% pulsazione caratteristica e coefficiente di smorzamento
omega=sqrt(det(a));
zeta=-trace(a)/(2*omega);
disp('pulsazione propria')
disp(omega)
disp('smorzamento')
disp(zeta)


return

% carica la distribuzione Monte Carlo di massa di E_st
load E_st
n_sample_rid=10000;
disp(['riduco i campioni a ',num2str(n_sample_rid)])
disp('press any key')
n_sample=n_sample_rid;
%massa_E_st{1}=massa_E_st{1}(1:n_sample_rid,:);
%E_cont=E_cont(1:n_sample,:);
%pause
%%% la densità di E_st non serve
%clear E_st
% massa totale di E_st nel tempo
% una riga per ogni sample
% una colonna per ogni istante temporale
%E_st =  massa_E_st{1};
E_st=massa_E_st{1};

%clear massa_E_st
clear massa_E_st

%time_downsample=time;


if n_sd>1
    error('questo codice funziona con un solo special disc')
end

% area complessiva delle incisure
a_inc=sum(inc(1,:)^2*inc(3,:)/2);
% diffusività efficaci (utilizzate solo se flag_gws è false)
a_tot=(1-theta_0)*pi*R^2+2*pi*R*sigma*epsilon_0+a_inc;
a_patent=2*pi*R*sigma*epsilon_0+a_inc;
kk_u_eff=kk_u* a_patent/a_tot;
kk_v_eff=kk_v* a_patent/a_tot;
% nota: per rendere questi risultati uguali a quelli del codice generale,
% occorrerebbe inoltre aggiungere la massa del disco speciale, che va a
% sommarsi a quella del primo nodo, sicché la matrice delle masse del
% sistema alle differenze finite non risulterebbe più scalare
% qui invece la massa del disco speciale è trascurata, e la matrice delle
% masse risulta scalare (quindi matrice identica)

% organizzazione del vettore delle incognite:
% valori di u da 1 a nz+1
% valori di v da nz+1+1 a nz+1+nz+1
% se flag_nonlin è true, u=[cGMP] e v=[Ca]
% se flag_nonlin è false, u=[cGMP]-[cGMP_ss] e v=[Ca]-[Ca_ss]

% dato iniziale
if flag_nonlin
    solin=[u_ss*ones(nz+1,1);v_ss*ones(nz+1,1)];
else
    solin=[zeros(nz+1,1);zeros(nz+1,1)];
end

if flag_theta
    % integra con il metodo theta
    % tutti i campioni contemporaneamente
    fprintf('\nIntegra nel tempo - Metodo Theta\n')
    
    % istanti di calcolo
    t_step=t_fin/n_step_t;
    time=(0:n_step_t)*t_step;

    if ~flag_gws
        % modello TWS: è presente la diffusione
        % matrice K che, moltiplicata per un vettore lungo nz+1, restituisce
        % l'opposto del laplaciano di u
        e = ones(nz+1,1);
        A = spdiags([-e +2*e -e], -1:1, nz+1, nz+1);
        % impone le condizioni di flusso nullo ai bordi
        A(1,2)=-2;
        A(nz+1,nz)=-2;
        % matrice di rigidezza globale della struttura
        % monta kk_u_eff/dz^2*A nel blocco (u,u) 
        % e kk_v_eff/dz^2*A nel blocco (v,v) 
        [I,J,V]=find(A);
        K=sparse([I;nz+1+I], [J;nz+1+J], [kk_u_eff/dz^2*V; kk_v_eff/dz^2*V]);

        % fattorizzazione LU della matrice risolvente del metodo theta
        % la matrice delle masse è l'identità
        [L,U]=lu(speye(2*(nz+1))+(theta*t_step)*K);
    end

    % replica solin per tutti i campioni
    % una riga per ogni punto di discretizzazione
    % una colonna per ogni campione
    solin=repmat(solin,1,n_sample);

    % corrente totale a tempo zero
    curr_tot(:,1)=current(solin, ...
        flag_gws, flag_nonlin, flag_curr_nonlin, ...
        R, nz, dz, Sigma_rod, ...
        u_ss, v_ss, ...
        j_cg_max, m_cg, K_cg, ...
        j_ex_sat, K_ex);
    
    % ciclo principale
    passo_downsample=1;
    rem_downsample=0;
    for passo=1:n_step_t
        
        % istante finale dello step temporale
        t=time(passo+1);
        
		% mostra lo stato di avanzamento
		fprintf('Diffusione nel cytosol: istante %6.4f di %6.4f\n',t,t_fin);

        % valore theta di E_st
        % è una riga, una componente per ogni campione
        f=(theta+rem_downsample)/downsample;
        E_st_th=(E_st(:,passo_downsample)*(1-f)+E_st(:,passo_downsample+1)*f)';
        
		% calcolo di u e v alla fine dello step di integrazione, in applicazione del theta-metodo
		
		% soliter sarà l'incognita a fine step
		% come valore di tentativo, prendo l'incognita ad inizio step
		soliter=solin;
		
		% iterazione sul punto fisso (soliter)
		iter=0;
        again=repmat(true,1,n_sample);
		while any(again),
			
            % stima dell'incognita all'istante t_theta
            solth=solin(:,again)*(1-theta)+soliter(:,again)*theta;
            
            % separa il vettore delle u dal vettore delle v
            u_th=solth(      1:nz+1 ,:);
            v_th=solth(nz+1+(1:nz+1),:);
		
            % assembla i termini noti
            [S_u,S_v]=sorgenti(u_th, v_th, E_st_th(again), ...
                flag_gws, flag_nonlin, saturation_effect, dz, ...
                epsilon_0, nu, V_cyto, ...
                k_hyd, PDE_s, k_st, ...
                alpha_max, alpha_min, m_cyc, k_cyc, ...
                B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
                j_ex_sat, K_ex, ...
                a, u_ss, H);
            % alloca in un unico vettore
            f_th=[S_u;S_v];

            % theta-metodo
            if flag_gws
                solfin(:,again) = solin(:,again) + t_step*f_th;
            else
                b=solin(:,again) + t_step*(f_th - (1-theta)*(K*solin(:,again)));
                solfin(:,again) = U \ (L \ b);
            end
            
            % rilassamento
            solfin(:,again)=(1-alpha)*soliter(:,again)+alpha*solfin(:,again);

            % clamp del calcio
            if flag_Ca_clamp
                % riporta la concentrazione del calcio in tutti i nodi 
                % al valore dello steady-state
                if flag_nonlin
                    solfin(nz+1+1:end,again)=v_ss*ones(nz+1,length(find(again)));
                else
                    solfin(nz+1+1:end,again)=zeros(nz+1,length(find(again)));
                end
            end
            
            % fine se raggiunta la convergenza: norma L^inf
            if norma_inf
                % norma L^inf
                again= (max(abs(solfin-soliter),[],1)./max(abs(solfin),[],1))>tol_fix;
            else
                % norma L^1
                again= (sum(abs(solfin-soliter),1)./sum(abs(solfin),1))>tol_fix;
            end
            
            % per il passo successivo
            soliter(:,again)=solfin(:,again);
            iter=iter+1;
    		fprintf('        iterazione %2i, campioni residui %4i\n',iter,length(find(again)));
            
		end % while punto fisso
        
        % per il passo temporale successivo
        solin=solfin;
    
        % salva i risultati sottocampionando
        rem_downsample=rem_downsample+1;
        if rem_downsample==downsample
            % aggiorna passo_downsample
            passo_downsample=passo_downsample+1;
            rem_downsample=0;

            % corrente totale
            curr_tot(:,passo_downsample)=current(solfin, ...
                flag_gws, flag_nonlin, flag_curr_nonlin, ...
                R, nz, dz, Sigma_rod, ...
                u_ss, v_ss, ...
                j_cg_max, m_cg, K_cg, ...
                j_ex_sat, K_ex);
        end

    end;

else
    % metodo ode
    fprintf('\nIntegra nel tempo - Metodo ode\n')
    % ciclo su tutti i campioni
    curr_tot=zeros(n_sample,length(time_downsample));
    for samp=1:n_sample
        fprintf('Campione aleatorio %4i di %4i\n',samp,n_sample);

        % integra nel tempo
        [time_downsample,sol] = ode45(@odefun, time_downsample, solin, [],...
            flag_gws, nz, dz, kk_u_eff, kk_v_eff, ...
            flag_nonlin, saturation_effect, ...
            epsilon_0, nu, V_cyto, ...
            k_hyd, PDE_s, k_st, ...
            alpha_max, alpha_min, m_cyc, k_cyc, ...
            B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
            j_ex_sat, K_ex, ...
            a, u_ss, H, ...
            E_st(samp,:), time_downsample, t_fin);
        % mette i tempi sulle colonne
        sol=sol';
        % calcola la corrente
        curr_tot(samp,:)=current(sol, ...
            flag_gws, flag_nonlin, flag_curr_nonlin, ...
            R, nz, dz, Sigma_rod, ...
            u_ss, v_ss, ...
            j_cg_max, m_cg, K_cg, ...
            j_ex_sat, K_ex);
    end
end % if flag_theta

% calcola il drop
curr_dark=curr_tot(:,1);
drop=100*(1-curr_tot./repmat(curr_dark,1,length(time_downsample)));
% informativa
[max_drop,I]=max(drop,[],2);
time_to_peak=(I-1)*t_fin/n_step_t*downsample;


save(nome, ...
    'flag_gws', 'flag_nonlin', 'saturation_effect', 'flag_curr_nonlin', 'flag_Ca_clamp', ...
    'drop', 'max_drop', 'time_to_peak', 'time_downsample')

figure(1)
hold on
plot(time_downsample, drop, 'g', 'Linewidth', 1)
set(gca,'Xlim',[0,t_fin])

return




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% calcola il secondo membro del sistema
% usata da ode
function f=odefun(t,x,...
    flag_gws, nz, dz, kk_u_eff, kk_v_eff, ...
    flag_nonlin, saturation_effect, ...
    epsilon_0, nu, V_cyto, ...
    k_hyd, PDE_s, k_st, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex, ...
    a, u_ss, H, ...
    E_st_sample, time, t_fin)

% mostra lo stato di avanzamento
fprintf('Diffusione nel cytosol: istante %6.4f di %6.4f\n',t,t_fin);

% calcola E_st_tot_sample all'istante t richiesto, con interpolazione lineare
E_st=interp1(time,E_st_sample,t,'linear');

% gradi di libeertà corrispondenti alle incognite u e v
gdl_u=      1:nz+1;
gdl_v=nz+1+(1:nz+1);

% estrae il vettore delle u e delle v dal vettore complessivo x
u=x(gdl_u); % cG
v=x(gdl_v); % Ca

% calcola il secondo membro del sistema
[S_u,S_v]=sorgenti(u, v, E_st, ...
    flag_gws, flag_nonlin, saturation_effect, dz, ...
    epsilon_0, nu, V_cyto, ...
    k_hyd, PDE_s, k_st, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex, ...
    a, u_ss, H);

% termine dovuto alla diffusione
S_u=S_u+ kk_u_eff/dz^2*[2*(u(2,:)-u(1,:)); diff(u,2); 2*(u(nz,:)-u(nz+1,:))];
S_v=S_v+ kk_v_eff/dz^2*[2*(v(2,:)-v(1,:)); diff(v,2); 2*(v(nz,:)-v(nz+1,:))];

% vettore da uguagliare a du/dt; dv/dt
f=[S_u;S_v];
return


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% calcola il secondo membro del sistema
function [S_u,S_v]=sorgenti(u, v, E_st, ...
    flag_gws, flag_nonlin, saturation_effect, dz, ...
    epsilon_0, nu, V_cyto, ...
    k_hyd, PDE_s, k_st, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex, ...
    a, u_ss, H)

% sorgenti non dipendenti dall'attivazione
if flag_nonlin
    % valore di dcG/dt
    S_u= -k_hyd*(PDE_s/(1/2*nu*epsilon_0))*u ...
          +(alpha_min+(alpha_max-alpha_min)*k_cyc^m_cyc./(k_cyc^m_cyc+v.^m_cyc));
    % valore di dCa/dt
    S_v= -j_ex_sat/(B_ca*F*V_cyto)*v./(K_ex+v) ...
        + j_cg_max*f_ca/(2*B_ca*F*V_cyto)*u.^m_cg./(K_cg^m_cg+u.^m_cg);
    % nota: per rendere questi risultati uguali a quelli del codice, 
    % la sorgente S_u va moltiplicata per: ((1-theta_0)*pi*R^2*H+2*pi*R^2*1/2*nu*epsilon_0)/V,
    % dove V=(1-theta_0)*pi*R^2*H+nu*epsilon_0*pi*R^2+sigma*epsilon_0*2*pi*R*H+sum(inc(3,:)*inc(1,:)^2/2)*H
    % mentre la sorgente S_v va moltiplicata per V_cyto/V
else
    % valore di dcG/dt
    S_u=a(1,1)*u+a(1,2)*v;
    % valore di dCa/dt
    S_v=a(2,1)*u+a(2,2)*v;
end

% termine dipendente dall'attivazione
% valore di u che entra in tale termine
if flag_nonlin
    u_actv=u_ss;
    if saturation_effect
        u_actv=u;
    end
else
    u_actv=u_ss;
    if saturation_effect
        u_actv=u_ss+u;
    end
end
% nel modello TWS solo la u nel primo nodo determina l'attivazione
% veceversa, nel modello GWS questo comando è ininfluente
u_actv=u_actv(1,:);
% calcola il termine dipendente dall'attivazione
actv=-k_st*E_st.*u_actv/V_cyto;
% nota: per rendere questi risultati uguali a quelli del codice, 
% il termine actv va moltiplicato per V_cyto/V
if ~flag_gws
    % well-stirred trasversale
    % la diesterasi è data come flusso in x=0
    % la condizione al bordo discretizzata nel primo nodo vale:
    % kk_u_eff*(u_2-u_0)/(2*dz)=(1/2)*k_st*E_st*u_actv/V_cyto*H
    % il fattore 1/2 è presente perché studiamo solo 1/2 rod, per simmetria
    % il fattore H deriva dal fatto che V_cyto=(1-theta_0)*pi*R^2*H
    % dalla condizione sopra si ricava u_0 e si sostituisce nel termine
    % dovuto alla diffusione,relativo al primo nodo:
    % kk_u_eff*(u_0-2*u_1-u_2)/dz^2
    % ne risulta il seguente termine non omogeneo:
    actv=(1/2)*actv*H * (2/dz);
end
% somma il termine dovuto all'attivazione
S_u(1,:)=S_u(1,:)+actv;
% la S_v non ha altri contributi
return



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function curr_tot=current(sol, ...
        flag_gws, flag_nonlin, flag_curr_nonlin, ...
        R, nz, dz, Sigma_rod, ...
        u_ss, v_ss, ...
        j_cg_max, m_cg, K_cg, ...
        j_ex_sat, K_ex)

% sol(punto,tempo), oppure (punto, campione)
% la storia della u è nelle righe      (1:nz+1)
% la storia della v è nelle righe nz+1+(1:nz+1)
u=sol(      1:nz+1 ,:);
v=sol(nz+1+(1:nz+1),:);
if ~flag_nonlin
    % somma lo steady-state
    u=u+u_ss;
    v=v+v_ss;
end
% storia delle correnti
if flag_curr_nonlin
    % la corrente è calcolate in modo esatto
    dens_curr_cGMP = j_cg_max/Sigma_rod*u.^m_cg./(K_cg^m_cg+u.^m_cg);
    dens_curr_Ca   = j_ex_sat/Sigma_rod*v./(K_ex+v);
else
    % la corrente è calcolate con una formula linearizzata
    dens_curr_cGMP = j_cg_max/Sigma_rod*( u_ss^m_cg/(K_cg^m_cg+u_ss^m_cg) + (K_cg^m_cg*m_cg*u_ss.^(m_cg-1)./(K_cg^m_cg+u_ss.^m_cg).^2).*(u-u_ss));
    dens_curr_Ca   = j_ex_sat/Sigma_rod*( v_ss/(K_ex+v_ss)                + (K_ex./(K_ex+v_ss).^2)                                    .*(v-v_ss));
end
% corrente totale, per tutti i campioni
if flag_gws
    curr_tot=(dens_curr_cGMP+dens_curr_Ca)*Sigma_rod;
else
    % moltiplica per 2 l'integrale in dz da 0 ad H/2
    curr_tot=2*trapz(dens_curr_cGMP+dens_curr_Ca)*dz*2*pi*R;
end
return




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function plot_fig
% dati E_st

% E_st_int=true;
E_st_int=false;
if E_st_int
    nome_E_st='\int_0^t E*(\tau)d\tau';
else
    nome_E_st='E*(t)';
end

% carica la distribuzione Monte Carlo di massa di E_st
load E_st_mouse
n_sample_rid=1000;
disp(['riduco i campioni a ',num2str(n_sample_rid)])
n_sample=n_sample_rid;
massa_E_st{1}=massa_E_st{1}(1:n_sample_rid,:);
pause
% la densità di E_st non serve
clear E_st
% massa totale di E_st nel tempo
% una riga per ogni sample
% una colonna per ogni istante temporale
E_st =  massa_E_st{1};

% dati ftttf
load ftttf
drop_ftttf        =drop;
max_drop_ftttf    =max_drop;
time_to_peak_ftttf=time_to_peak;
% scala dei tempi
[max_drop,I]=max(mean(drop_ftttf));
scala_ftttf=(I-1)*time_downsample(end)/length(time_downsample);

% dati ttttf
load ttttf
drop_ttttf        =drop;
max_drop_ttttf    =max_drop;
time_to_peak_ttttf=time_to_peak;
% scala dei tempi
[max_drop,I]=max(mean(drop_ttttf));
scala_ttttf=(I-1)*time_downsample(end)/length(time_downsample);

% dati ffttf
load ffttf
drop_ffttf        =drop;
max_drop_ffttf    =max_drop;
time_to_peak_ffttf=time_to_peak;
% scala dei tempi
[max_drop,I]=max(mean(drop_ffttf));
scala_ffttf=(I-1)*time_downsample(end)/length(time_downsample);

% % % % dati ftftf
% % % load ftftf
% % % drop_ftftf        =drop;
% % % max_drop_ftftf    =max_drop;
% % % time_to_peak_ftftf=time_to_peak;
% % % % scala dei tempi
% % % [max_drop,I]=max(mean(drop_ftftf));
% % % scala_ftftf=(I-1)*time_downsample(end)/length(time_downsample);

% % % % dati fttff
% % % load fttff
% % % drop_fttff        =drop;
% % % max_drop_fttff    =max_drop;
% % % time_to_peak_fttff=time_to_peak;
% % % % scala dei tempi
% % % [max_drop,I]=max(mean(drop_fttff));
% % % scala_fttff=(I-1)*time_downsample(end)/length(time_downsample);

% dati ftttt
load ftttt
drop_ftttt        =drop;
max_drop_ftttt    =max_drop;
time_to_peak_ftttt=time_to_peak;
% scala dei tempi
[max_drop,I]=max(mean(drop_ftttt));
scala_ftttt=(I-1)*time_downsample(end)/length(time_downsample);

% media
% tempi naturali
figure(1)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (mu of drop) vs. time')
a=plot(time_downsample,mean(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample,mean(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample,mean(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample,mean(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample,mean(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample,mean(drop_ftttt),'k','Linewidth',2);
% % % legend([a,b,c,d,e,f],'ftttf', 'ttttf','ffttf','ftftf','fttff','ftttt')
legend([a,b,c,f],'TWS', 'GWS','TWS, lin','TWS, clamp Ca','Location','NorthEast')
xlabel('time [s]')
ylabel('mu of drop')
axis([0 1 0 13])
print -depsc mu_tnat

% media
% tempi in unità time_to_peak
figure(2)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (mu of drop) vs. (time / time to peak)')
a=plot(time_downsample/scala_ftttf,mean(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample/scala_ttttf,mean(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample/scala_ffttf,mean(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample/scala_ftftf,mean(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample/scala_fttff,mean(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample/scala_ftttt,mean(drop_ftttt),'k','Linewidth',2);
% % % legend([a,b,c,d,e,f],'ftttf', 'ttttf','ffttf','ftftf','fttff','ftttt')
legend([a,b,c,f],'TWS', 'GWS','TWS, lin','TWS, clamp Ca','Location','NorthEast')
xlabel('time / time to peak')
ylabel('mu of drop')
axis([0 3 0 13])
print -depsc mu_tpeak


% std
% tempi naturali
figure(3)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (std of drop) vs. time')
a=plot(time_downsample,std(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample,std(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample,std(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample,std(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample,std(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample,std(drop_ftttt),'k','Linewidth',2);
legend([a,b,c,f],'TWS', 'GWS','TWS, lin','TWS, clamp Ca','Location','NorthEast')
xlabel('time [s]')
ylabel('std of drop')
axis([0 1 0 13])
print -depsc std_tnat

% std
% tempi in unità time_to_peak
figure(4)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (std of drop) vs. (time / time to peak)')
a=plot(time_downsample/scala_ftttf,std(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample/scala_ttttf,std(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample/scala_ffttf,std(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample/scala_ftftf,std(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample/scala_fttff,std(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample/scala_ftttt,std(drop_ftttt),'k','Linewidth',2);
legend([a,b,c,f],'TWS', 'GWS','TWS, lin','TWS, clamp Ca','Location','NorthEast')
xlabel('time / time to peak')
ylabel('std of drop')
axis([0 3 0 13])
print -depsc std_tpeak


% CV del drop(t)
% tempi naturali
figure(5)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (CV of drop) vs. time')
a=plot(time_downsample,std(drop_ftttf)./mean(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample,std(drop_ttttf)./mean(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample,std(drop_ffttf)./mean(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample,std(drop_ftftf)./mean(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample,std(drop_fttff)./mean(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample,std(drop_ftttt)./mean(drop_ftttt),'k','Linewidth',2);
% qui plotta anche il CV della E*
if E_st_int
    g=plot(time_downsample,std(cumtrapz(E_st,2))./mean(cumtrapz(E_st,2)),'r-*','Linewidth',2);
else
    g=plot(time_downsample,std(E_st,[],1)./mean(E_st,1),'y-*','Linewidth',2);
end
legend([a,b,c,f,g],'TWS', 'GWS','TWS, lin','TWS, clamp Ca',nome_E_st,'Location','SouthEast')
xlabel('time [s]')
ylabel('CV of drop')
axis([0 1 0 1])
print -depsc cv_tnat

% CV del drop(t)
% tempi in unità time_to_peak
figure(6)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (CV of drop) vs. (time / time to peak)')
a=plot(time_downsample/scala_ftttf,std(drop_ftttf)./mean(drop_ftttf),'r','Linewidth',2);
b=plot(time_downsample/scala_ttttf,std(drop_ttttf)./mean(drop_ttttf),'g','Linewidth',2);
c=plot(time_downsample/scala_ffttf,std(drop_ffttf)./mean(drop_ffttf),'b','Linewidth',2);
% % % d=plot(time_downsample/scala_ftftf,std(drop_ftftf)./mean(drop_ftftf),'c','Linewidth',2);
% % % e=plot(time_downsample/scala_fttff,std(drop_fttff)./mean(drop_fttff),'m','Linewidth',2);
f=plot(time_downsample/scala_ftttt,std(drop_ftttt)./mean(drop_ftttt),'k','Linewidth',2);
% qui plotta anche il CV della E*
if E_st_int
    g=plot(time_downsample/scala_ftttf,std(cumtrapz(E_st,2))./mean(cumtrapz(E_st,2)),'r-*','Linewidth',2);
else
    g=plot(time_downsample/scala_ftttf,std(E_st,[],1)./mean(E_st,1),'y-*','Linewidth',2);
end
legend([a,b,c,f,g],'TWS', 'GWS','TWS, lin','TWS, clamp Ca',nome_E_st,'Location','SouthEast')
xlabel('time / time to peak')
ylabel('CV of drop')
axis([0 3 0 1])
print -depsc cv_tpeak



% CV del \int_0^\t drop(\tau) d\tau
% tempi naturali
figure(7)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (CV of \int_0^t drop(\tau)d\tau) vs. time')
a=plot(time_downsample,std(cumtrapz(drop_ftttf,2))./mean(cumtrapz(drop_ftttf,2)),'r','Linewidth',2);
b=plot(time_downsample,std(cumtrapz(drop_ttttf,2))./mean(cumtrapz(drop_ttttf,2)),'g','Linewidth',2);
c=plot(time_downsample,std(cumtrapz(drop_ffttf,2))./mean(cumtrapz(drop_ffttf,2)),'b','Linewidth',2);
% % % d=plot(time_downsample,std(cumtrapz(drop_ftftf,2))./mean(cumtrapz(drop_ftftf,2)),'c','Linewidth',2);
% % % e=plot(time_downsample,std(cumtrapz(drop_fttff,2))./mean(cumtrapz(drop_fttff,2)),'m','Linewidth',2);
f=plot(time_downsample,std(cumtrapz(drop_ftttt,2))./mean(cumtrapz(drop_ftttt,2)),'k','Linewidth',2);
% qui plotta anche il CV della E*
if E_st_int
    g=plot(time_downsample,std(cumtrapz(E_st,2))./mean(cumtrapz(E_st,2)),'r-*','Linewidth',2);
else
    g=plot(time_downsample,std(E_st,[],1)./mean(E_st,1),'y-*','Linewidth',2);
end
legend([a,b,c,f,g],'TWS', 'GWS','TWS, lin','TWS, clamp Ca',nome_E_st,'Location','SouthEast')
xlabel('time [s]')
ylabel('CV of \int_0^t drop(\tau)d\tau')
axis([0 1 0 1])
print -depsc int_cv_tnat

% CV del \int_0^\t drop(\tau) d\tau
% tempi in unità time_to_peak
figure(8)
hold on
box on
set(gca,'Fontsize',14)
title('Mouse, (CV of \int_0^t drop(\tau)d\tau) vs. (time / time to peak)')
a=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_ftttf,2))./mean(cumtrapz(drop_ftttf,2)),'r','Linewidth',2);
b=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_ttttf,2))./mean(cumtrapz(drop_ttttf,2)),'g','Linewidth',2);
c=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_ffttf,2))./mean(cumtrapz(drop_ffttf,2)),'b','Linewidth',2);
% % % d=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_ftftf,2))./mean(cumtrapz(drop_ftftf,2)),'c','Linewidth',2);
% % % e=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_fttff,2))./mean(cumtrapz(drop_fttff,2)),'m','Linewidth',2);
f=plot(time_downsample/scala_ftttf,std(cumtrapz(drop_ftttt,2))./mean(cumtrapz(drop_ftttt,2)),'k','Linewidth',2);
% qui plotta anche il CV della E*
if E_st_int
    g=plot(time_downsample/scala_ftttf,std(cumtrapz(E_st,2))./mean(cumtrapz(E_st,2)),'r-*','Linewidth',2);
else
    g=plot(time_downsample/scala_ftttf,std(E_st,[],1)./mean(E_st,1),'y-*','Linewidth',2);
end
legend([a,b,c,f,g],'TWS', 'GWS','TWS, lin','TWS, clamp Ca',nome_E_st,'Location','SouthEast')
xlabel('time / time to peak')
ylabel('CV of \int_0^t drop(\tau)d\tau')
axis([0 3 0 1])
print -depsc int_cv_tpeak

% tabella
disp(['CV del max di E* e dell''integrale da 0 a ',num2str(time_downsample(end)),' s di E*'])
disp([std(max(E_st,[],2))./mean(max(E_st,[],2)), std(trapz(E_st,2))./mean(trapz(E_st,2))])

disp(['CV del max drop, e di \int_0^',num2str(time_downsample(end)),' s drop(tau) dtau, caso TWS'])
disp([std(max_drop_ftttf)./mean(max_drop_ftttf), std(trapz(drop_ftttf))./mean(trapz(drop_ftttf)) ]);

disp(['CV del max drop, e di \int_0^',num2str(time_downsample(end)),' s drop(tau) dtau, caso GWS'])
disp([std(max_drop_ttttf)./mean(max_drop_ttttf), std(trapz(drop_ttttf))./mean(trapz(drop_ttttf)) ]);

disp(['CV del max drop, e di \int_0^',num2str(time_downsample(end)),' s drop(tau) dtau, caso TWS, lin'])
disp([std(max_drop_ffttf)./mean(max_drop_ffttf), std(trapz(drop_ffttf))./mean(trapz(drop_ffttf)) ]);

disp(['CV del max drop, e di \int_0^',num2str(time_downsample(end)),' s drop(tau) dtau, caso TWS, Ca clamp'])
disp([std(max_drop_ftttt)./mean(max_drop_ftttt), std(trapz(drop_ftttt))./mean(trapz(drop_ftttt)) ]);


