% problema della creazione di R* ed E^* sui dischi speciali incisi
% metodo Monte Carlo semplificato
% tasso di creazione nu_RE E* da parte di R* dipendente solo dallo stato di R*
function [time, R_st_id, E_st_id]=diffusione_ipd_MC(...
    n_ipd_id , n_tri_id , p_ipd_id, t_ipd_id, ...
    n_sd, n_inc, inc, ...
    D_E_st, D_R_st, ...
    n_step_R, nu_RE, k_E, tau, equal_step, mode_time, ...
    n_Phi, Phi, ...
    peak_delta, R,  ...
    plot_pool , inspect, ...
    t_fin, n_step_t)

% output
% R_st_id(d)(passo,:,j) ed E_st_id(d)(passo,:) sono le densità di R* ed E* [molecole/(\mu m)^2]
% sul disco speciale d, all'istante temporale (passo-1)*t_step, 
% in tutti i nodi della mesh p_ipd_id del pivot disk inciso/raffittito
% l'indice j della R_st_id indica lo stato, j=1..n_step_R
% quindi la densità totale della R*, in qualsiasi stato ancora attivo, si
% ottiene sommando su j

% metodo ammissibile solo per flag_model==1 (3D)

% costruzione  delle matrici di massa e rigidezza, e delle collezioni dei vettori dei carichi degli elementi
% nota: le matrici K e M non sono moltiplicate per alcun coefficiente
% contengono gli integrali dei prodotti dei gradienti delle funzioni di forma (K) o dei prodotti delle funzioni di forma (M)
fprintf('\nAssembla le matrici K, M del problema di diffusione sul disco inciso\n');
[K, M] = assembla_ipd(n_ipd_id, n_tri_id, p_ipd_id, t_ipd_id);

% % % % area di un disco speciale discretizzato
% % % area_sd=ones(1,n_ipd_id )*M*ones(n_ipd_id ,1);

% la soluzione per la R* è organizzata nel cell array R_st_id
% la matrice R_st_id{d}(:,:,j) contiene la storia della soluzione del problema per R_st nel disco speciale d
% su ogni riga: la soluzione calcolata nell'istante temporale individuato dall'omonima riga di time
% l'indice j indica lo stato, j=1..n_step_R, in cui si trova la R*

% la soluzione per la E* è organizzata nel cell array E_st_id
% la matrice E_st_id{d} contiene la storia della soluzione del problema per E_st nel disco speciale d
% su ogni riga: la soluzione calcolata nell'istante temporale individuato dall'omonima riga di time


% chiamata dummy a find_tri, per determinare le matrici pt_min, pt_max
[tri, F, pt_min, pt_max]=find_tri(p_ipd_id, t_ipd_id, [0;0], R);

% geometria delle incisure
% lunghezza
l_inc=inc(1,:);
% posizione angolare
theta_inc=inc(2,:);
cos_theta_inc=cos(theta_inc);
sin_theta_inc=sin(theta_inc);
% normale in senso antiorario
n_theta_inc=[-sin_theta_inc; cos_theta_inc];


% random walk di R*, che man mano crea E*, con un tasso nu_RE dipendente
% unicamente dallo stato di R*

% inizializza i cell array R_st_id e E_st_id
R_st_id=cell(1,n_sd);
E_st_id=cell(1,n_sd);

% ciclo su tutti i dischi speciali
for d=1:n_sd
    fprintf('\nRandom walk per R* ed E* sul disco %i\n',d);

    % predispone il biliardo (circnferenza e tracce delle incisure 
    % per il plottaggio dei random walk di R* ed E*)
    if plot_pool
        % figura del random path R*
        % serve a random_update
        f13=figure(13);
        newplot(f13)
        title(['R* random walk, disc ',num2str(d)])
        xlabel('x [\mu m]')
        ylabel('y [\mu m]')
        zlabel('R* [mol/\mu m^2]')
        axis([-R*1.1 R*1.1 -R*1.1 R*1.1])
        set(f13,'Position',[-30 154 672 504]);
        view(2)
        hold on
        axis([-R,R,-R,R])
        axis equal
        % circonferenza
        angle_c= 0:2*pi/50:2*pi;
        plot(R*cos(angle_c), R*sin(angle_c),'r-','Linewidth',2)
        % incisure
        for cont=1:n_inc
            plot([(R-l_inc(cont))*cos(theta_inc(cont)), R*cos(theta_inc(cont))], ...
                [(R-l_inc(cont))*sin(theta_inc(cont)), R*sin(theta_inc(cont))],'r-','Linewidth',2)
        end

        % figura del random path E*
        % serve a random_update
        f14=figure(14);
        newplot(f14)
        title(['R* random walk, disc ',num2str(d)])
        xlabel('x [\mu m]')
        ylabel('y [\mu m]')
        zlabel('R* [mol/\mu m^2]')
        axis([-R*1.1 R*1.1 -R*1.1 R*1.1])
        set(f14,'Position',[650 154 672 504]);
        view(2)
        hold on
        axis([-R,R,-R,R])
        axis equal
        % circonferenza
        angle_c= 0:2*pi/50:2*pi;
        plot(R*cos(angle_c), R*sin(angle_c),'r-','Linewidth',2)
        % incisure
        for cont=1:n_inc
            plot([(R-l_inc(cont))*cos(theta_inc(cont)), R*cos(theta_inc(cont))], ...
                [(R-l_inc(cont))*sin(theta_inc(cont)), R*sin(theta_inc(cont))],'r-','Linewidth',2)
        end
    end

    % inizializza le matrici R_st_id{d} e E_st_id{d}
    % primo indice: istante temporale
    % secondo indice: nodo
    % terzo indice, solo per R_st, stato
    R_st_id{d}=zeros(n_step_t+1,n_ipd_id,n_step_R);
    E_st_id{d}=zeros(n_step_t+1,n_ipd_id);

    % passo di integrazione temporale
    t_step=t_fin/n_step_t;

    % istanti di calcolo
    time=(0:n_step_t)'*t_step;

    % cammino libero medio (in 2 dimensioni spaziali) di R* ed E*
    % in uno step temporale
    dx_R_st=sqrt(2*2*D_R_st*t_step);
    dx_E_st=sqrt(2*2*D_E_st*t_step);
    % media quadratica della gaussiana monodimensionale:
    % è calcolata in modo che la gaussiana
    % bidimensionale ottenuta come prodotto cartesiano
    % di due gaussiane monodomnsionali abbia media
    % quadratica dx_R_st (resp., dx_E_st)
    sigma_R_st=dx_R_st/sqrt(2);
    sigma_E_st=dx_E_st/sqrt(2);

    % posizioni iniziali delle fotoisomerizzazioni R* e loro numero
    n_R_st=n_Phi(d);
    rho=Phi{d}(1,:);
    alpha=Phi{d}(2,:);
    p_delta_R_st=[rho.*cos(alpha); rho.*sin(alpha)];
    % colori, che si ripetono
    c_R_st=[];
    while length(c_R_st)<n_R_st
        c_R_st=[c_R_st,'bgrcmyk'];
    end
    
    % inizialmente nessuna E*
    n_E_st=0;
    p_delta_E_st=zeros(2,n_E_st);
    c_E_st=[];

    % vettori riga di stato state_R_st e state_E_st
    % state_R_st ha tante componenti quante sono le rodopsine inizialmente accese
    % state_E_st non ha componenti (nessuna E* a tempo zero)
    % Indicano in che stato si trova ciascuna rodopsina inizialmente accesa, 
    % e ciascuna E* generata
    % Per R* i possibili stati sono espressi da un intero compreso tra 1 e n_step_R+1
    % state_R_st=1 per le rodopsine accese non ancora fosforilate
    % state_R_st=n_step_R+1 per le rodopsine spente
    % Per E* i possibili stati sono espressi 1 o 2
    % state_E_st=1 per le E* accese
    % state_E_st=2 per le E* spente
    % I vettori state_R_st e state_E_st sono aggiornati ad ogni passo temporale
    state_R_st=ones(1,n_R_st);
    state_E_st=ones(1,n_E_st);


    % probabilità che una R* subisca 0,1,2,...,n_step_R transizioni 
    % durante lo step di tempo t_step
    % ogni colonna si riferisce ad un differente stato iniziale
    % ogni riga si riferisce ad un differente stato finale
    prob_state_R_st=zeros(n_step_R,n_step_R);
    for i=1:n_step_R
        if equal_step
            [prob_state_R_st(i:n_step_R,i)]=prob_state(n_step_R-i+1,tau,true,t_step)';
        else
            [prob_state_R_st(i:n_step_R,i)]=prob_state(n_step_R-i+1,tau(i:end),false,t_step)';
        end
    end
    prob_state_R_st=cumsum(prob_state_R_st,1);

    % La E* è una GTPasi e per questo il suo spegnimento avviene in un unico step
    n_step_E=1;
    prob_state_E_st=prob_state(n_step_E,1/k_E,true,t_step);  
    prob_state_E_st=cumsum(prob_state_E_st)';
    
    % contatore delle E_st_generate
    n_E_st_gen=0;
    
    % è necessario un ciclo sui tempi
    for passo=1:n_step_t+1
        
        % densità nodali di R* e E*
        delta_R_st_nod=dirac(n_ipd_id, p_ipd_id, t_ipd_id, M, p_delta_R_st, peak_delta, R, pt_min, pt_max);
        delta_E_st_nod=dirac(n_ipd_id, p_ipd_id, t_ipd_id, M, p_delta_E_st, peak_delta, R, pt_min, pt_max);

        % aggiorna le matrici delle densità nodali R_st_id{d} e E_st_id{d}
        % importa separatamente ciascuna delta_R_st
        for j=1:n_R_st
            % in ogni istante di tempo, mette il contributo
            % della j-esima delta_R_st nel vettore nodale
            % corrispondente a state_R_st(passo,j)
            % tutte le R* rimaste sono efficaci (quelle in stato n_step_R+1 vengono rimosse)
            R_st_id{d}(passo,:,state_R_st(j))=R_st_id{d}(passo,:,state_R_st(j))+delta_R_st_nod(:,j)';
        end
        % importa le delta_E_st attive insieme, essendo indistinguibili
        active_E_st= (state_E_st<=1);
        E_st_id{d}(passo,:)=E_st_id{d}(passo,:)+sum(delta_E_st_nod(:,active_E_st),2)';
        
        if plot_pool
            figure(14)
        end
        % crea nuove E* a partire dalle R* esistenti
        % i tassi di creazione sono in nu_RE
        % quindi nell'intervallo t_step crea, in media, nu_RE*t_step nuove E*
        % usa una statistica poissoniana
        % tutte le R* sono in stato <=n_step_R
        prole=poissrnd(nu_RE(state_R_st)*t_step);
        % numero di E* attualmente attive
        n_E_st=n_E_st+sum(prole);
        % numero di E* generate in totale
        n_E_st_gen=n_E_st_gen+sum(prole);
        % aggiorna p_delta_E_st e state_E_st
        for j=1:n_R_st
            % nella posizione che la mamma R* ha a fine step
            p_delta_E_st=[p_delta_E_st, repmat(p_delta_R_st(:,j),1,prole(j))];
            % accese, appena nate
            state_E_st=[state_E_st, ones(1,prole(j))];
            % col colore della mamma
            c_E_st=[c_E_st, repmat(c_R_st(j),1,prole(j))];
            if plot_pool
                % mette un cerchietto nel luogo di nasita
                plot(p_delta_R_st(1,j),p_delta_R_st(2,j),'o','Linewidth',1,'Color',c_R_st(j))
            end
        end
        
        if plot_pool
            figure(13)
        end
        % calcola il vettore random di spostamento di R*
        delta_R_st=sigma_R_st*[randn(1,n_R_st); randn(1,n_R_st)];
        % esegue lo spostamento delta_R_st
        p_delta_R_st=random_update(p_delta_R_st, delta_R_st, R, n_inc, l_inc, theta_inc, n_theta_inc, ...
            plot_pool, c_R_st, '*');

        if plot_pool
            figure(14)
        end
        % calcola il vettore random di spostamento di E*
        delta_E_st=sigma_E_st*[randn(1,n_E_st); randn(1,n_E_st)];
        % esegue lo spostamento delta_E_st
        p_delta_E_st=random_update(p_delta_E_st, delta_E_st, R, n_inc, l_inc, theta_inc, n_theta_inc, ...
            plot_pool, c_E_st, '+');
        
        % Aggiorna i vettori state_R_st e state_E_st
        % se mode_time==0, lascia lo stato a 1
        if mode_time==2
            % estrae un numero tra 0 e 1 per ogni R* ed E*
            % e lo confronta con le componenti k del vettore prob_state_R_st e
            % prob_state_E_st, che sono le probabilità cumulate 
            % che la R*,E* abbiano subito k transizioni
            % in effetti, prob_state_R_st è una matrice: prende la
            % colonna corrispondente allo stato in cui attualmente si trova la R*

            % R*
            random_test_R_st=rand(1,n_R_st);
            % transizioni
            for s=1:n_R_st
                % a differenza di check_exp_R_st, non è necessario assicurarsi 
                % che la R* sia ancora viva, perché quelle spente vengono eliminate
                % dalle itruzioni di sotto
                state_R_st(s)=state_R_st(s)+sum(repmat(random_test_R_st(s),n_step_R-state_R_st(s)+1,1)...
                    -prob_state_R_st(state_R_st(s):n_step_R,state_R_st(s))>0,1);
            end

            % E*
            random_test_E_st=rand(1,n_E_st);
            state_E_st=state_E_st+sum(repmat(random_test_E_st,n_step_E,1)-prob_state_E_st>0,1);

            % uccide le R* che sono in stato >n_step_R
            alive_R_st=find(state_R_st<=n_step_R);
            p_delta_R_st=p_delta_R_st(:,alive_R_st);
            state_R_st=state_R_st(:,alive_R_st);
            c_R_st=c_R_st(:,alive_R_st);

            % numero di R* ancora attive
            n_R_st=size(p_delta_R_st,2);
            
            % uccide le E* che sono in stato >1
            alive_E_st=find(state_E_st<=n_step_E);
            p_delta_E_st=p_delta_E_st(:,alive_E_st);
            state_E_st=state_E_st(:,alive_E_st);
            c_E_st=c_E_st(:,alive_E_st);
            
            % numero di E* ancora attive
            n_E_st=size(p_delta_E_st,2);
        end
        
        fprintf('\nDisc %2i, step %3i of %3i: R* (orig, st 1, st 2, ...): %2i,',d,passo,n_step_t+1,n_Phi(d));
        for j=1:n_step_R
            fprintf(' %2i,',length(find(state_R_st==j)));
        end
        fprintf('    E* (gen, act): %3i, %3i',n_E_st_gen, n_E_st);
        
%         p_delta_R_st
%         state_R_st
%         c_R_st
%         p_delta_E_st
%         state_E_st
%         c_E_st
        presskey(inspect);
        
    end % ciclo sui tempi
    fprintf('\n');
    
    if plot_pool
        % nuove figure per ogni disco speciale
        close(13)
        close(14)
    end
    
end % ciclo sui dischi speciali

return
