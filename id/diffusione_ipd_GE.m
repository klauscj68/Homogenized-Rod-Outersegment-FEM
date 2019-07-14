% problema della creazione e diffusione di G^* ed E^*, in modo accoppiato, 
% sui dischi speciali incisi
function [sol_G_st_id, sol_E_st_id]=diffusione_ipd_GE(n_ipd_id , M, K, ...
                    cc_G_st, D_G_st, cc_E_st, D_E_st, k_E, k_GE, PDE_s, ...
                    V_G_st, sol_G_st_id_0, sol_E_st_id_0, ...
                    t_fin, n_step_t, metodo_id, solver, theta, alpha, tol_fix, inspect, tipo)
                
% Integrazione nel tempo dell'equazione di diffusione per 2 specie chimiche

% output
% sol_G_st_id e sol_E_st_id hanno su ogni riga la soluzione campionata 
% nell'istante del corrispondente elemento di time

% passo di integrazione temporale
t_step=t_fin/n_step_t;

% istanti di calcolo
time=(0:n_step_t)'*t_step;

% le subunità totali di E sono il doppio delle molecole di E
E_tot=PDE_s*2;


if metodo_id(1)==1
    % Integrazione con il metodo theta
    fprintf('\nIntegra nel tempo - Metodo Theta - Problema per %s\n',tipo)
    
	% inizializza
	sol_G_st_id=zeros(n_step_t+1, n_ipd_id);
	sol_E_st_id=zeros(n_step_t+1, n_ipd_id);
    % dato iniziale
	sol_G_st_id(1,:)=sol_G_st_id_0;
	sol_E_st_id(1,:)=sol_E_st_id_0;
	
    % solin è inizializzata al dato iniziale, in colonna prima G* poi E*
	solin=[sol_G_st_id_0';sol_E_st_id_0'];
    solfin=zeros(2*n_ipd_id,1);
    
    % ciclo principale
	tic

    % matrici delle masse e delle rigidezze del problema accoppiato
    [I_M,J_M,V_M]=find(M);
    [I_K,J_K,V_K]=find(K);
    % la matrice delle masse globale ha la struttura:
    % [cc_G_st*M  0 \\ 0  cc_E_st*M]
    M_gl=sparse([I_M;           n_ipd_id+I_M], ...
                [J_M;           n_ipd_id+J_M], ...
                [cc_G_st*V_M;   cc_E_st*V_M], 2*n_ipd_id, 2*n_ipd_id);
    % la matrice delle rigidezze globale ha la struttura:
    % [D_G*K+k_GE*E_tot*M  0 \\ -k_GE*E_tot*M  D_E*K+k_E*M]
    K_gl=sparse([I_K;        I_M;            n_ipd_id+I_K; n_ipd_id+I_M; n_ipd_id+I_M    ], ...
                [J_K;        J_M;            n_ipd_id+J_K; n_ipd_id+J_M; J_M             ], ...
                [D_G_st*V_K; k_GE*E_tot*V_M; D_E_st*V_K;   k_E*V_M;      -k_GE*E_tot*V_M ], 2*n_ipd_id, 2*n_ipd_id);
    
    % fattorizzazione
    % permuta se richiesto
    if (metodo_id(3)==1) && (metodo_id(2)~=0)
        fprintf('\nColumn approximate minimum degree permutation\n');
        p = colamd(M_gl+theta*t_step*K_gl);
    elseif (metodo_id(3)==2) && (metodo_id(2)~=0)
        error('Symmetric reverse Cuthill-McKee permutation not allowed for coupled disk diffusion');
    else
        fprintf('\nNo permutation\n');
        p=(1:n_ipd_id);
    end
    % fattorizza se richiesto
    switch metodo_id(2)
        case 1
            fprintf('\nInizio della fattorizzazione LU\n');
            [LL,UU]=lu(M_gl(:,p)+theta*t_step*K_gl(:,p));
        case 2
            error('Cholesky factorization not allowed for coupled disk diffusion');
        otherwise
            fprintf('\nNo factorization\n');
    end

    for passo=1:n_step_t
        % istante finale dello step temporale
        t=time(passo+1);
        
		% mostra lo stato di avanzamento
		fprintf('Diffusione sul disco: istante %6.4f di %6.4f\n',t,t_fin);

        % soliter sarà l'incognita a fine step
		% come valore di tentativo, prendo l'incognita ad inizio step
		soliter=solin;
		
		% iterazione sul punto fisso (soliter)
		done=0;
		iter=0;
		while ~done,
			
            % stima dell'incognita all'istante t_theta
            solth=solin*(1-theta)+soliter*theta;
            
            % separa il vettore delle G_st dal vettore delle E_st
            G_st_th=solth(1:n_ipd_id);
            E_st_th=solth(n_ipd_id+(1:n_ipd_id));
            
            % assembla il termine noto
            f_th=k_GE*(M*(G_st_th.*E_st_th));
            f_th_gl=[f_th;-f_th];
            f_th_gl(1:n_ipd_id)=f_th_gl(1:n_ipd_id)+ (1-theta)*V_G_st(passo,:)'+theta*V_G_st(passo+1,:)';
            
            % theta-metodo
            b=M_gl*solin + t_step*(f_th_gl - (1-theta)*(K_gl*solin));
            switch metodo_id(2)
                % metodo_id(2) non può essere 2
                case 0
                    solfin=(M_gl+theta*t_step*K_gl) \ b ;
                case 1
                    solfin(p)=UU\(LL \ b);
            end

            % rilassamento
            solfin=(1-alpha)*soliter+alpha*solfin;

            % fine se raggiunta la convergenza
            err=norm((solfin-soliter))/norm(solfin);
            done=(err<tol_fix);
            
            % per il passo successivo
            soliter=solfin;
            iter=iter+1;
    		fprintf('        iterazione %i\n',iter);
            
		end % while punto fisso
            
        % salva in sol
        sol_G_st_id(passo+1,:)=solfin(1:n_ipd_id);
        sol_E_st_id(passo+1,:)=solfin(n_ipd_id+(1:n_ipd_id));

        % per il passo temporale successivo
        solin=solfin;
    end;
	toc
    presskey(inspect);
else 
    % metodo ode
    fprintf('\nIntegra nel tempo - Metodo %s - Problema per %s\n',solver,tipo)

    tic

    % integra
    error('non ancora implementato')
    
	toc
    presskey(inspect);

end % if metodo

return
