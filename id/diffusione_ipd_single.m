% problema della creazione e diffusione di una singola specie sui dischi speciali incisi
function [sol]=diffusione_ipd_single(n, M, K, V, sol_0, ...
    t_fin, n_step_t, metodo_id, solver, theta, inspect, tipo)

% Integrazione nel tempo dell'equazione di diffusione
% n è la dimensione del problema
% M è la matrice delle masse, K la matrice delle rigidezze
% la matrice V ha su ogni riga il vettore dei carichi campionato nell'istante del corrispondente elemento di time

% output
% sol ha su ogni riga la soluzione campionata nell'istante del corrispondente elemento di time

% passo di integrazione temporale
t_step=t_fin/n_step_t;

% istanti di calcolo
time=(0:n_step_t)'*t_step;


if metodo_id(1)==1
    % Integrazione con il metodo theta
    fprintf('\nIntegra nel tempo - Metodo Theta - Problema per %s\n',tipo)
    
	% inizializza
	sol=zeros(n_step_t+1, n);
    % dato iniziale
	sol(1,:)=sol_0;
	
    % solin è inizializzata al dato iniziale
	solin=sol_0';
	solfin=sol_0';
    
    % valore iniziale del termine noto
    f_in=V(1,:)';

    % ciclo principale
	tic

    % fattorizzazione
    % permuta se richiesto
    if (metodo_id(3)==1) && (metodo_id(2)~=0)
        fprintf('\nSymmetric approximate minimum degree permutation\n');
        p = symamd(M+theta*t_step*K);
    elseif (metodo_id(3)==2) && (metodo_id(2)~=0)
        fprintf('\nSymmetric reverse Cuthill-McKee permutation\n');
        p = symrcm(M+theta*t_step*K);
    else
        fprintf('\nNo permutation\n');
        p= 1:n;
    end
    % fattorizza se richiesto
    switch metodo_id(2)
        case 1
            fprintf('\nInizio della fattorizzazione LU\n');
            [LL,UU]=lu(M(p,p)+theta*t_step*K(p,p));
        case 2
            fprintf('\nInizio della fattorizzazione di Cholesky\n');
            RR=chol(M(p,p)+theta*t_step*K(p,p));
        otherwise
            fprintf('\nNo factorization\n');
    end

    for passo=1:n_step_t
        
        % istante finale dello step temporale
        t=time(passo+1);
		fprintf('Diffusione sul disco: istante %6.4f di %6.4f\n',t,t_fin);

        % valore finale del termine noto
        f_fin=V(passo+1,:)';
        
        % valore theta del termine noto
        f_th=(1-theta)*f_in+theta*f_fin;

        % theta-metodo
        b=M*solin + t_step*(f_th - (1-theta)*K*solin);
        switch metodo_id(2)
            case 0
                solfin=(M+theta*t_step*K) \ b ;
            case 1
                solfin(p)=UU\(LL \ b(p));
            case 2
                solfin(p)=RR\(RR' \ b(p));
        end

        % salva in sol
        sol(passo+1,:)=solfin';

        % per il passo temporale successivo
        solin=solfin;
        f_in=f_fin;
    end;
	toc
    presskey(inspect);
else 
    % metodo ode
    fprintf('\nIntegra nel tempo - Metodo %s - Problema per %s\n',solver,tipo)

    tic

    % matrice delle masse
    options = odeset('Mass',M);
    
    % integra
	[time, sol] = feval(solver, @odefun, time, sol_0, options, K, V, n_step_t, t_step);
    
	toc
    presskey(inspect);

end % if metodo

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% routine usata da ode nell'integrare il sistema M\dot{x}+K{x}=V:
% ritorna V-K{x}
function y=odefun(t, x, K, V, n_step_t, t_step)

i_prec= fix(t/t_step)+1;
i_succ=i_prec+1;
if i_succ>n_step_t+1
    % istante finale
    f=V(n_step_t+1,:);
else
    % istante intermedio: interpolazione lineare
    alpha=rem(t,t_step)/t_step;
    f=V(i_prec,:)*(1-alpha)+V(i_succ,:)*alpha;
end
% calcola il secondo membro del sistema
y=-K*x+f';
fprintf(['\nt=',num2str(t)])

return
