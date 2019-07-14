% Funzione che fornisce un vettore con n componenti. 
% La i-esima componente contiene la probabilità che, nell'intervallo (0,t), 
% la rodopsina abbia subito esattamente i-1 transizioni e quindi si trovi nello stato i
function [p]=prob_state(n, mu, lambda, t)
p=zeros(length(t),n);

N_step=length(t);

% Integrazione delle equazioni differenziali con la matrice
% esponenziale


% Matrice dei coefficienti del sistema ode
M=diag(-(lambda+mu),0)+diag(lambda(1:end-1),-1);

    
P0=[1;zeros(n-1,1)];
    
    
p(1,:)=P0';
    
    for cont=2:N_step
        t_curr=t(cont);
        val=expm(M*t_curr)*P0;
        p(cont,:)=val';
    end

return
