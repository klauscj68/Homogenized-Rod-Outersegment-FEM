function check_exp_R_st
% verifica delle funzioni che danno la distribuzione statistica delle R* in ciascuno stato

close all

% aggiunge alla path la directory id
path('../id',path);

% parametri su cui è fatto il check

% numero di stati attivi della R_st (lo stato spento è n_step_R+1)
n_step_R=4;

% durate medie dei singoli stati
% step tutti uguali o tutti diversi
% equal_step=true;
equal_step=false;
% durate medie dei singoli passi
if equal_step
    % step tutti uguali
    k_R=2.8;
    tau=(1/k_R)/n_step_R;
else
    % step tutti diversi
    tau=[0.2, 0.1, 0.05, 0.25];
    if abs(det(vander(tau)))<1000*eps
        error('Error in the data: in the case of different steps, no two equal step are allowed')
    end
end

% istanti di caolcolo
n_step_t=10;
t_fin=1;
t_step=t_fin/n_step_t;

% numero di campioni (metodo Monte Carlo)
N_sample=100000;

% istanti di tempo su cui è fatto il check
time=(0:n_step_t)'*t_step;



% PRIMO METODO
% distribuzioni delle R* nei vari stati, come ottenute, tempo per tempo, dal limite statistico
% usato da diffusione_ipd_PDE, mode_time=3
[fraction_exppdf]=prob_state(n_step_R,tau,equal_step,time);



% SECONDO METODO
% calcolo delle stesse distribuzioni come media statistica degli stati ottenuti in modo random
% durate effettive dei singoli stati: una colonna per ogni esperimento
% usato da diffusione_ipd_PDE, mode_time=2
t=zeros(n_step_R,N_sample);
if equal_step
    t=exprnd(tau,n_step_R,N_sample);
else
    for i=1:n_step_R
        t(i,:)=exprnd(tau(i),1,N_sample);
    end
end
% tempi di transizione
t_transition=cumsum(t,1);
% matrice delle percentuali in ciascuno stato
fraction_MC1=zeros(n_step_t+1,n_step_R);
% è necessario un ciclo sui tempi
for passo=1:n_step_t+1
    % calcola il vettore state
    state=sum(time(passo)>t_transition,1)+1;
    % calcola la frazione di R* in ciascuno stato
    for i=1:n_step_R
        fraction_MC1(passo,i)=length(find(state==i))/N_sample;
    end
end


% TERZO METODO
% vettori riga di stato
% state ha tante componenti quante sono le rodopsine inizialmente accese
% Indica in che stato si trova ciascuna rodopsina inizialmente accesa
% Per R* i possibili stati sono espressi da un intero compreso tra 1 e n_step_R+1
% state_R_st=1 per le rodopsine accese non ancora fosforilate
% state_R_st=n_step_R+1 per le rodopsine spente
% il vettore state è aggiornato ad ogni passo temporale
% usato da diffusione_ipd_MC, mode_time=2
state=ones(1,N_sample);
% matrice delle percentuali in ciascuno stato
fraction_MC2=zeros(n_step_t+1,n_step_R);
% probabilità che una R* subisca 0,1,2,...,n_step_R transizioni durante lo step di tempo t_step
% ogni colonna si riferisce ad un differente stato iniziale
% ogni riga si riferisce ad un differente stato finale
prob=zeros(n_step_R,n_step_R);
for i=1:n_step_R
    if equal_step
        [prob(i:n_step_R,i)]=prob_state(n_step_R-i+1,tau,true,t_step)';
    else
        [prob(i:n_step_R,i)]=prob_state(n_step_R-i+1,tau(i:end),false,t_step)';
    end
end
prob=cumsum(prob,1);

% è necessario un ciclo sui tempi
for passo=1:n_step_t+1
    % calcola la frazione di R* in ciascuno stato
    for i=1:n_step_R
        fraction_MC2(passo,i)=length(find(state==i))/N_sample;
    end
    % transizioni
    prob_sample=rand(1,N_sample);
    for s=1:N_sample
        if state(s)<=n_step_R
            % ancora viva
            state(s)=state(s)+sum(repmat(prob_sample(s),n_step_R-state(s)+1,1)-prob(state(s):n_step_R,state(s))>0,1);
        end
    end
end




% figure
for i=1:n_step_R
    figure(i)
    hold on
    title(['State ',num2str(i)])
    plot(time,fraction_exppdf(:,i),'r','Linewidth',3)
    plot(time,fraction_MC1(:,i),'g','Linewidth',2)
    plot(time,fraction_MC2(:,i),'b','Linewidth',1)
    axis([0 max(time) 0 1])
end
figure(n_step_R+1)
hold on
title('Dead')
plot(time,1-sum(fraction_exppdf,2),'r','Linewidth',3)
plot(time,1-sum(fraction_MC1,2),'g','Linewidth',2)
plot(time,1-sum(fraction_MC2,2),'b','Linewidth',1)
axis([0 max(time) 0 1])

return
