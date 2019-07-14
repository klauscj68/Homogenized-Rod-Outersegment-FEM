% esempio: la distribuzione esponenziale

% numero di campioni
% rappresenta il numero di R* che sono state osservate nel loro spegnimento
n_sample=1000000;
% tasso di spegnimento di R*
k_R=1.65;
% vita media della R*
t_mean=1/k_R;
% risoluzione temporale per costruire l'istogramma
t_res=0.01;

% genera gli istanti di spegnimento, aventi media t_mean
t=exprnd(t_mean,n_sample,1);

% costruisce l'istogramma, con risoluzione temporale t_res
ind=fix(t/t_res)+1;
ind_max=max(ind);
% densità frequentistica
pdf=sparse(ind,ones(1,length(ind)),ones(1,length(ind)))/n_sample/t_res;
% istanti di tempo (centri degli intervallini)
time=((1:ind_max)-0.5)*t_res;

% funzione teorica
pdf_teorico=exppdf(time, t_mean);

% figura
close all
figure(1)
hold on
plot(time,pdf,'r','LineWidth',2)
plot(time,pdf_teorico,'g','LineWidth',2)
