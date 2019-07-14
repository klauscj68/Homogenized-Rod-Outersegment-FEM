% Assemblaggio del vettore dei carichi esterni per il problema omogeneizzato
function [f_th]=carichi(R, H, n_sez, n_sd, ...
    n_pd, n_os, ...
    n_gdl, gdl_vol, gdl_sd, gdl_os, ...
    epsilon_0, nu, ...
    k_hyd, PDE_s, k_st, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex, ...
    u_th, v_th, E_st_th, ...
    M_vol, M_sd, M_os, ...
    u_ss,v_ss)





% i carichi sono interpolati in accordo alla rappresentazione tramite le funzioni di forma

% parametri geometrici
theta_0=1/(1+nu);
Sigma_rod=2*pi*R*H;

% numero di campioni residui
n_sample=size(u_th,2);

% inizializza i vettori di accumulo per il calcolo del vettore dei carichi
% Ind_f_th contiene gli indici mentre f_th i corrispondenti valori
% alla fine si sommano tra loro i valori corrispondenti allo stesso indice con il comando sparse
n_f_th=n_pd*n_sez + n_pd*n_sd + n_os*n_sez;
Ind_f_th=zeros(n_f_th,1);
f_th=zeros(n_f_th,n_sample);
point_f_th=0;

% volume

% calcola la restrizione di u_th e v_th ai nodi interni
% gdl_vol ordinato in colonna
gdl_vol_col=reshape(gdl_vol,n_pd*n_sez,1);
% i nodi di volume corrispondono ai gdl in gdl_vol
u_th_vol=u_th(gdl_vol_col,:);
v_th_vol=v_th(gdl_vol_col,:);

% carico k_hyd*PDE_s*u e f1 (ciclasi) nel cilindro calcolato in corrispondenza della soluzione corrente
% fosfodiesterasi spenta
q_th =     (1-theta_0)*k_hyd*(PDE_s/(1/2*nu*epsilon_0))*(u_th_vol-u_ss);
% termine f1, che dipende da v nel problema per u (ciclasi)
q_th = q_th + (1-theta_0)*((alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v_ss^m_cyc)^2*m_cyc*v_ss^(m_cyc-1))*(v_th_vol-v_ss);
 
% assembla i carichi (problema per u: i termini noti sono sistemati nelle componenti da 1 a n_gdl)
% aggiorna il cumulo
Ind_f_th(point_f_th+1:point_f_th+n_pd*n_sez) = gdl_vol_col;
f_th(point_f_th+1:point_f_th+n_pd*n_sez,:) = M_vol*q_th;
% aggiorna il puntatore
point_f_th = point_f_th+n_pd*n_sez;

% dischi speciali

for d=1:n_sd
    % calcola la restrizione di u_th e v_th al disco speciale d
    u_th_sd=u_th(gdl_sd(:,d),:);
    v_th_sd=v_th(gdl_sd(:,d),:);

    % carico k_hyd*PDE_s*u e f1 (ciclasi) nel cilindro calcolato in corrispondenza della soluzione corrente
	% fosfodiesterasi spenta sulle due facce
	q_th =     2*k_hyd*PDE_s*(u_th_sd-u_ss);
	% fosfodiesterasi accesa
    %%%%%%%%%%%%%%%%%%q_th = q_th + k_st*E_st_th(:,:,d).*u_th_sd;
    q_th = q_th +  k_st*E_st_th(:,:,d)*u_ss;
	% termine f1, che dipende da v nel problema per u (ciclasi)
	q_th = q_th + 2*((alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v_ss^m_cyc)^2*m_cyc*v_ss^(m_cyc-1))*(1/2*nu*epsilon_0)*(v_th_sd-v_ss);

	% assembla i carichi (problema per u: i termini noti sono sistemati nelle componenti da 1 a n_gdl)
	% aggiorna il cumulo
    Ind_f_th(point_f_th+1:point_f_th+n_pd) = gdl_sd(:,d);
    f_th(point_f_th+1:point_f_th+n_pd,:) = M_sd*q_th;
    % aggiorna il puntatore
    point_f_th = point_f_th+n_pd;
end

% calcola la restrizione di u_th e v_th allo special disc
% gdl_os ordinato in colonna
gdl_os_col=reshape(gdl_os,n_os*n_sez,1);
% estrae la soluzione sui nodi di os
u_th_os=u_th(gdl_os_col,:);
v_th_os=v_th(gdl_os_col,:);

% termine g1, che dipende da v nel problema per v (scambiatore), calcolato in corrispondenza della soluzione corrente
q_th =     j_ex_sat/(Sigma_rod*B_ca*F)*K_ex/(K_ex+v_ss)^2*(v_th_os-v_ss);
% termine g2, che dipende da u nel problema per v (canali dipendenti dal cGMP), calcolato in corrispondenza della soluzione corrente
q_th = q_th - j_cg_max*f_ca/(2*Sigma_rod*B_ca*F)*m_cg*K_cg^m_cg*u_ss^(m_cg-1)/(K_cg^m_cg+u_ss^m_cg)^2*(u_th_os-u_ss);

% assembla i carichi (problema per v: i termini noti sono sistemati nelle componenti da n_gdl+1 a 2*n_gdl)
% aggiorna il cumulo
Ind_f_th(point_f_th+1:point_f_th+n_os*n_sez) = n_gdl+gdl_os_col;
f_th(point_f_th+1:point_f_th+n_os*n_sez,:) = M_os*q_th;
% aggiorna il puntatore
point_f_th = point_f_th+n_os*n_sez;

% Crea il vettore dei carichi
% f_th contiene in ogni colonna i carichi relativi ad un sample
% f_th(:) mette le colonne una sopra l'altra
% i corrsponenti indici di riga si ottengono replicando Ind_f_th, una
% replica sull'altra, per quanti sono i sameple
% gli indici di colonna saranno 1 (point_f_th volte) per il primo sample, 
% 2 (point_f_th volte) per il secondo sample, e così via 
righe=repmat(Ind_f_th,n_sample,1);
colonne=repmat(1:n_sample,point_f_th,1);
f_th=accumarray([righe, colonne(:)], f_th(:), [2*n_gdl, n_sample]);

% porta a secondo membro i termini noti
f_th=-f_th;

return
