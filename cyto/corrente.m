% corrente totale, per tutti i campioni
function [curr_tot]=corrente( ...
    n_sez, n_os, n_gdl, gdl_os, ...
    solfin, ...
    R, H, ...
    j_cg_max, m_cg, K_cg, ...
    j_ex_sat, K_ex, M_os, u_ss, v_ss)
% solfin contiene la soluzione ad un fissato tempo
% ogni colonna si riferisce ad un campione stocastico

% superficie laterale
Sigma_rod=2*pi*R*H;

% storia della soluzione sull'outer shell
% gdl_os ordinato in colonna
gdl_os_col=reshape(gdl_os,n_os*n_sez,1);
% estrae la soluzione sui nodi di os
u_os=solfin(      gdl_os_col,:);
v_os=solfin(n_gdl+gdl_os_col,:);

%Sigma_rod=1;

% storia dei flussi di corrente nei nodi
dens_curr_cGMP=j_cg_max/Sigma_rod*(u_os.^m_cg)./(K_cg^m_cg+u_os.^m_cg);
dens_curr_Ca  =j_ex_sat/Sigma_rod*v_os./(K_ex+v_os);

% % Profilo di corrente lungo z
% n_p=n_os/2;
% J_z=dens_curr_cGMP(n_p:n_os:end)+dens_curr_Ca(n_p:n_os:end);
% 
% u_os_z=v_os(n_p:n_os:end);

% integra sulla superficie laterale
curr_cGMP=(ones(1,n_os*n_sez)*M_os*dens_curr_cGMP);
curr_Ca  =(ones(1,n_os*n_sez)*M_os*dens_curr_Ca);





% corrente totale, per tutti i campioni
curr_tot=curr_cGMP+curr_Ca;




return
