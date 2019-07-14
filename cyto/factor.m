% assembla le matrici globali e predispone per l'integrazione nel tempo
function [M_vol, M_sd, M_os, M_gl, K_gl, LL, UU, RR, p]=...
    factor(R, H, n_sez, n_sd, ...
    n_pd, n_tri, n_os, n_inc, n_p_inc, ...
    os2pd, inc2pd, ...
    p_pd, t_pd, Z_s, Z_sd, inc, ...
    n_gdl, gdl_vol, gdl_sd, gdl_os, gdl_inc,...
    t_fin, n_step_t, metodo_cyto, theta, ...
    epsilon_0, nu, sigma, ...
    cc_u, kk_u, cc_v, kk_v)

% M_gl, K_gl, LL, UU, RR, p sono utilizzate per l'integrazione nel tempo
% M_vol, M_sd, M_os sono utilizzate per il calcolo dei carichi

% costruzione delle matrici di massa e rigidezza di ciascuna regione: vol, sd, os, inc
tic
[M_vol, K_vol, M_sd, K_sd, M_os, K_os, M_inc, K_inc]=...
    assembla(n_pd, n_tri, n_os, n_inc, n_p_inc, n_sez, ...
    p_pd, t_pd, Z_s, os2pd, inc2pd);
toc

% % % % controlla l'assemblaggio
% % % check_MK(R, H, n_sez, n_sd, ...
% % %     n_pd, n_tri, n_os, n_inc, n_p_inc, ...
% % %     os2pd, inc2pd, ...
% % %     p_pd, t_pd, Z_s, Z_sd, inc, ...
% % %     M_vol, K_vol, M_sd, K_sd, M_os, K_os, M_inc, K_inc);

% assembla le matrici delle rigidezze e delle masse per il problema omogeneizzato
% mettendo insieme i pezzi di vol, sd, os, inc
tic
[K_gl, M_gl]=assembla_gl(n_pd, n_os, n_inc, n_p_inc, ...
    n_sez, n_sd, ...
    n_gdl, gdl_vol, gdl_sd, gdl_os, gdl_inc, ...
    epsilon_0, nu, sigma, inc, ...
    cc_u, kk_u, cc_v, kk_v, ...
    M_vol, K_vol, M_sd, K_sd, M_os, K_os, M_inc, K_inc);
toc

% recupera un po' di memoria: le K non servono più, viceversa le M servono per i carichi
clear K_vol K_sd K_os K_inc M_inc

% passo di integrazione temporale
t_step=t_fin/n_step_t;

% fattorizzazione
if metodo_cyto(1)==1

    % permuta se richiesto
    tic
    if (metodo_cyto(3)==1) && (metodo_cyto(2)~=0)
        fprintf('\nSymmetric approximate minimum degree permutation\n');
        p = symamd(M_gl+(theta*t_step)*K_gl);
    elseif (metodo_cyto(3)==2) && (metodo_cyto(2)~=0)
        fprintf('\nSymmetric reverse Cuthill-McKee permutation\n');
        p = symrcm(M_gl+(theta*t_step)*K_gl);
    else
        fprintf('\nNo permutation\n');
        p=1:2*n_gdl;
    end
    toc

    % fattorizza se richiesto
    LL=[];
    UU=[];
    RR=[];
    tic
    switch metodo_cyto(2)
        case 1
            % fattorizzazione LU
            fprintf('\nFattorizzazione LU\n');
            [LL,UU]=lu(M_gl(p,p)+(theta*t_step)*K_gl(p,p));
        case 2
            % fattorizzazione di Cholesky
            fprintf('\nFattorizzazione di Cholesky\n');
            RR=chol(M_gl(p,p)+(theta*t_step)*K_gl(p,p));
        otherwise
            % non fattorizza
            fprintf('\nNo factorization\n');
    end
    toc

else
    % metodo ode
    % matrici vuote per evitare errore
    LL=[];
    UU=[];
    RR=[];
    p=[];
end

return
