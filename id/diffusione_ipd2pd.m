% mappa la soluzione del problema per E_st_id nei dischi speciali incisi sui dischi non incisi
function [E_st_sample]=diffusione_ipd2pd(n_pd, p_pd, n_pd_id, p_pd_id, t_pd_id, ...
    n_sd, n_inc, inc2pd_id, ipd2pd_id, ...
    n_ref_cyto, n_ref_id, n_step_t, R, E_st_id)

fprintf('\nProietta la diesterasi sulla mesh dei dischi speciali\n');

% cuce le incisure, assegnando a ciascun nodo del disco pivot, 
% la media dei valori della E* sui due labbri dell'incisura

% ciclo sulle incisure
for k=1:n_inc	
    % determina i nodi del labbro a theta-
    labbro_m=inc2pd_id{k}(2:end);
    % trova i corrispondenti nodi sul labbro a theta+
    labbro_p=n_pd_id+posizione(labbro_m,ipd2pd_id(n_pd_id+1:end));
    
	% cuce: ciclo sui dischi speciali
    for d=1:n_sd
        % assegna ai nodi di labbro_m la media dei valori corrispondenti a labbro_m e labbro_p
        E_st_id{d}(:,labbro_m)=( E_st_id{d}(:,labbro_m) + E_st_id{d}(:,labbro_p) )/2;
    end
end

% toglie dalle matrici E_st_id{d} le colonne corrispondenti ai nodi sdoppiati
for d=1:n_sd
    % conserva solo le prime n_pd_id colonne
    E_st_id{d}= E_st_id{d}(:,1:n_pd_id);
end

% proietta
if n_ref_cyto<n_ref_id
    % la mesh utilizzata per la soluzione del problema di diffusione sul disco inciso (generata con n_ref_id raffittimenti)
    % è più fine di quella usata per il problema omogeneizzato (generata con n_ref_cyto raffittimenti), sicché
	% E_st_id{d}, soluzione del problema di diffusione sul disco inciso, 
	% deve essere proiettata sulla mesh del pivot disc della diffusione volumica, fornendo E_st_sample{d}
	
	% inizializza E_st_sample(d), con le stesse righe di E_st_id{d}, ma con n_pd colonne
    E_st_sample=cell(1,n_sd);
    for d=1:n_sd
        E_st_sample(d)={zeros(size(E_st_id{d},1),n_pd)};
    end
	
    % il fine è interpolare la soluzione sulla mesh rifinita per calcolarla
    % in ogni nodo p_pd della mesh di arrivo
    % per questo, trova il triangolo della mesh rifinita cui appartiene ogni nodo
    % della mesh di arrivo, quindi combina linearmente la soluzione della
    % mesh rifinita con le funzioni di forma campionate nel nodo della mesh
    % di arrivo
    % in pratica, avendo già cucito E_st_id, usa la mesh del pivot disk non
    % inciso, ma con raffittimento n_ref_id
    % in altri termini, usa p_pd_id, t_pd_id invece di p_ipd_id, t_ipd_id
    % tri e F hanno n_pd colonne
    [tri,F]=find_tri(p_pd_id,t_pd_id,p_pd,R);
    
    % combinazione lineare della soluzione sulla mesh raffittita con le funzioni F_elem
    % lista dei nodi dei triangoli della mesh raffittita che contengono i
    % punti p_pd, in colonna
    nodi=t_pd_id(:,tri);
    nodi=nodi(:);
    % numero di campioni temporali
    % interpola linearmente i tre valori nodali di E_st_id{d}, per tutti i tempi
    for d=1:n_sd
        % soluzione calcolata nei nodi dei triangoli che contengono i punti p_pd, per tutti i tempi
        sol=E_st_id{d}(:,nodi);
        % moltiplica la soluzione sol' per le repliche nei tempi del vettore
        % F, organizzato in colonna
        % questo fornisce una matrice che ha sulle righe i nodi a tre a
        % tre, sulle colonne i tempi
        % questa matrice è riordinata in modo da avere solo 3 righe
        % la somma sulle 3 righe esegue la combinazione lineare dei 3
        % valori nodali nei 3 verici del triangolo
        % il vettore riga risultante è poi riorganizzato, a gruppi di n_pd
        % elementi per ogni istante temporale
        E_st_sample{d}=reshape(sum(reshape(repmat(F(:),1,n_step_t+1).*sol',3,n_pd*(n_step_t+1)),1),n_pd,n_step_t+1)';
    end

else
    % la mesh utilizzata per la soluzione del problema di diffusione sul disco inciso
    % è identica a quella usata per il problema omogeneizzato, sicché non è necessaria alcuna proiezione
    E_st_sample=E_st_id;
end % if

return
