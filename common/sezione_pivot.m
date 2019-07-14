% mesh del disco pivot, con e senza incisure
function [n_pd, n_int_pd, n_tri, n_os, n_p_inc, ...
        os2pd, inc2pd, pd2int_pd, p_pd, t_pd, ...
        n_ipd, p_ipd, p_ipd_bea, t_ipd, ipd2pd]=...
    sezione_pivot(R, n_inc, inc, ...
    taglia, n_ref, tol_R, tol_angle)

% il parametro locale n_ref:
% corrisponde a n_ref_cyto quando sezione_pivot è chiamata da genemesh, che genera la mesh del problema omogeneizzato;
% corrisponde a n_ref_id quando sezione_pivot è chiamata da diffusione_id, che genera la mesh del problema di diffusione di E* sul disco inciso

% contatori
% n_pd è il numero di punti nel disco pivot
% n_int_pd è il numero di punti del disco pivot che non appartengono alla buccia né ad incisure
% n_tri è il numero di triangoli in cui è partizionato il disco pivot
% n_os è il numero di punti sul perimetro del disco pivot
% n_p_inc è un vettore che contiene il numero di punti della traccia di ogni incisura sul disco pivot

% corrispondenze
% os2pd contiene gli indici dei nodi del bordo del disco pivot, in senso antiorario, nella numerazione del disco pivot
% inc2pd è un cell array: per ogni incisura contiene gli indici dei nodi, dal vertice al bordo, nella numerazione del disco pivot
% pd2int_pd associa ad ogni nodo sul disco pivot: 
%   zero, se il nodo è sulla buccia oppure su una incisura; 
%   il numero progressivo nella numerazione dei soli nodi del disco pivot non sulla buccia né su incisure, altrimenti

% mesh
% p_pd è il vettore delle coordinate nodali del disco pivot
% t_pd è la matrice delle incidenze del disco pivot

% mesh del disco inciso (per il problema della diffusione di E* sul disco)
% n_ipd è il numero di punti nel disco pivot inciso
% p_ipd è il vettore delle coordinate nodali del disco pivot inciso
% p_ipd_bea è il vettore delle coordinate nodali del disco pivot inciso con incisure beanti
% t_ipd è la matrice delle incidenze del disco pivot inciso
% ipd2pd associa ad ogni nodo del disco inciso l'indice corrispondente nel disco non inciso


% genera la mesh non incisa
fun='incisure';
[p_pd,e_pd,t_pd]=initmesh(fun,'Hmax',taglia); 
% n_ref raffittimenti regolari
for j=1:n_ref
    [p_pd,e_pd,t_pd]=refinemesh(fun,p_pd,e_pd,t_pd,'regular'); 
end

% numero di nodi nel disco pivot
n_pd=size(p_pd,2);
% numero di triangoli nel disco pivot
n_tri=size(t_pd,2);

% nodi sul perimetro
% pesca tutti i nodi a distanza >=R-tol_R dall'origine
os2pd=find(abs(p_pd(1,:)+i*p_pd(2,:))>=R-tol_R);
% numero di nodi sul perimetro
n_os=length(os2pd);
% calcola le anomalie
angoli=angle(p_pd(1,os2pd)+i*p_pd(2,os2pd));
% trasporta da (-pi,pi] a [0,2*pi)
negativi=find(angoli<0-tol_angle);
angoli(negativi)=angoli(negativi)+2*pi;
% ordina nel senso delle anomalie crescenti
[th,I]=sort(angoli);
os2pd=os2pd(I);









% inizializza i parametri delle incisure, per evitare problemi nel caso non siano presenti incisure
n_p_inc=[];
inc2pd(1)={[]};
% ciclo sulle incisure
for k=1:n_inc	
	% indici dei nodi sulla incisura k, compreso quello di vertice
    % pesca tutti i nodi che hanno anomalia inc(2,k) pari a quella dell'incisura (mod 2*pi), a meno della tolleranza, 
    % e che al tempo stesso hanno distanza dall'origine >=R-lunghezza_incisura inc(1,k), a meno della 
    rho=abs(p_pd(1,:)+i*p_pd(2,:));
    theta=angle(p_pd(1,:)+i*p_pd(2,:));
    nodi=find((abs(modulo_near(theta-inc(2,k),2*pi))<tol_angle) & (rho>=R-inc(1,k)-tol_R));
    % ordina nel senso del raggio vettore cresecente
    [rr,I]=sort(abs(p_pd(1,nodi)+i*p_pd(2,nodi)));
    inc2pd(k)={nodi(I)};
    % numero dei nodi dell'incisura
    n_p_inc(k)=length(nodi);
end



% determina pd2int_pd
% insieme dei nodi sulla buccia o su una qualsiasi incisura
% sono quelli in cui la mesh del disco speciale è cucita a quella dell'os e delle incisure
np_cuciti=os2pd;
for k=1:n_inc
    np_cuciti=[np_cuciti, inc2pd{k}];
end
% il negato di ismember ritorna 0 in corrispondenza di ogni nodo cucito, ed 1 in corrispondenza di ogni nodo non cucito
interni=~ismember(1:n_pd,np_cuciti);
% numero di nodi non cuciti (tanti quanti sono gli 1 nel vettore interni)
n_int_pd=sum(interni);
% il cumsum quindi conta i soli nodi non cuciti perché si incrementa ad ogni 1
% il .* azzera i termini che corrispondono ai nodi cuciti
pd2int_pd=cumsum(interni).*interni;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCO INCISO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inizializza il numero di nodi nel disco inciso
n_ipd=n_pd;
% inizializza le coordinate correnti dei nodi del disco inciso
p_ipd=p_pd;
% inizializza le coordinate fittizie per il disegno delle incisure beanti
p_ipd_bea=p_pd;
% inizializza la matrice delle incidenze del disco inciso
t_ipd=t_pd;
% inizializza ipd2pd con l'identità, ristretta ai nodi del disco non inciso
ipd2pd= 1:n_pd;
% ciclo sulle incisure, che sdoppia i nodi appartenenti a ciascuna incisura, vertice escluso
for k=1:n_inc	
    % nello sdoppiamento devono essere aggiunti n_p_inc(k)-1 nodi (il vertice non è sdoppiato)
    vertice=inc2pd{k}(1);
    labbro=inc2pd{k}(2:end);

    % crea nuovi nodi, in fondo alla matrice p_ipd, duplicando quelli in labbro
    p_ipd=[p_ipd, p_ipd(:,labbro)];
    % nodi per le incisure beanti, ottenuti con una rotazione dell'angolo di beanza inc(3,k) intorno al vertice dell'incisura
    r=exp(i*inc(3,k))*((p_ipd(1,labbro)+i*p_ipd(2,labbro))-(p_ipd(1,vertice)+i*p_ipd(2,vertice)))+(p_ipd(1,vertice)+i*p_ipd(2,vertice));
    p_ipd_bea=[p_ipd_bea, [real(r);imag(r)]];
    
    % aggiorna ipd2pd: i nodi aggiunti in fondo mappano sui loro archetipi
    ipd2pd=[ipd2pd, labbro];

    % modifica la matrice delle incidenze
    % pesca gli elementi a sinistra dell'incisura k che hanno almeno un nodo su essa
    % individua gli elementi a sinistra dell'incisura k perché si trovano nella regione k+1 (cfr incisure.m)
    elem=find((t_ipd(4,:)==k+1) & (ismember(t_ipd(1,:),labbro)|ismember(t_ipd(2,:),labbro)|ismember(t_ipd(3,:),labbro)));

    % ciclo sui 3 nodi di ciascun elemento pescato
    for n=1:3
        % elementi a sinistra dell'incisura k che hanno il nodo n su essa
        elem_n=elem(ismember(t_ipd(n,elem),labbro));
        % la nuova numerazione dei nodi n-esimi di questi elementi viene ottenuta ricercando la loro posizione 
        % dentro la lista dei nodi sdoppiati, che è poi shiftata del valore corrente di n_ipd
        t_ipd(n,elem_n)=posizione(t_ipd(n,elem_n),labbro)+n_ipd;
    end

    % aggiorna il numero di nodi
    n_ipd=n_ipd+n_p_inc(k)-1;
end

% lascia in t_pd solo le incidenze
t_pd=t_pd(1:3,:);
t_ipd=t_ipd(1:3,:);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% routine per calcolare il modulo con approssimazione all'intero più vicino
function f=modulo_near(x,y)
    f=x-round(x/y)*y;        
return
