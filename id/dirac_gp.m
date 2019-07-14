% calcola i valori nodali corrispondenti alle singole delta, 
% la j-esima localizzata in p_delta(:,j)
function delta_nod=dirac_gp(n_p, p, t, M, p_delta, peak_delta, R, pt_min, pt_max)

% p, t definiscono la mesh
% M è la matrice delle masse (nota: non moltiplicata per alcun coefficiente)
% p_delta sono le locazioni delle delta, ciascuna con massa unitaria
%   queste devono essere all'interno del cerchio di raggio R
% peak_delta è una flag per la modalità di costruzione della delta
%   vero: delta del tipo M\V;
%   falso: delta del tipo k*V
% pt_min, pt_max contengno informazioni sui rettangoli minimi contenenti i
%   triangoli, e servono per rendere più rapida la chiamata a find_tri

% la routine calcola il vettore V dei carichi corrispondente a ciascuna delle delta, 
% campionando in essa le funzioni di forma dei triangoli che le contengono
% quindi V ha per righe il numero di nodi, per colonne il numero delle delta , 
% cioè il numero di di colonne di p_delta

% se peak_delta=vero, 
% i valori nodali sono calcolati, in modo esatto, come M\V
% determina i valori nodali di una funzione lineare a tratti cui corrisponde lo stesso carico V
% questo approccio può dare valori negativi
% è preferifile, per il fatto che conserva il carico, quando le delta non
% diffondono

% se peak_delta=falso,
% i valori nodali sono presi come k*V, dove k è determinata dalla
% condizione di conservazione della massa
% questi valori nodali corrispondono ad una pagoda
% cui però corrisponde un vettore dei carichi in generale diverso da V
% questo approccio è preferifile per il calcolo del dato iniziale quando
% le delta diffondono

% numero delle delta
n_delta=size(p_delta,2);

% assembla il vettore V dei carichi corrispondente a ciascuna delta, 
% campionando in essa le funzioni di forma dei triangoli che le contengono
% per ogni locazione di fotoisomerizzazione, determina l'elemento 
% in cui essa si trova ed il corrispondente vettore di interpolazione F_elem
% determina gli elementi tri della mesh di triangoli 
% definita da p_ipd, t_ipd appartengono i punti p_delta
% e campiona in essi le corrispondenti funzioni di forma
[tri, F_elem]=find_tri(p, t, p_delta, R, pt_min, pt_max);



% calcola il vettore dei termini noti
% t(:,tri) contiene, per ciascuna locazione di fotoisomerizzazione,
% i tre nodi del triangolo che contiene tale punto
% trasforma in vettore, incolonnando i tre nodi di ciascun triangolo con
% quelli degli altri triangoli
nodi_tri=t(:,tri);
% F_elem(:) contiene, pper ciascuna locazione di fotoisomerizzazione,
% i valori delle funzioni di forma in esso, relative al triangolo che contiene tale punto
% mette il contributo di ciascuna delta in una colonna diversa
% questo è necessario perché le delta possono essere ciascuna in uno stato
% diverso, quindi i loro contributi vanno tenuti separati
colonne=repmat(1:n_delta,3,1);
V=accumarray([nodi_tri(:),colonne(:)],F_elem(:),[n_p,n_delta]);


if peak_delta
    % delta esatta: può dare valori negativi
    % determina i valori nodali di una funzione lineare a tratti cui corrisponde
    % lo stesso vettore V calcolato a partire dalle delta
    delta_nod=M\V;
else
    % approccio alternativo: delta a pagoda, aprossimato
    % delta_nod è un vettore proporzionale a V, in modo tale che ad esso
    % corrisponda massa unitaria 
    % (la massa corrispondente ad un insieme V di valori nodali è ones(1,n_ipd)*M*V)
    % questo approccio non fornisce densità negative
    delta_nod=V./repmat((ones(1,n_p)*M*V),n_p,1);
end


return
