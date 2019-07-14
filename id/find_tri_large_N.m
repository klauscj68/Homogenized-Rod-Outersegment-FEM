% determina gli elementi tri(i) della mesh di triangoli definita da p,t cui p_find(:,i) appartiene
% e campiona in p_find(:,i) le corrispondenti funzioni di forma F(:,i)
% questo codice è vettorializzato sui punti
% nota: se i punti da piazzare sono pochi, conviene usare
% find_tri, che è vettorializzato sui triangoli
function [tri, F, pt_min, pt_max]=find_tri_large_N(p, t, p_find, R, pt_min, pt_max)

% le matrici pt_min e pt_max determinano i rettangoli minimi con
% lati paralleli agli assi cartesiani, contenenti ciascun triangolo
% pt_min(:,j) è il vertice in basso a sinistra del rettangolo minimo relativo al triangolo j
% pt_max(:,j) è il vertice in alto a destra del rettangolo minimo relativo al triangolo j

% il calcolo di queste matrici è alquanto dispendioso
% però dipende dalla mesh e non da p_find
% quindi può essere eseguito una volta soltanto, passandole tra i dati le
% volte successive

% p_find è una matrice con 2 colonne; p_find(:,i) contiene le coordinate
% del punto da processare

% metodo:
% ogni triangolo è pensato contenuto in un rettangolo minimo con
% lati paralleli agli assi cartesiani
% viene eseguito un ciclo sui triangoli
% si fa prima una sgrossata, pescando i punti che cadono nel suo rettangolo
% minimo contenente, poi si lavora di fino coi semipiani

% tolleranza utilizzata
tol=10000*eps;

% tolleranza aggiuntiva sul bounding box, dovuta al fatto che l'unione dei
% triangoli è strettamente contenuta nel cerchio
tol_pap=tol+max([R-max(p(1,:)), R+min(p(1,:)), R-max(p(2,:)), R+min(p(2,:))]);

% determina se deve calcolarsi pt_min e pt_max
if nargin==4,

    % calcola pt_min e pt_max

    % le coordinate dei vertici di ciascun triangolo sono:
    % primi vertici: p(:,t(1,:))
    % secondi vertici: p(:,t(2,:))
    % terzi vertici: p(:,t(3,:));

    % determina i rettangoli minimi contenenti ciascun triangolo
    % minimo di ciascuna coordinata
    pt_min=zeros(2,size(t,2));
    pt_min(1,:)=min([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[],1);
    pt_min(2,:)=min([p(2,t(1,:));p(2,t(2,:));p(2,t(3,:))],[],1);
    % massimo di ciascuna coordinata
    pt_max=zeros(2,size(t,2));
    pt_max(1,:)=max([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[],1);
    pt_max(2,:)=max([p(2,t(1,:));p(2,t(2,:));p(2,t(3,:))],[],1);
end

% inizializza gli array di output
n_p_find=size(p_find,2);
tri=zeros(1,n_p_find);
rett=zeros(1,n_p_find);
F=zeros(3,n_p_find);

% ciclo su tutti triangoli
for cont=1:size(t,2)

    % individua i punti contenuti nel rettangolo minimo contenente il
    % triangolo cont

%     manipolo=find( (p_find(1,:)>pt_min(1,cont)-tol_pap)&(p_find(1,:)<pt_max(1,cont)+tol_pap) & ...
%                    (p_find(2,:)>pt_min(2,cont)-tol_pap)&(p_find(2,:)<pt_max(2,cont)+tol_pap) );

    % metodo short circuit
    manipolo=find(     (p_find(1,:)       >pt_min(1,cont)-tol_pap) );
    manipolo=manipolo( (p_find(1,manipolo)<pt_max(1,cont)+tol_pap) );
    manipolo=manipolo( (p_find(2,manipolo)>pt_min(2,cont)-tol_pap) );
    manipolo=manipolo( (p_find(2,manipolo)<pt_max(2,cont)+tol_pap) );

    % se manipolo è vuoto, passa oltre
    if ~isempty(manipolo)
        % attribuisce il rettangolo cont a tutti i punti di manipolo
        % questo vettore servirà per assegnare qualcosa ai punti che sono
        % rimasti senza triangolo (cioè che hanno il corrispondente tri a zero 
        % alla fine del ciclo sui triangoli)
        rett(manipolo)=cont;

        % individua i punti che sono veramente contenuti nel triangolo cont
        % (non solo nel rettangolo minimo che lo contiene)
        % il triangolo è orientato, quindi se (i,j,k) è una permutazione positiva,
        % (pj-pi)\cross(pk-pi))\dot verosre k è positivo
        % quindi p appartiene al triangolo se 
        % (pj-pi)\cross(p-pi)\dot versore k è positivo
        % cioè
        % \cross(p-pi)\dot [versore k \cross(pj-pi) ] è positivo
        % nota: [versore k \cross(pj-pi) ] è l'hodge di (pj-pi)
        % per le tre permutazioni positive
        % 1,2,3
        % 2,3,1
        % 3,1,2
        % primo vertice
        X1=p(1,t(1,cont));
        Y1=p(2,t(1,cont));
        % secondo vertice
        X2=p(1,t(2,cont));
        Y2=p(2,t(2,cont));
        % terzo vertice
        X3=p(1,t(3,cont));
        Y3=p(2,t(3,cont));
        % punti papabili: sono entro il rettangolo minimo contenente il
        % triangolo cont
        Xp=p_find(1,manipolo);
        Yp=p_find(2,manipolo);
        % ricerca quali dei punti di manipolo sono effettivamente nel triangolo cont
        % e riferisce alla numerazione in manipolo, che è quella globale
        manipolo_eff=manipolo( ((Y2-Y1).*(Xp-X1)-(X2-X1).*(Yp-Y1)<=tol) & ((Y3-Y2).*(Xp-X2)-(X3-X2).*(Yp-Y2)<=tol) & ((Y1-Y3).*(Xp-X3)-(X1-X3).*(Yp-Y3)<=tol) );
        if ~isempty(manipolo_eff)
            % attribuisce il triangolo cont a tutti i punti di manipolo_eff
            tri(manipolo_eff)=cont;
        end
    end
end

% verifica se ci sono punti che non hanno avuto assegnato il loro triangolo
orfani=find(tri==0);
tri(orfani)=rett(orfani);

% ciclo sui triangoli per determinare il vettore F che tiene conto di tutti
% i punti nel triangolo 
for cont=1:size(t,2)
    % posizione dei nodi del triangolo cont
    X=p(1,t(:,cont));
    Y=p(2,t(:,cont));
    % indice dei punti attribuiti al triangolo cont
    point=find(tri==cont);
    % posizione dei punti contenuti in tale triangolo
    Xp=p_find(1,point);
    Yp=p_find(2,point);
    area=-X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(1)*Y(2) + X(3)*Y(1) - X(3)*Y(2);
    F(:,point) = [X(2)*Y(3)-X(3)*Y(2)-Y(3)*Xp+Yp*X(3)+Xp*Y(2)-X(2)*Yp; ...
        -Xp*Y(1)+Y(3)*Xp-X(1)*Y(3)+X(1)*Yp-Yp*X(3)+X(3)*Y(1); ...
        Xp*Y(1)-Xp*Y(2)+X(1)*Y(2)+X(2)*Yp-X(2)*Y(1)-X(1)*Yp]/area;
end
