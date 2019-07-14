function [n_gdl, gdl_vol, gdl_sd, gdl_os, gdl_inc]=...
    gdl(n_pd, n_int_pd, n_os, n_inc, n_p_inc, ...
    n_sez, n_sd, ...
    os2pd, inc2pd, pd2int_pd, sd2sez, flag_model, inspect)

% gradi di libertà corrispondenti ai nodi del volume, dei dischi speciali, dell'outer shell, delle incisure

% ordine dei gradi di libertà nelle matrici globali:
% modello 3D (flag_model=1)
%   dapprima i nodi del volume, pari al numero delle sezioni per i nodi su un disco pivot (n_pd)
%   poi i nodi non cuciti dei dischi speciali, 
%   pari al numero dei dischi speciali per i nodi del disco pivot che non appartengono alla buccia né ad incisure
%   all'os e alle incisure non corrispondono gdl aggiuntivi, perché la soluzione sull'os e sulle incisure 
%   è ottenuta come traccia di quella nel volume
% modello well-stirred trasversale (flag_model=2)
%   un gdl per tutti i nodi di ciascuna sezione;
%   i gdl sono numerati nel senso delle z crescenti
% modello well-stirred globale (flag_model=3)
%   un solo gdl per tutti i nodi

% contatori
% n_pd è il numero di punti nel disco pivot
% n_int_pd è il numero di punti del disco pivot che non appartengono alla circonferenza esterna né ad incisure
% n_tri è il numero di triangoli in cui è partizionato il disco pivot
% n_os è il numero di punti sul perimetro del disco pivot
% n_p_inc è un vettore che contiene il numero di punti della traccia di ogni incisura sul disco pivot
% n_sez è il numero delle sezioni nel volume
% n_sd è il numero dei dischi speciali

% corrispondenze
% os2pd contiene gli indici dei nodi del bordo del disco pivot, in senso antiorario, nella numerazione del disco pivot
% inc2pd è un cell array: per ogni incisura contiene gli indici dei nodi, dal vertice al bordo, nella numerazione del disco pivot
% pd2int_pd associa ad ogni nodo sul disco pivot: 
%   zero, se il nodo è sulla buccia oppure su una incisura; 
%   il numero progressivo nella numerazione dei soli nodi del disco pivot non sulla buccia né su incisure, altrimenti
% sd2sez associa a ciascun disco speciale il numero di sezione cui esso è sovrapposto


% output
% n_gdl è il numero di gradi di libertà per ogni variabile (u e v)
% il numero totale dei gdl del problema è il doppio di n_gdl
% gdl_vol(j,k)    è il gdl del nodo di indice j nella numerazione del disco pivot, nella sezione k
% gdl_sd(j,k)     è il gdl del nodo di indice j nella numerazione del disco pivot, nel disco speciale k
% gdl_os(j,k)     è il gdl del nodo di indice j nell'outer shell (nella numerazione progressiva della buccia), nella sezione k
% gdl_inc{m}(j,k) è il gdl del nodo nell'incisura m, j^mo a partire dal vertice, nella sezione k

fprintf('\nDetermina la corrispondenza fra nodi e gradi di libertà\n');


switch flag_model
    case 1
        % modello 3D

        % contatore progrssivo dei gdl
        n_gdl=0;

        % gdl dei nodi del volume
        gdl_vol=n_pd*repmat(0:n_sez-1,n_pd,1)+ repmat((1:n_pd)',1,n_sez);
        % aggiorna i gdl
        n_gdl=n_gdl+n_sez*n_pd;

        % gdl dei nodi dei dischi speciali
        gdl_sd=zeros(n_pd,n_sd);
        for k=1:n_sd
            % per i nodi sulla buccia o sulle incisure, in cui pd2int_pd è zero, prende il gdl del corrispondente nodo della mesh del volume
            % per i nodi interni, in cui pd2int_pd non è zero, crea un nuovo gdl che legge da pd2int_pd
            % gdl che si avrebbe se tutti i nodi coincidessero con quelli del volume (cucitura)
            gdl_c=n_pd*(sd2sez(k)-1)+ (1:n_pd);
            % gdl che deriva dall'essere, i nodi interni, nuovi rispetto a quelli del volume
            gdl_nc=n_gdl+pd2int_pd;
            % combina in modo opportuno:
            % nelle posizioni in cui pd2int_pd è zero (nodo cucito), prende il gdl da gdl_c;
            % nelle posizioni in cui pd2int_pd è non zero (nodo nuovo), prende il gdl da gdl_nc;
            gdl_sd(:,k)=(gdl_c.*(pd2int_pd==0)+gdl_nc.*(pd2int_pd~=0))';
            % aggiorna n_gdl, sommando quelli appena aggiunti
            n_gdl=n_gdl+n_int_pd;
        end

        % gdl dei nodi nell'outer shell
        gdl_os=n_pd*repmat(0:n_sez-1,n_os,1)+repmat(os2pd',1,n_sez);
        % non ci sono gdl aggiuntivi

        % gdl dei nodi delle incisure
        % inizializza gdl_inc per evitare l'errore nel caso in cui non ci siano incisure
        gdl_inc(1)={[]};
        for m=1:n_inc
            gdl_inc(m)={n_pd*repmat(0:n_sez-1,n_p_inc(m),1)+repmat(inc2pd{m}',1,n_sez)};
        end
        % non ci sono gdl aggiuntivi
        
        
    case 2
        % modello well-stirred trasversale

        % i gdl sono pari al numero delle sezioni
        n_gdl=n_sez;
    
        % i nodi di una sezione corrispondono allo stesso gdl
        gdl_vol=repmat(1:n_sez,n_pd,1);

        % gdl dei nodi dei dischi speciali
        % i nodi di ciascun disco speciale hanno il gdl
        % corrispondente ai nodi sulla circonferenza in comune con l'outer shell
        % che a sua volta è quello dei nodi dell'interior alla stessa quota
        gdl_sd=repmat(sd2sez,n_pd,1);

        % gdl dei nodi nell'outer shell
        gdl_os=repmat(1:n_sez,n_os,1);

        % gdl dei nodi delle incisure
        % inizializza gdl_inc per evitare l'errore nel caso in cui non ci siano incisure
        gdl_inc(1)={[]};
        for m=1:n_inc
            gdl_inc(m)={repmat(1:n_sez,n_p_inc(m),1)};
        end


    case 3
        % modello well-stirred globale

        % un solo gdl per tutti i nodi
        n_gdl=1;

        % volume
        gdl_vol=ones(n_pd,n_sez);
        
        % dischi speciali
        gdl_sd=ones(n_pd,n_sd);
        
        % outer shell
        gdl_os=ones(n_os,n_sez);

        % incisure
        % inizializza gdl_inc per evitare l'errore nel caso in cui non ci siano incisure
        gdl_inc(1)={[]};
        for m=1:n_inc
            gdl_inc(m)={ones(n_p_inc(m),n_sez)};
        end
        
end % switch

% % stampe di controllo
% n_gdl, gdl_vol, gdl_sd, gdl_os
% for k=1:n_inc
%     'incisura',k
%     gdl_inc{k}
% end

% informativa
fprintf('\nMesh generale: %i gradi di libertà complessivi (u e v)\n', 2*n_gdl);
presskey(inspect);




return
