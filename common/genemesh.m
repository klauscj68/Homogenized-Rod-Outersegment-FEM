% costruisce la mesh
function [n_pd, n_int_pd, n_tri, n_os, n_p_inc, ...
         os2pd, inc2pd, pd2int_pd, p_pd, t_pd, ...
         n_ipd, p_ipd, p_ipd_bea, t_ipd, ipd2pd, ...
         n_pd_id, n_int_pd_id, n_tri_id, n_os_id, n_p_inc_id, ...
         os2pd_id, inc2pd_id, pd2int_pd_id, p_pd_id, t_pd_id, ...
         n_ipd_id, p_ipd_id, p_ipd_bea_id, t_ipd_id, ipd2pd_id, ...
         Z_s, sd2sez]=...
     genemesh(R, H, n_sez, flag_geom_sp, dz_0, ...
     n_sd, Z_sd, ...
     n_inc, inc, ...
     taglia, tol_R, tol_angle, n_ref_cyto, n_ref_id, ...
     plot_mesh, plot_num, inspect)
     
% output

% contatori
% n_pd è il numero di punti nel disco pivot
% n_int_pd è il numero di punti del disco pivot che non appartengono alla circonferenza esterna né ad incisure
% n_tri è il numero di triangoli in cui è partizionato il disco pivot
% n_os è il numero di punti sul perimetro del disco pivot
% n_p_inc è un vettore che contiene il numero di punti della traccia di ogni incisura sul disco pivot
% n_sez è il numero delle sezioni nel volume (assegnato in input)
% n_sd è il numero dei dischi speciali (assegnato in input)

% corrispondenze
% os2pd contiene gli indici dei nodi del bordo del disco pivot, in senso antiorario, nella numerazione del disco pivot
% inc2pd è un cell array: per ogni incisura contiene gli indici dei nodi, dal vertice al bordo, nella numerazione del disco pivot
% pd2int_pd associa ad ogni nodo sul disco pivot: 
%   zero, se il nodo è sulla buccia oppure su una incisura; 
%   il numero progressivo nella numerazione dei soli nodi del disco pivot non sulla buccia né su incisure, altrimenti
% sd2sez associa a ciascun disco speciale il numero di sezione cui esso è sovrapposto

% mesh
% p_pd è il vettore delle coordinate nodali del disco pivot
% t_pd è la matrice delle incidenze del disco pivot
% t_pd è la matrice delle incidenze del disco pivot inciso
% Z_s è il vettore delle quote delle sezioni rette

% a soli fini grafico-estetici in questa routine sono calcolati:
% p_vol, t_vol è la mesh di prismi retti a base triangolare nel volume
% p_sd, t_sd è la mesh di triangoli su uno qualsiasi dei dischi speciali
% p_os, t_os è la mesh di rettangoli sull'outer shell
% p_inc, t_inc sono strutture a cell array che contengono le mesh di rettangoli su ciascuna incisura
% nota: queste variabili non escono da genemesh, e non sono usate altrove

% nota bene:
% le variabili 
% n_pd, n_int_pd, n_tri, n_os, n_p_inc, ...
%          os2pd, inc2pd, pd2int_pd, p_pd, t_pd, ...
%          n_ipd, p_ipd, p_ipd_bea, t_ipd, ipd2pd
% escono anche nella forma 
%          n_pd_id, n_int_pd_id, n_tri_id, n_os_id, n_p_inc_id, ...
%          os2pd_id, inc2pd_id, pd2int_pd_id, p_pd_id, t_pd_id, ...
%          n_ipd_id, p_ipd_id, p_ipd_bea_id, t_ipd_id, ipd2pd_id
% che si riferisce alla mesh del disco con grado di raffittimento n_ref_id
% anziché n_ref_cyto
% le variabili senza suffisso _id servono per il problema della diffusione
% nel cytosol; le variabili col suffisso _id servono per il problema 
% della cascata sul disco


% numero di segmenti in cui è partizionato [0,H]
n_seg=n_sez-1;

% mesh del disco pivot, con e senza incisure
% nota bene: qui si passa n_ref_cyto, perché questo output serve per la mesh
% della diffusione hom
fprintf('\nGenera la mesh del disco inciso per il problema nel cytosol (raffittimento %4i)\n',n_ref_cyto);
[n_pd, n_int_pd, n_tri, n_os, n_p_inc, ...
        os2pd, inc2pd, pd2int_pd, p_pd, t_pd, ...
        n_ipd, p_ipd, p_ipd_bea, t_ipd, ipd2pd]=...
    sezione_pivot(R, n_inc, inc, ...
    taglia, n_ref_cyto, tol_R, tol_angle);





% informativa
fprintf('\nDati sulla mesh della diffusione hom');
fprintf('\nMesh della sezione trasversale: %4i nodi; %4i triangoli', n_pd, n_tri);
fprintf('\nMesh del volume:                %4i nodi; %4i prismi', n_pd*n_sez, n_tri*n_seg);
fprintf('\nMesh del outer shell:           %4i nodi; %4i rettangoli', n_os*n_sez, n_os*n_seg);
for k=1:n_inc
    fprintf('\nMesh dell''incisura %2i:          %4i nodi; %4i rettangoli', k, n_p_inc(k)*n_sez, (n_p_inc(k)-1)*n_seg);
end
fprintf('\n');

% % stampe di controllo
% n_pd
% n_int_pd
% n_tri
% n_os
% os2pd
% pd2int_pd
% p_pd
% t_pd, 
% n_ipd
% p_ipd
% p_ipd_bea
% t_ipd
% ipd2pd
% for k=1:n_inc
%     'incisura', k
%     'n_p_inc(k)', n_p_inc(k)
%     'inc2pd{k}', inc2pd{k}
% end

% crea la mesh del disco inciso con grado di raffittimento n_ref_id
% infatti il numero di raffittimenti del problema della diffusione di E^* 
% sul disco inciso può essere maggiore di quello del problema omogeneizzato
fprintf('\nGenera la mesh del disco inciso per il problema sul disco (raffittimento %4i)\n',n_ref_id);
[n_pd_id, n_int_pd_id, n_tri_id, n_os_id, n_p_inc_id, ...
        os2pd_id, inc2pd_id, pd2int_pd_id, p_pd_id, t_pd_id, ...
        n_ipd_id, p_ipd_id, p_ipd_bea_id, t_ipd_id, ipd2pd_id]=...
    sezione_pivot(R, n_inc, inc, ...
    taglia, n_ref_id, tol_R, tol_angle);



% informativa
fprintf('\nMesh del disco inciso: %i nodi; %i triangoli\n', n_ipd_id, n_tri_id);

% calcola le quote delle sezioni
fprintf('\nDetermina le quote delle sezioni\n');
[Z_s, sd2sez]=quote(H, n_sez, flag_geom_sp, dz_0, n_sd, Z_sd);






% sd2sez
% Z_sd
% 'Z_s(sd2sez)', Z_s(sd2sez)

% disegna la mesh del disco pivot non inciso
% con la traccia delle incisure
if plot_mesh
	figure(1)
    hold on
	pdemesh(p_pd,[],t_pd), axis equal
    for k=1:n_inc
        plot(p_pd(1,inc2pd{k}([1,end])),p_pd(2,inc2pd{k}([1,end])),'r','Linewidth',1.5)
    end
    if plot_num
		for cont=1:n_pd
            aa=text(p_pd(1,cont),p_pd(2,cont),num2str(cont));
            set(aa,'color','k');
		end
		for cont=1:n_tri
            bb=text(mean(p_pd(1,t_pd(1:3,cont))),mean(p_pd(2,t_pd(1:3,cont))),num2str(cont));
            set(bb,'color','g');
		end
    end
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
    axis([-R*1.1 R*1.1 -R*1.1 R*1.1])
	view(2)
end




% disegna la mesh dei dischi speciali con incisure beanti
if plot_mesh
    for d=1:n_sd
		figure(2)
        hold on
		% Mesh di triangoli nello special disk d
        vert=[p_ipd_bea; ones(1,n_ipd)*Z_sd(d)];
		l=patch('Vertices',vert','Faces',t_ipd');
		set(l,'facecolor','r')
		axis equal
		xlabel('x [\mu m]')
		ylabel('y [\mu m]')
		zlabel('z [\mu m]')
        axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 H])
		view(3)
    end
end

% mesh di prismi nel volume
% è ottenuta come prodotto cartesiano tra il disco pivot e il vettore delle quote

% coordinate nodali: x e y in p_ipd sono replicate in riga a blocchi per ogni sezione
% z in Z_s è replicata in riga elemento per elemento n_ipd volte 
% coordinate con incisure beanti, a fini estetici
p_vol_bea=[ reshape(p_ipd_bea(1,:)'*ones(1,n_sez),1,n_ipd*n_sez);...
        reshape(p_ipd_bea(2,:)'*ones(1,n_sez),1,n_ipd*n_sez);...
        reshape(ones(n_ipd,1)*Z_s,1,n_ipd*n_sez)];

% incidenze nel volume
t_vol=zeros(6,0);
for k=1:n_sez-1
    t_vol=[t_vol, [(k-1)*n_ipd+t_ipd; k*n_ipd+t_ipd]];
end

% disegna la mesh dei prismi
if plot_mesh
	figure(3)
	% Mesh di prismi nel cilindro
    % basi
	l=patch('Vertices',p_vol_bea','Faces',[t_vol([1 2 3],:)'; t_vol([4 5 6],:)']);
    % superficie laterale
	l=patch('Vertices',p_vol_bea','Faces',[t_vol([1 2 5 4],:)'; t_vol([2 3 6 5],:)'; t_vol([3 1 4 6],:)']);
	set(l,'facecolor','c')
	axis equal
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
	zlabel('z [\mu m]')
    axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 H])
	view(3)
end

% mesh di rettangoli nell'outer shell
% sono ottenute come prodotto cartesiano tra la buccia del disco pivot e il vettore delle quote

% coordinate nodali: x,y in p_os è ottenuta da os2pd, e replicata in riga a blocchi per ogni sezione
% z in Z_s è replicata in riga elemento per elemento n_os volte
% costruisce la mesh
% coordinate x,y,z
p_os=[ reshape(p_pd(1,os2pd)'*ones(1,n_sez),1,n_os*n_sez);...
        reshape(p_pd(2,os2pd)'*ones(1,n_sez),1,n_os*n_sez);...
        reshape(ones(n_os,1)*Z_s,1,n_os*n_sez)];
% incidenze
t_os=zeros(4,0);
for k=1:n_sez-1
    t_os=[t_os, [(k-1)*n_os+(1:n_os); (k-1)*n_os+[2:n_os, 1]; k*n_os+(1:n_os); k*n_os+[2:n_os, 1]] ];
end

% disegna la mesh dei settori sull'os
if plot_mesh
	figure(4)
	% Mesh di rettangoli nell'outer shell
	l=patch('Vertices',p_os','Faces',t_os([1 2 4 3],:)');
	set(l,'facecolor','g')
	axis equal
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
	zlabel('z [\mu m]')
    axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 H])
	view(3)
end

% mesh di rettangoli nelle incisure
% sono ottenute come prodotto cartesiano tra la traccia dell'incisura sul disco pivot e il vettore delle quote

% coordinate nodali: rho in p_inc è ottenuta inc2pd, e replicata in riga a blocchi per ogni sezione
% z in Z_s è replicata in riga elemento per elemento n_p_inc volte
p_inc=cell(1,n_inc);
t_inc=cell(1,n_inc);
for j=1:n_inc
    % costruisce la mesh
    nodi=inc2pd{j};
    n_nodi=n_p_inc(j);
    % coordinate rho,z
	rho=abs(p_ipd(1,nodi)+i*p_ipd(2,nodi));
	p_inc(j)={[ reshape(rho'*ones(1,n_sez),1,n_nodi*n_sez);...
                reshape(ones(n_nodi,1)*Z_s,1,n_nodi*n_sez)]};
    % incidenze nell'incisura
    t_inc(j)={zeros(4,0)};
	for k=1:n_sez-1
        t_inc(j)={[t_inc{j}, [(k-1)*n_nodi+(1:n_nodi-1); (k-1)*n_nodi+(2:n_nodi); k*n_nodi+(1:n_nodi-1); k*n_nodi+(2:n_nodi)] ]};
	end
end

% disegna la mesh delle incisure
if plot_mesh
    for j=1:n_inc
		figure(5)
        hold on
		% Mesh di rettangoli nella incisura j, di anomalia inc(2,j)
        vert=[p_inc{j}(1,:)*cos(inc(2,j)); p_inc{j}(1,:)*sin(inc(2,j)); p_inc{j}(2,:)]; 
		l=patch('Vertices',vert','Faces',t_inc{j}([1 2 4 3],:)');
		set(l,'facecolor','g')
		axis equal
		xlabel('x [\mu m]')
		ylabel('y [\mu m]')
		zlabel('z [\mu m]')
        axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 H])
		view(3)
    end
end

% disegna la mesh incisa, con incisure beanti
if plot_mesh
	figure(6)
    hold on
	pdemesh(p_ipd_bea_id,[],t_ipd_id), axis equal
    for k=1:n_inc
        plot(p_pd_id(1,inc2pd_id{k}([1,end])),p_pd_id(2,inc2pd{k}([1,end])),'r','Linewidth',1.5)
    end
    if plot_num
        for cont=1:n_ipd_id
            aa=text(p_ipd_bea_id(1,cont),p_ipd_bea_id(2,cont),num2str(cont));
            set(aa,'color','k');
        end
        % il numero di elementi non varia nell'effettuare le incisioni
        for cont=1:n_tri
            bb=text(mean(p_ipd_bea_id(1,t_ipd(1:3,cont))),mean(p_ipd_bea_id(2,t_ipd(1:3,cont))),num2str(cont));
            set(bb,'color','g');
        end
    end
	xlabel('x [\mu m]')
	ylabel('y [\mu m]')
    axis([-R*1.1 R*1.1 -R*1.1 R*1.1])
	view(2)
end


presskey(inspect);




return
