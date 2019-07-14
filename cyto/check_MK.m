% test per verificare l'assemblaggio delle matrici di massa e rigidezza
function check_MK(R, H, n_sez, n_sd, ...
    n_pd, n_tri, n_os, n_inc, n_p_inc, ...
    os2pd, inc2pd, ...
    p_pd, t_pd, Z_s, Z_sd, inc, ...
    M_vol, K_vol, M_sd, K_sd, M_os, K_os, M_inc, K_inc)

% dati sulle incisure
 l_inc=inc(1,:);

fprintf('\n\nErrori relativi\n\n');


% MESH VOL
% vanno commentati in assembla i termini di assemblaggio di sd, os, inc
% vanno commentati in questo file i check su sd, os, inc

% coordinate nodali: x e y in p_pd sono replicate in riga a blocchi per ogni sezione
% z in Z_s è replicata in riga elemento per elemento n_pd volte 
p_vol=[ reshape(p_pd(1,:)'*ones(1,n_sez),1,n_pd*n_sez);...
        reshape(p_pd(2,:)'*ones(1,n_sez),1,n_pd*n_sez);...
        reshape(ones(n_pd,1)*Z_s,1,n_pd*n_sez)];

% verifica dei volumi tramite matrice delle masse
a=pi*R^2*H;
vec=ones(n_pd*n_sez,1);
err=(vec'*M_vol*vec-a)/a;
fprintf('volume cilindro da M_vol  %8.4f\n',err)

% verifica dei momenti di inerzia tramite matrice delle masse
a=pi*R^4*H/4;
vec=p_vol(2,:)';
err=(vec'*M_vol*vec-a)/a;
fprintf('Ix cilindro da M_vol  %8.4f\n',err)

a=pi*R^4*H/4;
vec=p_vol(1,:)';
err=(vec'*M_vol*vec-a)/a;
fprintf('Iy cilindro da M_vol  %8.4f\n',err)

a=pi*R^2*H^3/3;
vec=p_vol(3,:)';
err=(vec'*M_vol*vec-a)/a;
fprintf('Iz cilindro da M_vol  %8.4f\n',err)

% Gradiente di una costante nel cilindro tramite matrice di rigidezza
vec=ones(n_pd*n_sez,1);
err=vec'*K_vol*vec;
fprintf('integrale del gradiente di una costante nel cilindro da K_vol  %8.4f\n',err)

% Vettore lineare con una delle tre coordinate
a=pi*R^2*H;
vec=p_vol(1,:)';
err=(vec'*K_vol*vec-a)/a;
fprintf('volume cilindro da K_vol  (f lineare in x)%8.4f\n',err)

a=pi*R^2*H;
vec=p_vol(2,:)';
err=(vec'*K_vol*vec-a)/a;
fprintf('volume cilindro da K_vol  (f lineare in y)%8.4f\n',err)

vec=p_vol(3,:)';
err=vec'*K_vol*vec;
fprintf('integrale del gradiente di una f lineare in z nel cilindro da K_vol  %8.4f\n',err)


% Vettore quadratico in una delle tre coordinate
a=pi*R^4*H/4;
vec=1/2*p_vol(2,:)'.^2;
err=(vec'*K_vol*vec-a)/a;
fprintf('Ix cilindro da K_vol  (f quadratica in y) %8.4f\n',err)

a=pi*R^4*H/4;
vec=1/2*p_vol(1,:)'.^2;
err=(vec'*K_vol*vec-a)/a;
fprintf('Iy cilindro da K_vol  (f quadratica in x) %8.4f\n',err)

vec=1/2*p_vol(3,:)'.^2;
err=vec'*K_vol*vec;
fprintf('integrale del gradiente_X di una f quadratica in z nel cilindro da K_vol  %8.4f\n',err)


% riga vuota
fprintf('\n')




% MESH SPECIAL DISK
% vanno commentati in assembla i termini di assemblaggio di vol, os, inc
% vanno commentati in questo file i check su vol, os, inc

for d=1:n_sd	
	% verifica dell'area tramite matrice delle masse
	a=pi*R^2;
	vec=ones(n_pd,1);
	err=(vec'*M_sd*vec-a)/a;
	fprintf('superficie del disco %i da M_sd  %8.4f\n',d,err)
	
    % verifica dei momenti di inerzia tramite matrice delle masse
	a=pi*R^4/4;
	vec=p_pd(2,:)';
	err=(vec'*M_sd*vec-a)/a;
	fprintf('Ix del disco %i da M_sd  %8.4f\n',d,err)
	
	a=pi*R^4/4;
	vec=p_pd(1,:)';
	err=(vec'*M_sd*vec-a)/a;
	fprintf('Iy del disco %i da M_sd  %8.4f\n',d,err)
	
	a=pi*R^2*Z_sd(d)^2;
	vec=Z_sd(d)*ones(n_pd,1);
	err=(vec'*M_sd*vec-a)/a;
	fprintf('Iz del disco %i da M_sd  %8.4f\n',d,err)
	

	% Gradiente di una costante nel disco tramite matrice di rigidezza
	vec=ones(n_pd,1);
	err=vec'*K_sd*vec;
	fprintf('integrale del gradiente di una costante nel disco %i da K_sd %8.4f\n',d,err)
	
	% Vettore lineare con una delle tre coordinate
	a=pi*R^2;
	vec=p_pd(1,:)';
	err=(vec'*K_sd*vec-a)/a;
	fprintf('superficie del disco %i da K_sd  (f lineare in x)%8.4f\n',d,err)
	
	a=pi*R^2;
	vec=p_pd(2,:)';
	err=(vec'*K_sd*vec-a)/a;
	fprintf('superficie del disco %i da K_sd  (f lineare in y)%8.4f\n',d,err)
	
	vec=Z_sd(d)*ones(n_pd,1);
	err=vec'*K_sd*vec;
	fprintf('integrale del gradiente_x di una f lineare in z nel disco %i da K_sd  %8.4f\n',d,err)

    % Vettore quadratico in una delle tre coordinate
	a=pi*R^4/4;
	vec=1/2*p_pd(2,:)'.^2;
	err=(vec'*K_sd*vec-a)/a;
	fprintf('Ix disco %i da K_sd  (f quadratica in y) %8.4f\n',d,err)
	
	a=pi*R^4/4;
	vec=1/2*p_pd(1,:)'.^2;
	err=(vec'*K_sd*vec-a)/a;
	fprintf('Iy disco %i da K_sd  (f quadratica in x) %8.4f\n',d,err)
	
	vec=1/2*Z_sd(d)*ones(n_pd,1).^2;
	err=vec'*K_sd*vec;
	fprintf('integrale del gradiente di una f quadratica in z nel disco %i da K_sd  %8.4f\n',d,err)

    % riga vuota
    fprintf('\n')

end




% MESH OUTER SHELL
% vanno commentati in assembla i termini di assemblaggio di vol, sd, inc
% % vanno commentati in questo file i check su vol, sd, inc

% coordinate nodali: x,y in p_os è ottenuta da os2pd, e replicata in riga a blocchi per ogni sezione
% z in Z_s è replicata in riga elemento per elemento n_os volte
% costruisce la mesh
% coordinate x,y,z
p_os=[ reshape(p_pd(1,os2pd)'*ones(1,n_sez),1,n_os*n_sez);...
        reshape(p_pd(2,os2pd)'*ones(1,n_sez),1,n_os*n_sez);...
        reshape(ones(n_os,1)*Z_s,1,n_os*n_sez)];

% verifica dell'area tramite matrice delle masse
a=2*pi*R*H;
vec=ones(n_os*n_sez,1);
err=(vec'*M_os*vec-a)/a;
fprintf('area outer shell da M_os  %8.4f\n',err)

% verifica dei momenti di inerzia tramite matrice delle masse
a=2*pi*R^3*H/2;
vec=p_os(2,:)';
err=(vec'*M_os*vec-a)/a;
fprintf('Ix outer shell da M_os  %8.4f\n',err)

a=2*pi*R^3*H/2;
vec=p_os(1,:)';
err=(vec'*M_os*vec-a)/a;
fprintf('Iy outer shell da M_os  %8.4f\n',err)

a=H^3*2*pi*R/3;
vec=p_os(3,:)';
err=(vec'*M_os*vec-a)/a;
fprintf('Iz outer shell da M_os  %8.4f\n',err)


% Gradiente di una costante nell'outer shell tramite matrice di rigidezza
vec=ones(n_os*n_sez,1);
err=vec'*K_os*vec;
fprintf('integrale del gradiente di una costante nell''ouetr da K_os  %8.4f\n',err)

% calcolo dell'integrale di \nabla_s x, \nabla_s y, \nabla_s z, sull'outer shell tramite matrice di rigidezza
a=R*H*pi;
vec=p_os(1,:)';
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s x|^2 sull''outer shell da K_os %8.4f\n',err)

a=R*H*pi;
vec=p_os(2,:)';
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s y|^2 sull''outer shell da K_os %8.4f\n',err)

a=2*pi*R*H;
vec=p_os(3,:)';
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s z|^2 sull''outer shell da K_os %8.4f\n',err)

% Vettore quadratico in una delle tre coordinate
a=pi*R^3*H/4;
vec=1/2*p_os(2,:)'.^2;
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s x^2|^2 sull''outer shell da K_os %8.4f\n',err)

a=pi*R^3*H/4;
vec=1/2*p_os(1,:)'.^2;
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s y^2|^2 sull''outer shell da K_os %8.4f\n',err)

a=H^3*2*pi*R/3;
vec=1/2*p_os(3,:)'.^2;
err=(vec'*K_os*vec-a)/a;
fprintf('int|nabla_s z^2|^2 sull''outer shell da K_os %8.4f\n',err)

% riga vuota
fprintf('\n')



% MESH INCISURE
% vanno commentati in assembla i termini di assemblaggio di vol, sd, os
% vanno commentati in questo file i check su vol, sd, os

p_inc=cell(1,n_inc);
for m=1:n_inc
    % costruisce la mesh
    nodi=inc2pd{m};
    nnodi=n_p_inc(m);
    % coordinate rho,z
	rho=abs(p_pd(1,nodi)+i*p_pd(2,nodi));
	p_inc(m)={[ reshape(rho'*ones(1,n_sez),1,nnodi*n_sez);...
                reshape(ones(nnodi,1)*Z_s,1,nnodi*n_sez)]};

	% integrale di 1*thickness(r_c) tramite la matrice delle masse
	a=H*l_inc(m)^2/2;
	vec=ones(n_p_inc(m)*n_sez,1);
    err=(vec'*M_inc{m}*vec-a)/a;
	fprintf('int 1*thickness(r_c) dr dz sull''incisura %i da M_inc{m} %8.4f\n',m,err)
	
	% integrale di r_c^2*thickness(r_c) tramite la matrice delle masse
	a=H*l_inc(m)^4/4;
	vec=p_inc{m}(1,:)'-rho(1);
	err=(vec'*M_inc{m}*vec-a)/a;
	fprintf('int r_c^2*thickness(r_c) dr dz sull''incisura %i da M_inc{m} %8.4f\n',m,err)
	
	% integrale di z^2*thickness(r_c) tramite la matrice delle masse
	a=H^3/3*l_inc(m)^2/2;
	vec=p_inc{m}(2,:)';
	err=(vec'*M_inc{m}*vec-a)/a;
	fprintf('int z^2*thickness(r_c) dr dz sull''incisura %i da M_inc{m} %8.4f\n',m,err)
	
	% Gradiente di una costante nell'incisura tramite matrice di rigidezza
	vec=ones(n_p_inc(m)*n_sez,1);
	err=vec'*K_inc{m}*vec;
	fprintf('integrale del gradiente di una costante nell''incisura %i da K_inc{m} %8.4f\n',m,err)
	
	% |nabla_B r_c|^2 nell'incisura tramite matrice di rigidezza
	a=H*l_inc(m)^2/2;
	vec=p_inc{m}(1,:)'-rho(1);
	err=(vec'*K_inc{m}*vec-a)/a;
	fprintf('|nabla_B r_c|^2 nell''incisura  %i da K_inc{m} %8.4f\n',m,err)

	% |nabla_B z|^2 nell'incisura tramite matrice di rigidezza
	a=H*l_inc(m)^2/2;
	vec=p_inc{m}(2,:)';
	err=(vec'*K_inc{m}*vec-a)/a;
	fprintf('|nabla_B z|^2 nell''incisura  %i da K_inc{m} %8.4f\n',m,err)
	
	% |nabla_B r_c^2|^2 nell'incisura tramite matrice di rigidezza
	a=4*H*l_inc(m)^4/4;
	vec=(p_inc{m}(1,:)'-rho(1)).^2;
	err=(vec'*K_inc{m}*vec-a)/a;
	fprintf('|nabla_B r_c^2|^2 nell''incisura  %i da K_inc{m} %8.4f\n',m,err)

	% |nabla_B z^2|^2 nell'incisura tramite matrice di rigidezza
	a=4*H^3/3*l_inc(m)^2/2;
	vec=p_inc{m}(2,:)'.^2;
	err=(vec'*K_inc{m}*vec-a)/a;
	fprintf('|nabla_B z^2|^2 nell''incisura  %i da K_inc{m} %8.4f\n',m,err)

    % riga vuota
	fprintf('\n')
	
end

return
