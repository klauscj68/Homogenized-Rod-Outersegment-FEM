function [x,y]=incisure(bs,s) 
% disco con incisure

% legge i dati sulla geometria
[R, H, n_sez, flag_geom_sp, dz_0, ...
        n_sd, Z_sd, ...
        n_inc, inc]=data;
    
% estrae da inc le informazioni sulle incisure
% lunghezze
l_inc=inc(1,:);
% posizioni angolari
theta_inc=inc(2,:);

% nel caso ci fossero meno di tre incisure, ne aggiunge di fittizie al fine di ottenerne tre
if n_inc==0
    % tre incisure a 0,2*pi/3,4*pi/3, ciascuna lunga R/3
    n_inc=3;
    l_inc=R/3*ones(1,3);
    theta_inc=2*pi/3*[0, 1, 2];
elseif n_inc==1
    % due incisure aggiuntive a 2*pi/3 gradi dall'unica presente, lunghe quanto questa
    n_inc=3;
    l_inc=l_inc(1)*ones(1,3);
    theta_inc=theta_inc(1)+2*pi/3*[0, 1, 2];
elseif n_inc==2
    % una incisura aggiuntiva che biseca quelle presenti
    % controlla che l'angolo dalla prima alla seconda sia minore del piatto, 
    % sicché l'incisura aggiunta, a pi+la media delle posizioni angolari delle altre, biseca l'angolo maggiore di pi
    if theta_inc(2)-theta_inc(1)>pi
        disp([theta_inc(1),theta_inc(2)]);
        error('l''angolo dalla prima alla seconda incisura, in senso antiorario, deve essere minore di pi')
    end
    n_inc=3;
    l_inc=[l_inc, mean(l_inc)];
    theta_inc=[theta_inc, mean(theta_inc)+pi];
end

% Calcolo numero di bordi, costituiti dai lati del poligono centrale, le incisure, e gli archi esterni
nbs=3*n_inc;

% Numero di segmenti
n_seg=2*n_inc;

% Numero lati poligono centrale
n_pol=n_inc;

% ritorna nbs
if nargin==0  
  x=nbs;   
  return 
end

% Calcolo della matrice dl
% sono numerati per primi i lati del poligono centrale, poi i segmenti-incisura, infine gli archi di circonferenza
dl=zeros(4,nbs);

% Tutti i segmenti si parametrizzano con parametro variabile tra 0 e 1
dl(1:2,1:n_seg)=[zeros(1,n_seg); ones(1,n_seg)];

% Tutti gli archi si parametrizzano come porzioni di una circonferenza da 0 a 2*pi
dl(1:2,n_seg+1:nbs)=[theta_inc; [theta_inc(2:end), theta_inc(1)+2*pi]];

% Assegnazione delle etichette per le zone delimitate dalle linee di bordo
% poligono centrale = zona 1
% settore fra l'incisura i e quella i+1 = zona i+1

% lati del poligono centrale = a sin zona 1, a dx zona i+1
dl(3:4,1:n_pol)=[ones(1,n_pol); 2:n_pol+1];

% incisure = a sin zona i+1, a destra zona i, salvo la prima incisura, che ha a dx la zona n_pol+1
dl(3:4,n_pol+1:n_seg)=[ 2:n_pol+1 ; [n_pol+1, 2:n_pol ]];

% archi = a sin la zona i+1, a destra niente
dl(3:4,n_seg+1:nbs)=[ 2:n_pol+1 ; zeros(1,n_pol)];

% ritorna dl
if nargin==1   
  x=dl(:,bs);   
  return 
end 

% calcolo di x(s),y(s)
[m,n]=size(s);
x=zeros(m,n);
y=zeros(m,n);
[mm,nn]=size(bs); 
if mm==1 && nn==1,   
  bs=bs*ones(m,n); % expand bs 
elseif mm~=size(s,1) || nn~=size(s,2),   
  error('bs must be scalar or of same size as s'); 
end 

for ii=1:m
    for jj=1:n
        
        if bs(ii,jj)<=n_pol
            % lati del poligono interno
            ind_i=bs(ii,jj);
            ind_f=mod(bs(ii,jj),n_pol)+1;
            theta_i=theta_inc(ind_i);
            theta_f=theta_inc(ind_f);
            r_i=(R-l_inc(ind_i))*[cos(theta_i), sin(theta_i)];
            r_f=(R-l_inc(ind_f))*[cos(theta_f), sin(theta_f)];
            x(ii,jj)=(r_f(1)-r_i(1))*s(ii,jj)+r_i(1);
            y(ii,jj)=(r_f(2)-r_i(2))*s(ii,jj)+r_i(2);
            
        elseif bs(ii,jj)<=n_seg
            % incisure
            ind=bs(ii,jj)-n_pol;
            theta=theta_inc(ind);
            r_i=(R-l_inc(ind))*[cos(theta),sin(theta)];
            r_f=R             *[cos(theta),sin(theta)];
            x(ii,jj)=(r_f(1)-r_i(1))*s(ii,jj)+r_i(1);
            y(ii,jj)=(r_f(2)-r_i(2))*s(ii,jj)+r_i(2);
        
        else
            % archi esterni
            x(ii,jj)=R*cos(s(ii,jj));
            y(ii,jj)=R*sin(s(ii,jj));

        end
        
    end
end

return
