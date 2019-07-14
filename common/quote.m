% quote delle sezioni della mesh del volume
function [Z_s, sd2sez]=quote(H, n_sez, flag_geom_sp, dz_0, n_sd, Z_sd)

% output

% Z_s quote delle sezioni della mesh del volume
% sd2sez contiene il numero di sezione cui corrisponde ciascun disco speciale


% la spaziatura longitudinale puè essere, in dipendenza da flag_geom_sp:
% lineare, con spaziature costanti
% logaritmica, con progressione geometrica delle spaziature allontanandosi da ciascun disco speciale, con primo interspazio dz_0

% il numero di sezioni nella discretizzazione:
% nel caso di spaziatura costante (flag_geom_sp=false), deve essere maggiore o uguale di 1+2*n_sd,
% nel caso di spaziatura logaritmica (flag_geom_sp=true) deve essere maggiore o uguale di 1+4*n_sd,

% dz_0 è il passo della discretizzazione longitudinale in adiacenza ad un disco speciale (usato solo se flag_geom_sp==true)
% se troppo piccolo, le sezioni si affollano intorno ai dischi speciali

% numero di segmenti in cui è partizionato [0,H]
n_seg=n_sez-1;

% zona di influenza del disco speciale i: 
% da punto intermedio fra due ds (o base inf) a punto intermedio fra due ds (o base sup)
% cioè da Z_sd(i)-lm(i) a Z_sd(i)+lp(i)
lm=[Z_sd(1), diff(Z_sd)/2];
lp=[diff(Z_sd)/2, H-Z_sd(end)];

% a partire da ciascun disco speciale, verso sotto e verso sopra nella sua area di influenza, 
% vengono disposti segmenti tutti uguali (progressione aritmentica) o di lunghezze in progressione geometrica

% incognite
% r ragione della progressione, aritmetica o geometrica
% nm, np: numero di segmenti sopra e sotto ciascun disco speciale

% equazioni:
% il numero di pezzi è n_seg
% la somma della progressione sotto e sopra ciascun disco speciale è rispettivamente lm o lp

% valore di primo tentativo: ragione e numero di segmenti in proporzione alle lunghezze
if flag_geom_sp
    r0=1.1;
else
    r0=H/n_seg;
end
x0=[r0, lm/H*n_seg, lp/H*n_seg];

% risolve
x = fsolve(@fun, x0, optimset('fsolve'), n_sd, n_seg, lm, lp, dz_0, flag_geom_sp); 

% approssima nm e np ad un valore intero
nm=round(x(2:n_sd+1));
np=round(x(n_sd+2:2*n_sd+1));

% si assicura che il numero complessivo di pezzi sia n_seg
% eventuale differenza
dp=n_seg-(sum(nm)+sum(np));
% attribuisce la differenza alla zona più lunga
if max(lm)>=max(lp)
    % la zona più lunga è quella inferiore
    [X,I]=max(lm);
    nm(I(1))=nm(I(1))+dp;
else
    % la zona più lunga è quella superiore
    [X,I]=max(lp);
    np(I(1))=np(I(1))+dp;
end
% controllo sul minimo numero di pezzi
if ((~flag_geom_sp) && (min([nm,np])<1)) || ((flag_geom_sp) && (min([nm,np])<2))
    disp([nm,np]);
    error('errore in quote: aumentare il numero di sezioni')
end

% calcola la ragione effettiva per ogni zona e crea le z
% infatti la ragione potrebbe essere variata a causa dell'arrotondamento sul numero di elementi nella zona
Z_s=[];
r0=x(1);
for d=1:n_sd
    % ragione effettiva nella zona di influenza del disco speciale d, sotto e sopra di questo
    rm = fzero(@fun_loc, r0, optimset('fsolve'), nm(d), lm(d), dz_0, flag_geom_sp);
    rp = fzero(@fun_loc, r0, optimset('fsolve'), np(d), lp(d), dz_0, flag_geom_sp);
    % lunghezze degli elementi
    if flag_geom_sp
        lenm=logspace(log10(dz_0),log10(dz_0*rm^(nm(d)-1)),nm(d));
        lenp=logspace(log10(dz_0),log10(dz_0*rp^(np(d)-1)),np(d));
    else
        lenm=rm*ones(1,nm(d));
        lenp=rp*ones(1,np(d));
    end
	% coordinate dei nodi
	v=cumsum([fliplr(lenm), lenp]);
    % porta la sezione del disco speciale alle quota giusta
    v=v-v(nm(d))+Z_sd(d);
    % aggiunge
    Z_s=[Z_s, v];
end
% prima sezione a quota nulla
Z_s=[0, Z_s];

% calcola sd2sez
% numero cumulativo di segmenti che precedono dischi speciali o sezioni a metà fra dischi speciali
v=cumsum(reshape([nm; np],1,2*n_sd));
sd2sez=1+v(1:2:2*n_sd);

% % controllo che le quote dei ds tornino giuste
% Z_s(sd2sez),Z_sd

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun(x, n_sd, n_seg, lm, lp, dz_0, flag_geom_sp)
% incognite, nell'ordine:
% r: è la lunghezza costante di tutti i segmenti, se flag_geom_sp=false;
%    è la ragione della progressione geometrica delle lunghezze (con primo
%    termine dz_0), se flag_geom_sp=true
% nm: numero di segmenti sotto il disco speciale
% np: numero di segmenti sopra il disco speciale

% componenti di x
r=x(1);
nm=x(2:n_sd+1);
np=x(n_sd+2:2*n_sd+1);

% output
y=zeros(1,2*n_sd+1);

% prima equazione:
% il numero totale di segmenti è n_seg
y(1)=sum(nm)+sum(np)-n_seg;

% ulteriori 2*n_sd equazioni:
% la somma della progressione delle lunghezze dei segmenti sotto è lm(i), poi
% la somma della progressione delle lunghezze dei segmenti sopra è lp(i)
if flag_geom_sp
	y(2:n_sd+1)       =dz_0*(r.^nm-1)/(r-1)-lm;
	y(n_sd+2:2*n_sd+1)=dz_0*(r.^np-1)/(r-1)-lp;
else
	y(2:n_sd+1)       =nm*r-lm;
	y(n_sd+2:2*n_sd+1)=np*r-lp;
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=fun_loc(r, n, ltot, dz_0, flag_geom_sp)
% incognita r:
% è la lunghezza costante di tutti i segmenti, se flag_geom_sp=false;
% è la ragione della progressione geometrica delle lunghezze (con primo
% termine dz_0), se flag_geom_sp=true

% equazione: somma della progressione delle lunghezze pari a ltot
if flag_geom_sp
    y=dz_0*(r^n-1)/(r-1)-ltot;
else
    y=n*r-ltot;
end

return
