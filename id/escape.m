% pesca i punti che tentano di evadere dal cerchio di raggio R
% ed i punti che intersecano qualche incisura
% per essi, individua il punto di intersezione 
% e specchia il cammino residuo rispetto alla normale
function [esc, p_touch, delta_mir]=escape(p, delta, R, n_inc, l_inc, theta_inc, n_theta_inc)

% esc contiene la lista dei punti che nel segmento da p a p+delta
% intersecano un'incisura o la circonferenza

% p_touch contiene le coordinate del punto di intersezione

% delta_mir è il vettore di spostamento che residua da percorrere, a
% partire da p_touch, specchiato rispetto alla normale

% nota: se il segmento da p a p+delta tocca più di una linea, tutte le
% quantità di sopra si riferiscono alla prima intersezione


% inizializza le matrici di lavoro 
% qui sono affastellati, sulle successive colonne, 
% prima i punti che evadono dalla circonferenza, 
% poi quelli che toccano la prima incisura, poi la seconda ecc.
% naturalmente, un punto che interseca linee multiple appare più volte
esc_dup=zeros(1,0);
t_touch_dup=zeros(1,0);
p_touch_dup=zeros(2,0);
delta_mir_dup=zeros(2,0);


% Tolleranza sulle eguaglianze fra numeri reali
tol=1000*eps;

% posizione finale, se non ci fossero riflessioni
p_new=p+delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% punti che fuoriescono dalla circonferenza
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evasi
esc_c=find(dot(p_new,p_new,1)>R^2-tol);

% calcolo di t_touch_c e p_touch_c
% i punti si spostano con la legge p(t)=p+t*delta, con 0\le t \le 1
% per i punti che tentano di evadere, trova il parametro \bar t a cui p(t) 
% è sulla circonferenza
% risolve l'equazione ||p(t)||^2=R^2
% cioè delta^2 t^2 + 2*delta*p t + ||p||^2-R^2=0
% punti interessati
p_esc_c=p(:,esc_c);
% rispettivi spostamenti
delta_esc_c=delta(:,esc_c);
% calcolo del t_esc_c
psq_esc_c=dot(p_esc_c,p_esc_c,1);
deltap_esc_c=dot(delta_esc_c,p_esc_c,1);
deltasq_esc_c=dot(delta_esc_c,delta_esc_c,1);
t_touch_c=(-deltap_esc_c+sqrt(deltap_esc_c.^2-deltasq_esc_c.*(psq_esc_c-R^2)))./deltasq_esc_c;

% punti in cui gli evasi toccano la circonferenza
p_touch_c=p_esc_c+repmat(t_touch_c,2,1).*delta_esc_c;

% vettore che resta da percorrere a partire da p_touch_c
delta_res_c=delta_esc_c.*(1-repmat(t_touch_c,2,1));

% versori normali esterni alla circonferenza nei punti p_touch_c
n_touch_c=p_touch_c./repmat(sqrt(dot(p_touch_c,p_touch_c,1)),2,1);

% vettore residuo di spostamento, specchiato rispetto al versore tangente
% l'operatore di specchiatura è I-2 n\otimes n
delta_mir_c=delta_res_c-2*repmat(dot(delta_res_c,n_touch_c,1),2,1).*n_touch_c;

% assembla nei vettori complessivi
esc_dup      =[esc_dup,       esc_c];
t_touch_dup  =[t_touch_dup,   t_touch_c];
p_touch_dup  =[p_touch_dup,   p_touch_c];
delta_mir_dup=[delta_mir_dup, delta_mir_c];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% punti che intersecano un'incisura
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % WARNING
% % % % % % % nel codice seguente si assume che j sia \sqrt(-1)
% % % % % % % quindi j NON deve essere stato utilizzato prima

% ciclo sulle incisure
for i=1:n_inc
    % il segmento da p a p_new=p+delta è p(t) = (1-t)p + t p_new
    % il versore normale all'incisura è n
    % c'è intersezione con la retta sovrapposta all'incisura se 
    % l'equazione n\dot p(t)=0 ha soluzione per t\in(0,1)
    % cioè se n\dot p e n\dot p_new sono discordi
    F0=dot(repmat(n_theta_inc(:,i),1,size(p,2)),p,1);
    F1=dot(repmat(n_theta_inc(:,i),1,size(p,2)),p_new,1);
    % qui la diseguaglianza va con <-tol, perché la ondizione F0=0 
    % (e ciò capita a partire da p_touch dopo una specchiatura)
    % non deve essere riconosciuta passibile di specchiatura, altrimenti il
    % vettore delta_mir comincia a oscillare
    % PERO' questo controllo consente a punti che hanno attraversato
    % un'incisura, rimanendo vicinissimi alla stessa, di sfuggire alla riflessione 
    esc_inc=find(F0.*F1<-tol);
    
    if  ~isempty(esc_inc)
        % gli elementi di esc_inc attraversano la retta di anomalia theta_inc
        % resta da vedere se tagliano l'incisura, che si estende su
        % \rho\in(R-l_inc,R)
        % calcolo il valore di t\in(0,1) per cui n\dot p(t)=0
        t_touch_inc=F0(esc_inc)./(F0(esc_inc)-F1(esc_inc));
        % punti di toccaggio (sulla linea che contiene l'incisura)
        p_touch_inc=p(:,esc_inc).*repmat(1-t_touch_inc,2,1) ...
            +p_new(:,esc_inc).*repmat(t_touch_inc,2,1);
        % vettore che resta da percorrere a partire da p_touch_inc
        delta_res_inc=delta(:,esc_inc).*(1-repmat(t_touch_inc,2,1));
        % vettore normale all'incisurea, replicato per tutti i punti
        n_touch_inc=repmat(n_theta_inc(:,i),1,length(esc_inc));
        % vettore residuo di spostamento, specchiato rispetto al versore tangente
        % l'operatore di specchiatura è I-2 n\otimes n
        delta_mir_inc=delta_res_inc-2*repmat(dot(delta_res_inc,n_touch_inc,1),2,1).*n_touch_inc;
        % raggio vettore e anomalia dei punti di toccaggio
        rho_p_touch_inc=abs(p_touch_inc(1,:)+j*p_touch_inc(2,:));
        angle_p_touch_inc=angle(p_touch_inc(1,:)+j*p_touch_inc(2,:));
        % effettivo toccaggio se \rho\in(R-l_inc,R) e se l'anomalia differisce
        % di multipli di 2*pi da theta_inc
        esc_inc_eff=find( (rho_p_touch_inc>=R-l_inc(i)) & (abs(modulo_near(angle_p_touch_inc-theta_inc(i),2*pi))<tol) );
        
        % assembla nei vettori complessivi
        esc_dup      =[esc_dup,       esc_inc(esc_inc_eff)];
        t_touch_dup  =[t_touch_dup,   t_touch_inc(esc_inc_eff)];
        p_touch_dup  =[p_touch_dup,   p_touch_inc(:,esc_inc_eff)];
        delta_mir_dup=[delta_mir_dup, delta_mir_inc(:,esc_inc_eff)];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sintesi (superiore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per ogni nodo che compare più di una volta, prende i dati relativi al
% t_touch più piccolo

if isempty(esc_dup)
    n_gruppi_sort=0;
    % inizializza matrici di uscita vuote
    esc=zeros(1,n_gruppi_sort);
    p_touch=zeros(2,n_gruppi_sort);
    delta_mir=zeros(2,n_gruppi_sort);
else
    % ordina esc, per avere vicine le copie
    % I contiene le posizioni dentro esc_dup
    [esc_sort,I]=sort(esc_dup);
    % essendo esc_sort=esc_dup(I), permuta analogamente anche 
    % t_touch_dup, p_touch_dup, delta_mir_dup
    t_touch_sort=t_touch_dup(I);
    p_touch_sort=p_touch_dup(:,I);
    delta_mir_sort=delta_mir_dup(:,I);

    % scorre esc_sort, raggruppandolo idealmente per gruppi di copie multiple
    % pesca gli indici entro esc_sort in cui sono presenti variazioni nella
    % lista, dovute al passaggio da un gruppo ad un altro
    variaz_sort=[0,find(diff(esc_sort)),length(esc_sort)];
    % pesca le lunghezze dei singoli gruppi, come diff degli indici in cui c'è
    % una variazione
    gruppi_sort=diff(variaz_sort);
    % numero di gruppi
    n_gruppi_sort=length(gruppi_sort);
    % inizializza matrici di uscita
    esc=zeros(1,n_gruppi_sort);
    p_touch=zeros(2,n_gruppi_sort);
    delta_mir=zeros(2,n_gruppi_sort);
    % ciclo su tutti i gruppi
    for g=1:n_gruppi_sort
        % gruppo g, che occupa gli indici da variaz_sort(g)+1 a variaz_sort(g+1) in esc_sort
        % nell'ambito di questo gruppo di copie multiple, sceglie il minimo t_touch
        [temp,pos]=min(t_touch_sort(variaz_sort(g)+1:variaz_sort(g+1)));
        % nuovo punto, in posizione variaz_sort(g)+pos nella lista sort complessiva
        esc(g)=esc_sort(variaz_sort(g)+pos);
        p_touch(:,g)=p_touch_sort(:,variaz_sort(g)+pos);
        delta_mir(:,g)=delta_mir_sort(:,variaz_sort(g)+pos);
    end
end
% amen

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% routine per calcolare il modulo con approssimazione all'intero più vicino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=modulo_near(x,y)
    f=x-round(x/y)*y;        
return
