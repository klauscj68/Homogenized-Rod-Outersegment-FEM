% aggiorna le posizioni delle particelle
function p_new=random_update(p, delta, R, n_inc, l_inc, theta_inc, n_theta_inc, ...
    plot_pool, c, mark)

% input:
% p         posizioni all'inizio dello step
% delta     vettore spostamento da percorrere
% R         raggio della rirconferenza
% n_inc, l_inc, theta_inc, n_theta_inc  dati sulle incisure: 
%           n_theta_inc è la normale alla linea di rimbalzo nel punto di rimbalzo
% plot_pool flag per il plottaggio del biliardino
% c         colori dei vari path
% mark      carattere scritto nel punto di rimbalzo

% aggiorna le posizioni, come primo tentativo, assumendo che non ci sia alcun rimbalzo
p_new=p+delta;

% pesca i punti che, partendo da p con spostamento delta, 
% tentano di evadere dal cerchio di raggio R o intersecano le incisure
% restituisce la loro lista (esc), i punti di impatto (p_touch) 
% ed il vettore di spostamento che manca  da percorrere, 
% specchiato rispetto alla normale al punto di impatto (delta_mir)
% p_touch e delta_mir hanno length(esc) colonne
[esc, p_touch, delta_mir]=escape(p, delta, R, n_inc, l_inc, theta_inc, n_theta_inc);

if plot_pool
    % disegna il random path
    plot_path_main(p, p_new, p_touch, esc, c, mark);
end

iter=0;
while ~isempty(esc)
    iter=iter+1;

    % aggiorna le posizioni dei punti che al passo precedente hanno
    % sbattuto, ipotizzandom quale primo tentativo, che non sbatteranno ancora
    p_new(:,esc)=p_touch+delta_mir;

    % pesca i punti che, partendo da p_touch con spostamento delta_mir, 
    % tentano di evadere dal cerchio di raggio R o intersecano le incisure
    % restituisce la loro lista (esc_new), i punti di impatto (p_touch_new) 
    % ed il vettore di spostamento che manca da percorrere, 
    % specchiato rispetto alla normale al punto di impatto (delta_mir_new)
    % nota: esc_new è nella numerazione relativa a p_touch (1..length(esc))
    % p_touch_new e delta_mir_new hanno length(esc_new) colonne
    [esc_new, p_touch_new, delta_mir_new]=escape(p_touch, delta_mir, R, n_inc, l_inc, theta_inc, n_theta_inc);

    if plot_pool
        % disegna il random path
        plot_path_main(p_touch, p_new(:,esc), p_touch_new, esc_new, c(esc), mark);
    end

    % prepara per il ciclo successivo
    p_touch=p_touch_new;
    delta_mir=delta_mir_new;
    % gli indici esc_new, come ottenuti dal find precedente, si riferiscono
    % alla numerazione locale di esc (cioè relativa all'ordine in p_touch)
    % sono qui convertiti alla numerazione globale 1..N
    esc=esc(esc_new);
    
% pause
end % while
% iter
return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_path(p, p_new, c)
    % particelle
    for i=1:size(p,2)
        plot([p(1,i),p_new(1,i)],[p(2,i),p_new(2,i)],'Linewidth',1,'Color',c(i))
    end
    drawnow;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_path_main(p, p_new, p_touch, esc, c, mark)
    % disegna il random path
    % p         punti di origine
    % p_new     posizioni finali dei punti
    % p_touch   punti di impatto, nella numerazione  1..length(esc)

    % punti che sbattono: esc
    % punti che non sbattono
    ok=setdiff(1:size(p,2),esc);
    
    % punti che non sbattono: si arriva fino alla fine
    plot_path(p(:,ok),p_new(:,ok),c(ok));

    % punti che sbattono: si arriva fino al p_touch
    plot_path(p(:,esc),p_touch,c(esc));
    
    % mark nel punto di impatto
    plot(p_touch(1,:),p_touch(2,:),['k',mark]);

    % se qualcuno sbatte, beep
    if ~isempty(esc)
        beep
    end
return
