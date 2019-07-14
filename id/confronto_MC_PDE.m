% confronto_MC_PDE
function confronto_MC_PDE(n_ipd_id, p_ipd_id, n_tri_id, t_ipd_id, ...
     n_sd, n_inc, inc, ...
     cc_R_st, D_R_st, ...
     n_Phi, Phi, ...
     tau, mode_time, ...
     peak_delta, R, ...
     inspect, plot_pool, plot_niagara, ...
     t_fin, n_step_t, metodo_id, solver)

% un solo disco speciale
if n_sd~=1
    error('in confronto_MC_PDE è ammesso un solo disco speciale');
end
d=1;
% una sola locazione
if n_Phi~=1
    error('in confronto_MC_PDE è ammessa una sola locazione');
end
j=1;


% coordinate polari della fotoisomerizzazione
rho=Phi{d}(1,j);
theta=Phi{d}(2,j);
% coordinate cartesiane della fotoisomerizzazione
Xp=rho*cos(theta);
Yp=rho*sin(theta);

% il confronto fra la soluzione PDE-FEM e quella MC
% è fatto importando entrambe le soluzioni sui vertici di una grid

% passo del grid
dx_grid=R/75;
% coordinate del grid
X= -R+dx_grid/2:dx_grid:R-dx_grid/2;
% numero di maglie lungo un asse
n_passi_grid=length(X);

% movimento random
% passo temporale
dt=t_fin/n_step_t;

% cammino libero medio (in 2 dimensioni spaziali)
dx=sqrt(2*2*D_R_st*dt);

% numero di particelle
% N=4;
N=100000;
% N=1000000;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% costruzione  delle matrici di massa e rigidezza, e delle collezioni dei vettori dei carichi degli elementi
% nota: le matrici K e M non sono moltiplicate per alcun coefficiente
% contengono gli integrali dei prodotti dei gradienti delle funzioni di forma (K) o dei prodotti delle funzioni di forma (M)
[K, M] = assembla_ipd(n_ipd_id, n_tri_id, p_ipd_id, t_ipd_id);

% chiamata dummy a find_tri, per determinare le matrici pt_min, pt_max
[tri, F, pt_min, pt_max]=find_tri(p_ipd_id, t_ipd_id, [0;0], R);

% vettore dato iniziale
delta_nod=dirac(n_ipd_id, p_ipd_id, t_ipd_id, M, [Xp;Yp], peak_delta, R, pt_min, pt_max);

% istanti di calcolo
time=(0:n_step_t)'*dt;

% assenza di assorbimrnto
R_st_id_diff=diffusione_ipd_single(n_ipd_id, cc_R_st*M, D_R_st*K, ...
    zeros(n_step_t+1,n_ipd_id), delta_nod', ...
    t_fin, n_step_t, metodo_id, solver, theta, inspect, ['R*, disco ',num2str(d)]);

% importa nel grid la soluzione FEM, ad ogni istante temporale
R_st_grid_FEM=zeros(n_step_t+1,n_passi_grid^2);
% istante iniziale
[R_st_grid_matr,tn,a2,a3]=tri2grid(p_ipd_id,t_ipd_id,R_st_id_diff(1,:)',X,X);
R_st_grid_FEM(1,:)=R_st_grid_matr(:);
for passo=2:n_step_t+1
    R_st_grid_matr=tri2grid(p_ipd_id,t_ipd_id,R_st_id_diff(passo,:)',tn,a2,a3);
    R_st_grid_FEM(passo,:)=R_st_grid_matr(:);
    % figura surf
    figure(15)
    h=surf(X,X,R_st_grid_matr);
    axis([-R,R,-R,R,0,1])
    shading interp
    view(3)
    light('Position',[ 0 -5  1000],'Style','infinite');
    colormap cool;
%     colormap copper;
%     colormap hot;
%     shading flat
    drawnow;
%     pause
end

% genera la griglia
[XX,YY] = meshgrid(X, X);
p_grid=[XX(:)';YY(:)'];

save temp_FEM

% % butta via i vertici fuori del cerchio
% % ind_grid_cerchio è una riga
% ind_grid_cerchio=find(p_grid(1,:).^2+p_grid(2,:).^2<=R^2);

% % trova i tri per ogni punto del grid
% % la matrice F_elem, ordinata in colonna, dovra essere moltiplicata
% % punto-punto per il vettore delle particelle in ogni nodo, opportunamente triplicato
% % questa operazione è fatta da random_update
% [tri, F_elem, pt_min, pt_max]=find_tri_large_N(p_ipd, t_ipd, p_grid(:,ind_grid_cerchio), R);
% 
% % elenco dei nodi di ciascun triangolo relativo a ciascun punto
% % serve in random_update per l'assemblaggio del vettore dei carichi
% nodi_tri_ipd_id=t_ipd(:,tri);
% valori fittizi, quando usa l'importazione sul grid, per evitare errore
% nella chiamata
nodi_tri_ipd_id=0;
ind_grid_cerchio=0;
F_elem=0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO MONTE CARLO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ogni giro è uguale al precedente
rand('state',13)
randn('state',13)

% geometria delle incisure
% lunghezza
l_inc=inc(1,:);
% posizione angolare
theta_inc=inc(2,:);
cos_theta_inc=cos(theta_inc);
sin_theta_inc=sin(theta_inc);
% normale in senso antiorario
n_theta_inc=[-sin_theta_inc; cos_theta_inc];

% predispone il biliardo (circnferenza e tracce delle incisure 
% per il plottaggio dei random walk di R* ed E*)
if plot_pool
    % figura del random path
    % serve a random_update
    f13=figure(13);
    newplot(f13)
    title('Random walk')
    xlabel('x [\mu m]')
    ylabel('y [\mu m]')
    axis([-R*1.1 R*1.1 -R*1.1 R*1.1])
%     set(f13,'Position',[-30 154 672 504]);
    view(2)
    hold on
    axis([-R,R,-R,R])
    axis equal
    % circonferenza
    angle_c= 0:2*pi/50:2*pi;
    plot(R*cos(angle_c), R*sin(angle_c),'r-','Linewidth',2)
    % incisure
    for cont=1:n_inc
        plot([(R-l_inc(cont))*cos(theta_inc(cont)), R*cos(theta_inc(cont))], ...
            [(R-l_inc(cont))*sin(theta_inc(cont)), R*sin(theta_inc(cont))],'r-','Linewidth',2)
    end
end

% tutte le particelle concentrate in un punto (delta iniziale)
p=repmat([Xp;Yp],1,N);

% posizioni inziali uniformi
% p=zeros(2,N);
% % butta via le particelle fuori del cerchio
% fine=0;
% ind=[1:N];
% while ~fine
%     % genera i punti ancora da allocare (dentro ind)
%     p(:,ind)=(rand(2,length(ind))-1/2)*2*R;
%     % vede quali punti, dentro ind, sono fuori del cerchio
%     ind_new=find((p(1,ind).^2+p(2,ind).^2)>=R^2);
%     ind=ind(ind_new);
%     fine=isempty(ind);
% end
% rho=R*sqrt(rand(1,N));
% th=rand(1,N)*2*pi;
% p=[rho.*cos(th); rho.*sin(th)];

% istanti di spegnimento delle singole particelle
% vita media della R*
t_mean=1/sum(tau);
% genera gli istanti di spegnimento, aventi media t_mean
if mode_time==1
    % spegnimento esponenziale
    t_shutoff=exprnd(t_mean,N,1);
else
    % sempre accese (non considera mode_time==2)
    t_shutoff=ones(1,N)*t_fin+1000*eps;
end

% % % % per la figura del path
% % % p_old=p;

% soluzione MC, ad ogni istante temporale
R_st_grid_MC=zeros(n_step_t+1,n_passi_grid^2);
% istante iniziale
dens=full(sparse(fix((R+p(1,:)')/dx_grid)+1, fix((R+p(2,:)')/dx_grid)+1, ones(N,1), n_passi_grid, n_passi_grid)/N/dx_grid^2);
R_st_grid_MC(1,:)=dens(:);

%n_step_t=10*n_step_t;
% n_step_t=5;
% particelle sopravviventi a tempo zero
N_survivors=N;

% colori (per plot_pool)
c=[];
while length(c)<N
    c=[c,'bgrcmyk'];
end

tic
for i=1:n_step_t
    % mostra lo stato di avanzamento
    fprintf('Diffusione Montecarlo: passo %4.0f di %4.0f\n',i,n_step_t);
    
%     % genera gli spostamenti random
%     % angolo di spostamento distribuito con legge uniforme in [0,2*pi)
%     th=rand(1,N_survivors)*2*pi;
% 
%     % ampiezza di spostamento costante
%     amp=dx*ones(1,N_survivors);
% 
%     % ampiezza di spostamento distribuita con legge uniforme nel cerchio di
%     % raggio 2*dx, in modo che il valore quadratico medio sia dx
%     % le ampiezze vanno distribuite con densità crescente con \rho
%     amp=2*dx*sqrt(rand(1,N_survivors));
% 
%     % ampiezza di spostamento distribuita con legge gaussiana, di media mu e vaianza sigma
%     % in modo che il valore quadratico medio \mu^2+\sigma^2 sia pari a dx^2
%     sigma=0.1;
%     if sigma>=dx
%         error('la sigma non può superare dx')
%     end
%     mu=sqrt(dx^2-sigma^2);
%     amp=sigma*randn(1,N_survivors)+mu;
% 
%     % vettore di spostamento
%     delta=repmat(amp,2,1).*[cos(th); sin(th)];

    % ampiezza di spostamento distribuita con legge gaussiana
    % bidimensionale di media nulla e varianza dx^2
    % la gaussiana 2D è ottenuta come prodotto cartesiano di gaussiane 1D
    % (X e Y) indipendenti ed identicamente distribuite, ciascuna con varianza sigma^2
    % quindi la varianza della gaussiana 2D è 2\sigma^2=dx^2
    sigma=dx/sqrt(2);
    delta=sigma*[randn(1,N_survivors); randn(1,N_survivors)];
    
    if plot_pool
        figure(13)
    end

    % aggiorna le posizioni delle particelle
    p=random_update(p, delta, R, n_inc, l_inc, theta_inc, n_theta_inc, ...
        plot_pool, c, '*');    

    % uccide le particelle che hanno t_shutoff<i*n_step_t
    survivors=find(t_shutoff>i*dt);
    N_survivors=length(survivors);
    p=p(:,survivors);
    t_shutoff=t_shutoff(survivors);

    % importa nel grid
    % ogni paticella vale 1: da ciò l'ones alla fine
    % gli indici si ottengono col fix
    dens=full(sparse(fix((R+p(1,:)')/dx_grid)+1, fix((R+p(2,:)')/dx_grid)+1, ones(N_survivors,1), n_passi_grid, n_passi_grid)/N/dx_grid^2);
    R_st_grid_MC(i+1,:)=dens(:);

    if plot_niagara
        % figura della cascata
        plot_niagara_main(p, R, dx_grid, n_passi_grid, nodi_tri_ipd_id, ind_grid_cerchio, F_elem, M, n_ipd_id, p_ipd_id, t_ipd_id);
    end

    presskey(inspect);
end
toc

save temp_MC

% confronto MC-FEM
% norma L^\inf(L^1)
norma1=zeros(n_step_t+1,1);
massa_FEM=zeros(n_step_t+1,1);
massa_MC=zeros(n_step_t+1,1);
giratore_FEM=zeros(n_step_t+1,1);
giratore_MC=zeros(n_step_t+1,1);
inerzia_FEM=zeros(n_step_t+1,1);
inerzia_MC=zeros(n_step_t+1,1);
for i=1:n_step_t+1
    ind=find(~isnan(R_st_grid_FEM(i,:)));
    norma1(i)=norm(R_st_grid_MC(i,ind)-R_st_grid_FEM(i,ind))/norm(R_st_grid_FEM(i,ind),1);
    massa_FEM(i)=sum(R_st_grid_FEM(i,ind))*dx_grid^2;
    inerzia_FEM(i)=sum(R_st_grid_FEM(i,ind).*((p_grid(1,ind)-Xp).^2+(p_grid(2,ind)-Yp).^2))*dx_grid^2;
    giratore_FEM(i)=sqrt(inerzia_FEM(i)/massa_FEM(i));
    massa_MC(i)=sum(R_st_grid_MC(i,ind))*dx_grid^2;
    inerzia_MC(i)=sum(R_st_grid_MC(i,ind).*((p_grid(1,ind)-Xp).^2+(p_grid(2,ind)-Yp).^2))*dx_grid^2;
    giratore_MC(i)=sqrt(inerzia_MC(i)/massa_MC(i));
end
figure(15)
plot(time,norma1,'r-','Linewidth',2)
axis([0 t_fin 0 2])

figure(16)
hold on
plot(time,giratore_FEM,'r-','Linewidth',2)
plot(time,giratore_MC,'g-','Linewidth',2)
axis([0 t_fin 0 5])

figure(17)
hold on
plot(time,massa_FEM,'r-','Linewidth',2)
plot(time,massa_MC,'g-','Linewidth',2)
axis([0 t_fin 0 1])

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_niagara_main(p, R, dx_grid, n_passi_grid, nodi_tri_ipd_id, ind_grid_cerchio, F_elem, M, n_ipd, p_ipd, t_ipd)

%     % numero di particelle
%     N=size(p,2);
%     % ogni paticella vale 1: da ciò l'ones alla fine
%     % gli indici si ottengono col fix
%     num_grid=full(sparse(fix((R+p(1,:)')/dx_grid)+1, fix((R+p(2,:)')/dx_grid)+1, ones(N,1), n_passi_grid, n_passi_grid));
%     % seleziona solo i punti della grid che sono entro il cerchio
%     % qui il num_grid è pensato in colonna num_grid(:)
%     num_grid=num_grid(ind_grid_cerchio);
% 
%     % triplica ogni elemento di num_grid
%     % num_gridè una riga, perché ind_grid_cerchio è una riga
%     num_grid_col=repmat(num_grid,3,1);
%     num_grid_col=num_grid_col(:);
%     % calcola il vettore dei termini noti
%     % nodi_tri_ipd_id(:) contiene, per ciascun punto della grid dentro il cerchio, 
%     % i tre nodi del triangolo che contiene tale punto
%     % i punti della grid nel cerchio sono ind_grid_cerchio
%     % F_elem(:) contiene, per ciascun punto della grid dentro il cerchio, 
%     % i valori delle funzioni di forma in esso, relative al triangolo che contiene tale punto
%     % num_grid_col contiene il numero di particelle attribuite a ciascun
%     % punto della grid dentro il cerchio, replicato a tre a tre
%     V=sparse(nodi_tri_ipd_id(:),ones(3*length(ind_grid_cerchio),1),F_elem(:).*num_grid_col,n_ipd,1);
%     % convert from load vector to nodal values for densities
%     dens=V/ones(1,n_ipd)*M*V;
%     dens=M\V;

%     % importa le particelle nella mesh di triangoli
%     % l'idea è di attribuire a ciascun nodo della mesh di triangoli 
%     % una densità media in un suo intorno circolare di raggio r_intorno
%     % ottenuta come numero di particelle presenti in questo intorno diviso
%     % per la sua area
%     % nel caso di nodi che distano meno di r_intorno dal bordo del cerchio
%     % di raggio R, l'area dell'intorno è diminuita dell'area della zona 
%     % dell'intorno di raggio r_intorno che cade fuori del cerchio di raggio R
%     r_intorno=R/20;
%     dens=zeros(n_ipd,1);
%     area_infl=pi*r_intorno^2*ones(n_ipd,1);
%     % punti della mesh vicini al bordo
%     ind_bordo=find(p_ipd(1,:).^2+p_ipd(2,:).^2>(R-r_intorno)^2);
%     h=R-sqrt(p_ipd(1,ind_bordo).^2+p_ipd(2,ind_bordo).^2);
%     theta=asin(h/r_intorno);
%     area_infl(ind_bordo)=pi*r_intorno^2/2+r_intorno*h.*cos(theta)+r_intorno^2*theta;
%     % ciclo su tutti i nodi della mesh
%     for n=1:n_ipd
%         % trova i punti in p a distanza da p_ipd(:,n) minore di r_intorno
%         dens(n)=length(find((p_ipd(1,n)-p(1,:)).^2+(p_ipd(2,n)-p(2,:)).^2<=r_intorno))/area_infl(n);
%     end
% 
%     % figura
%     f1=figure(14);
%     newplot(f1)
%     title('density')
%     % chromatic scale RGB
%     colori=[1-dens/max(dens), dens/max(dens), zeros(n_ipd,1)];
%     % Mesh of triangles
%     vert=[p_ipd', dens];
%     patch('Vertices',vert,'Faces',t_ipd','FaceVertexCData',colori,...
%           'FaceColor','interp','EdgeColor','interp');
%     xlabel('x [\mu m]')
%     ylabel('y [\mu m]')
%     zlabel('dens [mol/\mu m^2]')
%     axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 5e5])
%     view(3)
%     drawnow;
% %     pause


    % importa nel grid
    X= -R+dx_grid/2:dx_grid:R-dx_grid/2;
    n_passi_grid=length(X);
    N=size(p,2);
    % ogni paticella vale 1: da ciò l'ones alla fine
    % gli indici si ottengono col fix
    dens=full(sparse(fix((R+p(1,:)')/dx_grid)+1, fix((R+p(2,:)')/dx_grid)+1, ones(N,1), n_passi_grid, n_passi_grid)/N/dx_grid^2)';

    % figura surf
    figure(14)
    h=surf(X,X,dens);
    axis([-R,R,-R,R,0,1])
    shading interp
    view(3)
    light('Position',[ 0 -5  1000],'Style','infinite');
    colormap cool;
%     colormap copper;
%     colormap hot;
%     shading flat
    drawnow;
return
