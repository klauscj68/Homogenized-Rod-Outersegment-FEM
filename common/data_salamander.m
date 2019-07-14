% Dati
% Data
function [R, H, n_sez, flag_geom_sp, dz_0, ...
    n_sd, Z_sd, ...
    n_inc, inc, ...
    taglia, tol_R, tol_angle, ...
    n_ref_cyto, n_ref_id, metodo_cyto, metodo_id, ...
    solver, theta, alpha, tol_fix, norma_inf, peak_delta, ...
    t_fin, n_step_t, downsample, ...
    plot_mesh, plot_num, plot_niagara, plot_pool , inspect, ...
    u_tent, v_tent, tol_stat, ...
    flag_model, flag_model_disc, flag_Ca_clamp, ...
    mode_time, mode_space, n_step_R, mu, lambda, ...
    nu_RE, nu_RG, cc_R_st, D_R_st, ...
    n_Phi, Phi, random_location, ...
    k_GE, cc_G_st, D_G_st, ...
    k_hyd, PDE_s, k_st, ...
    cc_E_st, D_E_st, k_E, ...
    MC_disk, ...
    n_sample, ...
    epsilon_0, nu, sigma, ...
    cc_u, kk_u, cc_v, kk_v, ...
    alpha_max, alpha_min, m_cyc, k_cyc, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg, ...
    j_ex_sat, K_ex]=data



% PARAMETRI GEOMETRICI
% GEOMETRICAL PARAMETERS

% raggio del bastoncello
% rod radius
R=5.5; % [\mu m]
% altezza del bastoncello
% rod height
H=22.4; % [\mu m]

% numero di sezioni trasversali nella discretizzazione
% nel caso di spaziatura costante (flag_geom_sp=false), deve essere maggiore o uguale di 1+2*n_sd, (n_sd è il numero dei dischi speciali)
% nel caso di spaziatura logaritmica (flag_geom_sp=true) deve essere maggiore o uguale di 1+4*n_sd,
% number of cross sections to be used in the discretization
% if a linear spacing is used (flag_geom_sp=false), n_sez must be greater than or equal to 1+2*n_sd (n_sd in the number of special discs)
% if a log spacing is used (flag_geom_sp=true), n_sez must be greater than or equal to 1+4*n_sd
% n_sez=7;
% n_sez=25;
% n_sez=30;
 n_sez=31;

% tipo di spaziatura longitudinale 
% true=progressione geometrica degli spazi fra sezioni trasversali consecutive, a partire dai dischi speciali;
% false=spaziature costanti fra sezioni trasversali consecutive
% type of spacing between cross sections
% true=geometric progression of interspaces between cross sections, starting a each special disc
% false=constant spacing between cross sections
flag_geom_sp=true;

% passo della discretizzazione longitudinale in adiacenza ad un disco speciale (usato solo se flag_geom_sp==true)
% può essere preso pari a (1+nu)*epsilon_0, ma così le sezioni si affollano intorno ai dischi speciali
% step of longitudinal space discretization near a special disc
% may be taken equal to (1+nu)*epsilon_0, but in this way discretization becomes crowded near special discs
dz_0=29e-3; % []


% DATI SULL'ATTIVAZIONE
% DATA ON THE ACTIVATION PROCESS

% numero dei dischi speciali
% number of special discs
n_sd=1;
% n_sd=2;

% quote dei dischi speciali: la lunghezza del vettore Z_sd deve essere n_sd
% position of the special discs wrt the basis: the length of Z_sd must be n_sd
Z_sd=[H/2]; % [\mu m]
% Z_sd=[H/2-H/400, H/2+H/400]; % [\mu m]


% DATI SULLE INCISURE
% DATA ON INCISURES

% Numero di incisure
% number of incisures
n_inc=0;
inc=zeros(3,0);
% n_inc=1;
% n_inc=4;

% I dati sulle incisure sono allocati nella matrice inc, con n_inc colonne
% ogni colonna ha tre componenti, che sono la lunghezza dell'incisura; la posizione angolare; l'apertura
% la lunghezza di ciascuna incisura deve essre minore del raggio
% le anomalie a cui sono posizionate le incisure devono essere tutte positive 
% e crescenti in senso antiorario al crescere dell'indice di incisura (anche oltre 2*pi se necessario)
% se sono presenti solo due incisure, l'angolo dalla prima alla seconda deve essere minore del piatto
% l'angolo di apertura (beanza) delle incisure, assunte essere settori circolari, è misurato rispetto al verice dell'incisura
% the data on incisures are in the matrix inc, with n_inc colums
% each column has three compoments, referring to the current incisure: length; angular position; angular amplitude
% the length of each incisure must be less than the radius
% the angles must be all positive and increasing in counterclockwise direction (may be greater that 2*pi if necessary)
% if there are only two incisures, the angle from the first one to the second one must be less than pi
% the angular amplitude of each incisure, regarded as a circular sector, is measured with respect to its vertex
% inc=zeros(3,0);
% inc=[ [4.5; 0; pi/50] ]; % [ \mu m; rad; rad] on each column
% inc=[ [3; 0*pi; pi/500],[3; 1/5*pi; pi/500]]; % [ \mu m; rad; rad] on each column
% inc=[ [3.5; 40/180*pi; pi/50], [4; 80/180*pi; pi/50], [5; 170/180*pi; pi/50] , [3; 290/180*pi; pi/50] ]; % [ \mu m; rad; rad] on each column


% % %% GENERATES EQUALLY SPACED INCISURES (comment if you want to put 0 incisures)
% n_inc=23;
% % width of a single incisure
% l_b=15e-3;
% % total area of incisures
% A_tot=0.80;
% % length of a single incisure
% l_r=2*A_tot/(n_inc*l_b);
% % angle of a single incisure
% theta_inc=l_b/l_r;
% % generate the variable inc
% inc=[l_r*ones(1,n_inc); 0:2*pi/n_inc:2*pi*(1-1/n_inc); theta_inc*ones(1,n_inc)];


% PARAMETRI DI DISCRETIZZAZIONE
% DISCRETIZATION PARAMETERS

% taglia degli elementi prodotti da initmesh (prima di invocare refinemesh)
% size of the elements generated by initmesh (before invoking refinemesh)
taglia=R/4; % [\mu m]

% tolleranza lineare dei setacci nella ricerca dei nodi
% due punti sono considerati coincidenti quando sono a distanza minore di tol_R
% linear tolerance of the grid for searching nodes
% two nodes are assumed to be the same if their distance is less than tol_R
tol_R=taglia/1000; % [\mu m]

% tolleranza angolare dei setacci nella ricerca dei nodi
% angular tolerance of the grid for searching nodes
tol_angle=2*pi/1000; % [rad]

% numero di raffittimenti per la mesh del problema omogeneizzato (numero di volte in cui è chiamato refinemesh)
% number of refinements of the mesh for the homogenized problem (numer of times refinemesh is invoked)
n_ref_cyto=0;

% numero di raffittimenti per la mesh del di diffusione sul disco inciso (numero di volte in cui è chiamato refinemesh)
% deve essere maggiore o uguale a n_ref_cyto
% number of refinements of the mesh for the diffusion problem on the incised disc
% it must be greater than or equal to n_ref_cyto
n_ref_id=1;



% DATI RELATIVI ALLA INTEGRAZIONE NEL TEMPO
% DATA FOR TIME INTEGRATION

% integration method
% first component: 0=ode method (slow); 1=theta method
% second component (only for theta): 0=no factorization; 1=LU factorization; 2=cholesky factorization
% third component (only for theta and factorizations): 0=no permutation; 1=Symmetric minimum degree permutation; 2=Symmetric reverse Cuthill-McKee permutation
% integration method for the homogenized problem

% homogenized problem
metodo_cyto=[1, 1, 1];

% diffusion problem on the incised discs
metodo_id=[1, 1, 1];

% tipo di solutore per il metodo ode
% type of built-in solver invoked if ode method is chosen
solver='ode45';

% parametro theta nell'integrazione temporale
% theta-parameter of the theta-method
theta=0.5; % []

% parametro di rilassamento nell'iterazione sul punto fisso: deve essere 0<alpha<=1 (alpha=1 per eliminare il rilassamento)
% relaxation parameter in the fixed-point iteration: it must be 0<alpha<=1   (alpha=1 means no relaxation)
alpha=1; % []

% errore relativo sul punto fisso nel prolema homogeneizzato
% relative error on the fixed-point iteration of the homogenized problem 
tol_fix=1e-5; % []

% tipo di norma utilizzata per la convergenza (L^\infty (se true), oppure L^1 (se false))
% type of norm used to estimate convergence: (L^\infty (if true), or L^1 (if false))
norma_inf=false;

% flag per la modalità di costruzione della delta
%   vero: delta del tipo M\V: esatta, può dare densità negative
%   falso: delta del tipo k*V: approssimata, ma dà densità positive
% flag for the construntion method of delta
%   true: type M\V: exact, but negative densities may arise
%   false: type k*V: approximate, but always yield positive densities
peak_delta=true;

% durata della simulazione
% time horizon
% t_fin=0.1; % [s]
t_fin=3; % [s]
% t_fin=8; % [s]

% numero di passi di integrazione
% number of integration steps
% n_step_t=50;
% n_step_t=130;
 n_step_t=1000;
% n_step_t=10;
% n_step_t=1000;

% sottocampiona: registra la soluzione ad un istante temporale ogni downsample istati temporali
% utile per ridurre la memoria necessaria per registrare sol
% ricostruisce la soluzione interpolando linearmente
% downsampling: save the solution at one time step each downsample time steps
% useful to reduce the amount of memory required to store sol
% reconstruct by linear interpolation
downsample=1;

% downsample deve essere divisore di n_step_t: cambia di conseguenza n_step_t
% downsample must divide n_step_t: accordingly, change n_step_t
n_step_t=downsample*ceil(n_step_t/downsample);

% FLAG DI DISEGNO E CONTROLLO
% PLOT AND FLOW-CONTROL OPTIONS

% flag per il disegno della mesh
% plotting meshes if true
plot_mesh=false;

% flag per il disegno dei numeri sulla mesh
% drawing nodal and element number on meshes if true
plot_num=false;

% flag per il disegno della cascata delle particelle nella diffusione 
% monte carlo o delle soluzione E* e R* dalla pde
% drawing waterfall of particles if true (monte carlo) or distribution of E*
% and R* solutions
plot_niagara=false;

% flag per il disegno del biliardo (random path)
% drawing pool and beep if true (random path) 
plot_pool =false;

% flag per l'esecuzione senza pause (se inspect=true, l'utente ha maggior controllo del flusso)
% continuous execution if false
% requesting to press keys during execution if =true, for the sake of controlling outputs
inspect=false;


% DATI INIZIALI DI TENTATIVO PER LO STEADY-STATE
% INITIAL VALUES TO BE USED IN SEARCHING FOR THE STEADY-STATE

% valore iniziale di tentativo del cGMP
% trial cGMP value
u_tent=3;  % [\mu M]

% valore iniziale di tentativo  del Ca2+ 
% trial Ca2+ value
v_tent=1; % [\mu M]

% tolleranza sui flussi da rendere nulli per calcolare lo stato stazionario
% maximum unbalancing of fluxes allowed in the steady-state 
tol_stat=1e-8; % [(\mu m) (\mu M)/s = 10^(-9) mole/(m^2 s)]


% SCELTA DEL MODELLO NEL CITOSOL
% flag_model=1 modello 3D omogeneizzato, flag_model=2 modello well
% stirred con diffusione lungo z, flag_model=3 modello ccompletamente
% well stirred
flag_model=1;

% flag_Ca_clamp=false: modello standard; flag_model=true modello con clamp del calcio
% flag_Ca_clamp=false: standard model; flag_model=true: Calcium clamp model
flag_Ca_clamp=false;

% SCELTA DEL MODELLO SU DISCO
% CHOICE OF THE DIFFUSION MODEL ON THE DISC
% flag_model_disc=1 modello R*-->E*
% flag_model_disc=2 modello R*-->G*-->E*
flag_model_disc=2;


% DATI SULLA R*
% DATA ON R*

% costante di decadimento di R*
% decay rate of R*
% k_R=2.8; % [s^{-1}]
k_R=2.5;

% modalità di spegnimento della R*: 
%   0=sempre accesa: in tal caso k_R è ignorato;
%   1=spegnimento brusco in n_R_step passi, con istanti di transizione 
%     deterministici specificati in Phi{d}(3:end,j)
%   2=spegnimento brusco, con istanti di transizione aleatori, 
%     ottenuti da una distribuzione, somma di n_step_R intervalli, 
%     ciascuno distribuito con legge esponenziale con media assegnata da tau
%     gli intervalli aleatori sono in genere differenti per differenti R*
%     pur presentando le stesse distribuzioni statistiche
%   3=decadimento esponenziale: 
%     posto n=n_step_R, 
%     la massa della R* che ha subito i-1 transizioni, i=1..n, 
%     e quindi è attiva in stadio R*_i, è:
%     [(n*t/tau)^{i-1}/(i-1)!] * exp(-n*t/tau) 
%     se gli step di spegnimento hanno tutti la stessa durata (equal_step==true)
%     altrimenti, vedi check_exp_R_st.m
% type of shuting off of R*
%   0=always on: k_R is not used; 
%   1=sudden shuting off in n_step_R steps, with deterministic 
%     transition times, specified in Phi{d}(3:end,j)
%   2=sudden shuting off, at a random instant obtained as the sum of 
%     n_step_R random intervals, each distributed according to an
%     exponential random law with average given by tau
%     the random instants are different for different R*
%     though their statistical distributions are identical
%   3=exponential decay:
%     set n=n_step_R, 
%     the mass of R* which underwent i-1 transitions, i=1..n, 
%     and thus is active in the i-th condition, R*_i, is:
%     [(n*t/tau)^{i-1}/(i-1)!] * exp(-n*t/tau), 
%     provided that the steps take, in the mean, the same time (equal_step==true)
%     otherwise, see check_exp_R_st.m
% mode_time=0;
% mode_time=1;
% mode_time=2;
 mode_time=3;

% modalità di movimento della R*
%   0=fissa alla locazione specificata da Phi{d}(:,j), per il disco d,
%     locazione j
%   1=random path con diffusività D_R_st a partire dalla locazione
%     Phi{d}(:,j) o generata a random, in dipendenza di random_location
%   2=soluzione dell'equazione del calore con diffusività D_R_st 
%     e dato iniziale delta alla locazione Phi{d}(:,j)
% type of spreading of R*
%   0=fixed at location specified by Phi{d}(:,j), for disc d,
%     location j
%   1=random path with diffusivity D_R_st starting from location
%     Phi{d}(:,j), or randomly generated, depending on random_location
%   2=solution of the heat equation with diffusivity D_R_st 
%     and initial datum delta at location Phi{d}(:,j)
mode_space=0;
% mode_space=1;
% mode_space=2;

% numero di passi di shutoff della R*
% questo parametro è usato se mode_time==1,2,3
% number of steps required for a complete shut off of R*
% this parameter is used when mode_time==1,2,3
n_step_R=7;



% Generazione delle stringe lambda, mu e nu_RG date dalle formule di Lixin

lambda_0=2;
mu_0=10;
K_v=0.41;

nu_RG_0=185;

%mu=[0.0000; 0.0000; 0.0000; .7500; 9.0000; 9.0000; 9.0000];
mu=[0 0 0 1 1 1 1];
%mu=[0.1 0.3 0.6 0.8 1 1 1];

mu=mu(1:n_step_R)*mu_0;
%mu=zeros(1,n_step_R);
lambda=lambda_0*(n_step_R-(1:n_step_R));
nu_RG=nu_RG_0*exp(-K_v*((1:n_step_R)-1));


%%%%%%%%%%%%%%%%%%%%%
    if length(lambda)~=n_step_R|length(mu)~=n_step_R|length(nu_RG)~=n_step_R
        error('lambda must have n_step_R components')
    end
% % %     if abs(det(vander(tau)))<1000*eps
% % %         error('Error in the data: in the case of different steps, no two equal step are allowed')
% % %     end



% Coefficiente nu_{RE}, pari al tasso di produzione di E* per ogni R*
% si tratta di un vettore colonna lungo n_step_R, in quanto ad ogni passo di
% fosforilazione di R* corrisponde un diverso tasso di produzione di E*
% usato quando flag_model_disc==1
% rate of production of E* per R* per second
% this is a column vector with n_step_R components: indeed, at each
% phosphorylation level, the R* has a different catalytic activity
% used when flag_model_disc==1
% nu_RE=183*[1]; % [s^{-1}]
% nu_RE=183*[1; 1/2]; % [s^{-1}]
nu_RE=nu_RG; % [s^{-1}]
% nu_RE=4/3*183*[1; 1/2]; % [s^{-1}]
if length(nu_RE)~=n_step_R
    error('nu_RE must have n_step_R rows')
end



% capacità areale della R*
% areal capacity of R*
cc_R_st=1; % []

% coefficiente di diffusione della R*
% diffusion coefficient of R*
D_R_st=0.7;   % [(\mu m)^2/s]
% D_R_st=10;   % [(\mu m)^2/s]
% D_R_st=0;   % [(\mu m)^2/s]

% numero di locazioni di fotoisomerizzione su ciscun disco speciale: la lunghezza del vettore n_Phi deve essere n_sd
% su ciascun disco speciale deve esserci almeno una locazione di fotoisomerizzazione
% number of position where fotoisomerizations occur on each special disc: the length of n_Phi must be n_sd
% on each special disk at least one position must be present
n_Phi=[1];
% n_Phi=[3];
% n_Phi=[4,1];
if length(n_Phi)~=n_sd
    error('n_Phi must have n_sd components')
end

% le posizioni e i tempi di spegnimento di ciascuna locazione di fotoisomerizzazione, 
% per ciascun disco speciale, sono nella cell structure Phi
% Essa deve avere n_sd componenti Phi{s}, s=1..n_sd, una per ogni disco speciale
% Phi{s} è una matrice con 2+n_step_R righe e n_Phi(s) colonne, 
% una per ogni locazione di fotoisomerizzazione
% su ogni colonna, sono presenti:
%   rho; theta, che sono le coordinate polari della locazione di fotoisomerizzazione
%   nota: (rho,theta) sono utilizzati solo se random_location==false
%   t_1, t_2, ... t_{n_step_R} che sono gli istanti di transizione da uno
%   stato al successivo degli n_step_R stati possibili
%   nota: questi istanti di transizione sono utilizzati solo se
%   mode_time=1 (spegnimento improvviso deterministico)
    
% position and transition times for each location of photoisomerizaztion, 
% for each special disk, are in the cell structure Phi
% Phi must have n_sd components Phi{s}, s=1..n_sd, one for each special disc
% Phi{s} is a matrix with 2+n_step_R rows and n_Phi(s) columns, 
% one for each location of photoisomerization
% on each column there are: 
%   rho; theta, which are the polar coordinates of that location; 
%   note: (rho,theta) are used only when random_location==false
%   t_1, t_2, ... t_{n_step_R}, which are the transition times from one
%   state to the next of the n_step_R admissible states
%   note: these transition times are used only when mode_time=1 
% (deterministic sudden shutoff)
Phi=cell(1,n_sd);

Phi(1)={[1*R; pi; cumsum(1./(lambda+mu))']};
if length(Phi)~=n_sd
    error('Phi must have n_sd components')
end
if length(Phi)~=n_sd
    error('Phi must have n_sd components')
end
for s=1:n_sd
    if size(Phi{s},1)~=2+n_step_R
        error('Phi(s) must have 2+n_step_R rows')
    end
    if size(Phi{s},2)~=n_Phi(s)
        error('Phi(s) must have n_Phi(s) columns')
    end
end

% se true, la posizione iniziale delle fotoisomerizzazioni è scelta
% randomm, differente per ogni campione
% if true, the initial positions of photoisomeritazions is randomly chosen,
% different for each sample
random_location=false;


% DATI SULLA G*
% DATA ON G*

% tasso di produzione di E* per ogni G* per densità unitaria (#/\mu m^2) di E spenta
% production rate of E* per each G* per unit density (#/\mu m^2)of basal E
% k2 in shraiman paper
k_GE=1; % [\mu m^2/s]

% capacità areale della G*
% areal capacity of G*
cc_G_st=1; % []

% coefficiente di diffusione della G*
% diffusion coefficient of G*
D_G_st=1.5+0*0.7;   % [(\mu m)^2/s]
% D_G_st=1.5+0.7;%+1.5;  
% D_G_st=1.5;%+1.5;


% PARAMETRI MICROSTRUTTURALI
% MICROSTRUCTURAL PARAMETERS

% spessore dei dischi
% disc thickness
epsilon_0=14e-3; %  [\mu m]

% rapporto fra spessore degli interdischi e spessore dei dischi
% ratio between interdisc and disc thickness 
nu=14/14; % []

% rapporto fra spessore dell'outer shell e spessore dei dischi
% ratio between outer-shell and disc thickness
sigma=15/14; % []


% DATI SULLA E*
% DATA ON E*

% Tasso di idrolisi del cGMP da parte della PDE attiva al buio (spenta)
% rate of hydrolysis of cGMP by dark-activated PDE
k_hyd=7e-5;  % [(\mu m)^3/(molecole s)]

% Densità superficiale di PDE totali
% nota: la densità di subunit è il doppio di PDE_s
% superficial density of total PDE
% note: subunit density is 2*PDE_s
PDE_s=100;  % [molecole/(\mu m)^2]

% % attivazione della PDE spenta = k_hyd*[PDE]_s = flusso
% gamma_0= k_hyd*PDE_s;  %  [\mu m s^(-1)]

% Costante di idrolisi del cGMP da parte della PDE attivata
% rate of hydrolysis of cGMP by light-activated PDE
k_st=0.5; % [(\mu m)^3/molecole s]

% capacità areale della E*
% areal capacity of E*
cc_E_st=1; % []

% coefficiente di diffusione della E*
% diffusion coefficient of E*
% D_E_st=0.7+1.5+0.8;   % [(\mu m)^2/s]

D_E_st=0.8;%+0.5;   % [(\mu m)^2/s]

% D_E_st=5;  % [(\mu m)^2/s]
% D_E_st=10;  % [(\mu m)^2/s]

% costante di decadimento di E*
% decay rate of E*
% kh in shraiman paper
% k_E=0.64;%0.58;%(0.58+0.76)/2; % [s^{-1}]
k_E=0.58;
% nota: Shraiman dà 10
%  k_E=0.2;%(0.58+0.76)/2; % [s^{-1}]


% TIPO DI SOLUZIONE: PDE vs MONTECARLO E*
% TYPE OF SOLUTION: PDE vs MONTECARLO E*

% tipo di soluzione del problema sui dischi incisi:
%   vero: Monte Carlo sia per R* che per E*; in tal caso, ignora mode_space:
%         nello spazio c'è random walk sia per R* che per E* (mode_space==1)
%         le particelle sono sempre accese o lo spegnimento è sempre random (mode_time==0,2)
%         attivo solo se flag_model_disc==1
%   falso:diffusione con spegnimento esponenziale per E*;
%         il trattamento di R* dipende da mode_space e mode_time
% solution type for the incised-disc problem
%   true: Monte Carlo both for R* and for E*; in this case, neglect
%         mode_space: always random walk both for R* and E* (mode_space==1)
%         particles are always on or random shut off (mode_time==0,2)
%         only allowed when flag_model_disc==1
%   false:diffusion (PDE) for E* with exponential shut off
%         mode_space and mode_time determine the tretment of R*
MC_disk=false;

% se le R* diffondono, non possono avere spegnimento brusco, né
% deterministico, né random
% se le R* fanno un random walk, non possono avere spegnimento esponenziale
% se MC_disk==true, deve essere mode_space==1, mode_time==0,2, flag_model==1, flag_model_disc==1
% if R* diffuse, they cannot suddenly shut off, neither randomly, neither
% deterministically
% if R* do a random walk, they cannot exponentially fade
% if MC_disk==true, it must be mode_space==1 and mode_time==0,2, flag_model==1, flag_model_disc==1
if ((mode_space==2)&&((mode_time==1)||(mode_time==2))) || ((mode_space==1)&&(mode_time==3)) || ...
        ( (MC_disk) && ((mode_space~=1)||((mode_time~=0)&&(mode_time~=2))||(flag_model~=1)||(flag_model_disc~=1)))
    error('combination of MC_disk, mode_space, mode_time, flag_model not allowed')
end


% numero di campioni statistici
% number of random samples
n_sample=1;
if (n_sample>1) && ~((mode_space==1)||(mode_time==2)||(random_location))
    error('random samples not necessary for a deterministic setting')
end

% COEFFICIENTI DI DIFFUSIONE NEL CITOSOL
% CITOSOL DIFFUSION COEFFICIENTS

% capacità volumica del cGMP
% volumic cGMP capacity
cc_u=1; % []

% coefficiente di diffusione volumico del cGMP 
% volumic cGMP diffusion coefficient
kk_u=160;%(50+196)/2;%100;  % [(\mu m)^2/s] 

% capacità volumica del Ca2+
% volumic Ca2+ capacity
cc_v=1; % []

% coefficiente di diffusione volumico del Ca2+ 
% volumic Ca2+ diffusion coefficient
kk_v=15;  % [(\mu m)^2/s] 


% DATI SULLA CICLASI
% DATA ON CYCLASE
 
% Tasso massimo di sintesi di cGMP
% maximum rate of production of cGMP
alpha_max=50; % [\mu M s^{-1}]

% Tasso minimo di sintesi di cGMP
% minumum rate of production of cGMP
alpha_min=0.02*alpha_max; % [\mu M s^{-1}]

% Esponente di Hill
% Hill coefficient
m_cyc=2; % []

% Concentrazione di calcio corrispondente alla metà del tasso massimo di attivazione
% half-maximum-activation concentration
k_cyc=0.135; % [\mu M]

% gamma_1=(alpha_max-alpha_min)*k_cyc^m_cyc*(1/2*nu*epsilon_0); 
% 
% beta_1=k_cyc; % [\mu M]
% 
% m=m_cyc;
% 
% cf1=alpha_min*(1/2*nu*epsilon_0); % [\mu M \mu m s^{-1}]


% DATI COMUNI AI CANALI DEL CALCIO
% COMMON DATA ON CALCIUM CHANNELS

% Buffer del Ca2+, assunto costante
% buffer of Ca2+, here assumed to be constant
B_ca=20; % []

% Costante di Faraday
% Faraday constant
F=96500/1e21; % [C / (\mu M (\mu m)^3)]


% DATI SUL CANALE cGMP-OPERATO
% DATA ON THE cGMP-OPERATED CHANNEL

% Corrente di scambio massima
% maximum exchange current
j_cg_max=7000e-12; % [A]

% Frazione di corrente attivata dal c_GMP costituita da ioni Ca2+
% fraction of cGMP-activated current given by Ca2+
f_ca=0.17; % []

% Esponente di Hill
% Hill coefficient
m_cg=2.5; % []

% Concentrazione di cGMP alla metà dell'apertura dei canali
% half-maximum-activation concentration, here assumed to be constant
K_cg=20;%32; % [\mu M]

% c_2=j_cg_max*f_ca/(Sigma_rod*B_ca*F*2); % [\mu M \mu m s^{-1}]
% 
% d_2=K_cg;
% 
% k=m_cg;
% 


% DATI SULLO SCAMBIATORE
% DATA ON Ca2+ EXCHANGER

% Corrente di scambio alla saturazione
% exchange current at saturation
j_ex_sat=17e-12; % [A]

% Concentrazione di calcio alla metà dell'apertura dei canali
% half-maximum-activation concentration
K_ex=1.5; %[\mu M]

% c_1=j_ex_sat/(Sigma_rod*B_ca*F); % [\mu M \mu m s^{-1}]
% 
% d_1=K_ex; % [\mu M]

return
