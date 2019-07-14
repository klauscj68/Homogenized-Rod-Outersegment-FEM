% problema della creazione e diffusione di E^* sui dischi speciali incisi
function [time, E_st]=diffusione_ipd_MC_A(R, n_pd, n_tri, p_pd, t_pd, n_sd,  ...
    t_fin, n_step_t, plot_mesh, inspect)

% output
% E_st{d}(j,:) is the density of E* [moleculs/(\mu m)^2] on the special disc d, at the time instant (j-1)*t_step, 
% in all the nodes of the mesh p_pd of the pivot disk of the homogenized problem

fprintf('\nInterpolating from Monte Carlo simulations\n');

if n_sd~=1
    error('only one special disc allowed at the moment from Monte Carlo simulations')
end
    
% load PDELattice501_28g;
load PDE_Ri75Rj75RDR1.65KE0.58pGp0.035AVG.mat;
PDELattice=PDELatticeAVG;


[n_A_mesh,ignore,n_A_step]=size(PDELattice);
ds_A=(2*R)/(n_A_mesh);


maxD=1.5; %um per second (awful thing to put absolute constants away from the data.m)

dt_A=(2*R)^2/(maxD*n_A_mesh^2);
t_end=(n_A_step-1)*dt_A;

% intermediate V, at Audi's times
V_int=zeros(n_A_step,n_pd);

% Audi's mesh
X=-R+ds_A/2:ds_A:R-ds_A/2;
Y=-R+ds_A/2:ds_A:R-ds_A/2;
T=0:dt_A:t_end;

% cycle over time
for j=1:n_A_step
    % calculate V_int(j,:)
    [I,J,W]=find(PDELattice(:,:,j));
    % numer of nonzero components
    n_nz=length(W);
    % cycle over non zero components
    for k=1:n_nz
        % compute the load due to W(k) deltas (E*) located at ((I(k)-0.5)*ds_A; (J(k)-0.5)*ds_A)

        % cartesia components
        Xp=X(I(k));      % -R+(I(k)-0.5)*ds_A;
        Yp=Y(J(k));      % -R+(J(k)-0.5)*ds_A;

        % find an element tri of the mesh of triangles defined by p_pd,t_pd to which (Xp,Yp) belong
        % and sample at (Xp,Yp) the corresponding shape function
        [tri,F_elem]=find_tri(p_pd,t_pd,[Xp;Yp],R);
        
        % indexes of the nodes of the triangle
        ind=t_pd(:,tri);
        % moltiply by the number of E* at point Xp,Yp and sum on the accumulation vector W(j,:)
        V_int(j,ind)=V_int(j,ind)+W(k)*F_elem';
    end
    % numer of Audi's E* at time step j
    sum(sum(PDELattice(:,:,j)));
    % number of our E* at time step j
    sum(V_int(j,:));
end


% interpolating in time
time=0:t_fin/n_step_t:t_fin;

% final V, at our times
V=zeros(n_step_t+1,n_pd);

% cycle over nodal points
for n=1:n_pd
    V(:,n)=interp1(T, V_int(:,n), time)';
end


% find the nodal values of a piece-wise linear function 
% giving rise to the same load vector computed before

% compute mass matrix
M=sparse(n_pd,n_pd);
for t=1:n_tri
    ind=t_pd(:,t);
    % nodal coordinates
    X=p_pd(1,ind);
    Y=p_pd(2,ind);
    % element matrix
    [M_elem]=tri_M(X,Y);
    % assemble
    M(ind,ind)=M(ind,ind)+M_elem;
end



% convert from load vector to nodal values for densities
E_st{1}=(M\V')'; % delta esatta: può dare valori negativi
% % % E_st{1}=(V'./repmat(ones(1,n_ipd)*M*V',n_pd,1))'; % delta
% approssimata da pagoda

% total numer of E* in time
E_st_tot=sum((M*E_st{1}')',2);
% ones(1,n_pd)*M*E_st{1}'

% stampa la storia delle diesterasi totali, per controllo
fprintf('\nTotal numer of E* in time\n')
massa=E_st_tot';
disp([time; massa]);
fprintf('\nGiratore polare di inerzia delle E*\n');
rho_sq=ones(n_step_t+1,1)*(p_pd(1,:).^2+p_pd(2,:).^2);
inerzia=ones(1,n_pd)*M* (rho_sq.*E_st{1})' ;
% esclude il primo istante, cui corrisponde massa di E* pari a zero
giratore=[0 sqrt(inerzia(2:end)./massa(2:end))];
disp([time; giratore]);
figure(13)
hold on
plot(time, giratore)


if plot_mesh && inspect
	figure(10)
	hold on
	plot(T,sum(V_int,2),'g')
% 	plot(time,sum(V,2),'b')
	plot(time,E_st_tot,'r')
	presskey(inspect);

    % draw the E* density

    % maximum E* density
    maxEst=max(max(E_st{1}));
    
    % one picture every 4 time steps
    for t=1:4:n_step_t+1
		f11=figure(11);
        newplot(f11)
        title('E* density')
        % chromatic scale RGB
        colori=[(1-E_st{1}(t,:)'/maxEst), E_st{1}(t,:)'/maxEst, zeros(n_pd,1)];
		% Mesh of triangles
        vert=[p_pd; E_st{1}(t,:)];
		patch('Vertices',vert','Faces',t_pd','FaceVertexCData',colori,...
              'FaceColor','interp','EdgeColor','interp');
        xlabel('x [\mu m]')
		ylabel('y [\mu m]')
		zlabel('E* [mol/\mu m^2]')
        axis([-R*1.1 R*1.1 -R*1.1 R*1.1 0 maxEst])
        view(3)
	    presskey(inspect);
    end
	presskey(inspect);
end

% save audi
% area_for_each_grid_point=ones(1,n_pd)*M
% plot(T,reshape(sum(sum(PDELattice,1),2),1,n_A_step),time,E_st_tot')

