% Filename: PDE_Solver.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, April 2020
% This function generates the solutions to Equations (8) and (9) with the 
% following parameters, and plot the density profiles at  t=20,30,40
% L_u = 0, L_v = L = 60, L(0) = 30, beta = 1, alpha = 0.5
% D = 1, lambda = 1, kappa_u = 1.2195, kappa_v = 0.5, dt = 0.001, 
% dxi = deta = 0.00025.
% A second figure is generated with the position of the moving boundary 
% in function of time L(t).

% relative density, Equation (13)
D = 1;
% relative proliferation rate, Equation (13)
lambda = 1;
% kappa_u, Equation (13)
kappa_u = 1.2195;
% kappa_v, Equation (13)
kappa_v = 0.5;
% \Delta t, Equation (16) and (17)
dt = 0.001;
% \Delta \xi, Equation (16)
dxi = 0.00025;
% \Delta \eta, Equation (17)
deta = 0.00025;
% tolerance \epsilon used in Newton-Raphson algorithm
tol = 1e-08;
% total length of the domain
L = 60;
% total time
total_time = 40;
% total number of steps
ts = round(total_time/dt+1);

% initial position of the moving boundary s0=L(0)
s0 = 30;
% next position of the moving boundary L(t)
s_t = s0;
% current position of the moving boundary L(t)
s_tp = s_t;
% array of L(t) for all time steps
st_array = zeros(1,ts); 

% spatial domains dicretised, one for each population
xi = 0:dxi:1;
eta = 1:deta:2;
% number of nodes in each spatial domain
nodes_xi = size(xi,2);
nodes_eta = size(eta,2);

% initialisation of variables used in Newton-Raphson algorithm
% correction of densities at each iteration
delta_u = ones(1,nodes_xi);
delta_v = ones(1,nodes_eta);

% function F
Fu = zeros(1,nodes_xi);
Fv = zeros(1,nodes_eta);
% current densities u(x,t) and v(x,t) at time t = current time
u_p = zeros(1,nodes_xi);
v_p = zeros(1,nodes_eta);

% coefficients a b c of the tridiagonal matrix
% Jacobian
% J(u)
coeffA_u = zeros(1,nodes_xi);
coeffB_u = zeros(1,nodes_xi);
coeffC_u = zeros(1,nodes_xi);
% J(v)
coeffA_v = zeros(1,nodes_eta);
coeffB_v = zeros(1,nodes_eta);
coeffC_v = zeros(1,nodes_eta);

%colors for plots
colors = [0 230/255 187/255; 232/255 206/255 0;];

% Initialisation of densities u(x,0) and v(x,0)
% \beta and \alpha from Equations (14) and (15)
beta = 1; 
alpha = 0.5;
% phi(x) in Equation (14)
for i = 1:nodes_xi
    x = xi(i)*s0;
    if x < beta
        u_p(i) = alpha;
    else
        if x < s0
            u_p(i) = alpha*(1-(x-beta)/(s0-beta));
        else
            u_p(i) = 0;
        end
    end
end
% psi(x) in Equation (15)
for i = 1:nodes_eta
    x = (eta(i)-1) * (L-s0) + s0;
    if x > L-beta
        v_p(i) = alpha;
    else
       if x  > 0
            v_p(i) =  alpha*(x-s0)/(L-s0-beta);
       else
            v_p(i) = 0;
       end
    end
end

% next time step densities u(x,t) and v(x,t)
u = u_p;
v = v_p;

figure
hold on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
%% Newton-Raphson algorithm - main loop

% for all time steps
for j = 1:ts
    %current time
    t = j * dt

    st_array(j) = s_t;
    condition = 1;

    % while the tolerance is not reached
    while (condition)

        % Equations (18) and (19)
        % boundary conditions for u(x,t)
        % at xi = 0
        coeffA_u(1,1) = 0.0;
        coeffB_u(1,1) = -1.0;
        coeffC_u(1,1) = 1.0;
        Fu(1) = -1.0*(u(1,2) - u(1,1));
        % at xi = 1
        coeffA_u(1,nodes_xi) = 0;
        coeffB_u(1,nodes_xi) = 1.0;
        coeffC_u(1,nodes_xi) = 0;
        Fu(1,nodes_xi) = -1.0 * (u(1,nodes_xi) );

        % J(u) delta u = -F(u)
        % Equation (16)
        for i = 2:nodes_xi-1
           coeffA_u(1,i) = 1.0/(dxi^2*s_t^2) - (i-1)*dxi/s_t * (s_t-s_tp)/(2*dt*dxi);
           coeffB_u(1,i) = - 2.0/(dxi^2*s_t^2)  - 1.0/dt + 1.0 - 2.0*u(1,i);
           coeffC_u(1,i) = 1.0/(dxi^2*s_t^2) + (i-1)*dxi/s_t * (s_t-s_tp)/(2*dt*dxi);
           Fu(1,i) = -(u(1,i+1) - 2*u(1,i) + u(1,i-1))/(dxi^2*s_t^2) ...
               - (i-1)*dxi/s_t * (u(1,i+1) - u(1,i-1)) * (s_t-s_tp)/(2*dt*dxi) ...
               + (u(1,i)-u_p(1,i)) / dt - u(1,i)*(1.0-u(1,i));
        end 
        delta_u = tridia(coeffA_u, coeffB_u, coeffC_u, Fu, nodes_xi);

        % correction of u
        for i = 1:nodes_xi
            u(1,i) = u(1,i) + delta_u(1,i);
        end

        % boundary conditions for v(x,t)        
        % Equations (18) and (19))
        % at eta = 1
        coeffA_v(1,1) = 0;
        coeffB_v(1,1) = 1.0;
        coeffC_v(1,1) = 0.0;
        Fv(1) = -1.0*(v(1,1));
        % at eta = 2
        coeffA_v(1,nodes_eta) = -1.0;
        coeffB_v(1,nodes_eta) = 1.0;
        coeffC_v(1,nodes_eta) = 0;
        Fv(1,nodes_eta) = -1.0*(v(1,nodes_eta)-v(1,nodes_eta-1)); 

        % Here L is the size of both domain
        M = L-s_t;

        % J(v) delta v = -F(v)
        % Equation (17)
        for i = 2:nodes_eta-1
           coeffA_v(1,i) = D/(deta^2*M^2) - (2-((i-1)*deta+1))/M * (s_t-s_tp)/(2.0*dt*deta);
           coeffB_v(1,i) = - 2.0*D/(deta^2*M^2)  - 1.0/dt + lambda * (1.0-2.0*v(1,i));
           coeffC_v(1,i) = D/(deta^2*M^2) + (2-((i-1)*deta+1))/M * (s_t-s_tp)/(2*dt*deta);
           Fv(1,i) = -D*(v(1,i+1) - 2*v(1,i) + v(1,i-1))/(deta^2*M^2) ...
               - (2-((i-1)*deta+1))/M * (v(1,i+1) - v(1,i-1)) ...
               * (s_t-s_tp)/(2*dt*deta) + (v(1,i)-v_p(1,i)) / dt - lambda*v(1,i)*(1-v(1,i)) ;
        end   

        delta_v = tridia(coeffA_v, coeffB_v, coeffC_v, Fv, nodes_eta);

        %correction of v
        for i = 1:nodes_eta
            v(1,i) = v(1,i) + delta_v(1,i);
        end

        % verify if tolerance is reached
        if (norm(delta_u,Inf) <= tol && norm(delta_v,Inf) <= tol)
           condition = 0;
        end
        % Stefan condition used to calculate L(t) from Equation (20)
        s_t = s_tp + dt*(-kappa_u*(u(1,nodes_xi)-u(1,nodes_xi-1))/(dxi*s_t) - kappa_v*(v(1,2)-v(1,1))/(deta*M));
    end
    % updating current L(t)
    s_tp = s_t;
    % updating current u(x,t) and v(x,t)
    u_p = u;
    v_p = v;

    % plotting curves at t = 20,30,40
    if (t == 20 || t == 30 || t == 40)
        plot(xi(1:nodes_xi)*s_t, u, 'k-','Color',colors(1,1:3), 'LineWidth',2 ,'DisplayName', '$u(x,20)$');

        plot((eta(1:nodes_eta)-1)*(L-s_t)+s_t, v,'r-','Color',colors(2,1:3),'LineWidth',2,'DisplayName', '$v(x,20)$' );;
    end
end
%%
% compute wave speed
c = (st_array(end) - st_array(end-1))/dt;
% display c
textc = strcat(strcat('$c = ',num2str(round(c,2))), '$');
text(40,0.7,textc,'interpreter','latex','fontsize',18)

ylabel('$u(x,t)$ \ \ $v(x,t)$','interpreter','latex','fontsize',18);
xlabel('$x$','interpreter','latex','fontsize',18);
box on
hold off

figure
% plotting L(t)
plot((1:ts)*dt, st_array,'Color','m','LineWidth',2,'DisplayName', 'L(t)' );
ylabel('$L(t)$','interpreter','latex','fontsize',18);
xlabel('$t$','interpreter','latex','fontsize',18);
xlim([0 40]);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
box on
%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are the three diagonals of the matrix A. N is the size of 
% the vector solution x.
function x = tridia(a,b,c,d,N)
    x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end

    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end 
