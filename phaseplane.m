% Filename: phaseplane.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, April 2020
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2020) 
% A sharp-front moving boundary model for malignant invasion.
% The script contains:
%   - two calls to the function pptwopop to generate 
%     Figures 3(b) and (c).
%	- the function pptwopop
%   - the function HeunsSolver

% Generating Figure 3(b)
% Calling the function pptwopop(c,D,lambda,dz,z_begin,z_end) with the
% following inputs: c = 0.2, D = 1, lambda = 1, dz = 0.001, z_begin = 0,
% z_end = 50
pptwopop(0.2,1,1,0.001,0,50);
% Generating Figure 3(c)
% Calling the function pptwopop(c,D,lambda,dz,z_begin,z_end) with the
% following inputs: c = 0.2, D = 1, lambda = 1, dz = -0.001, z_begin = 50,
% z_end = 0
pptwopop(0.2,1,1,-0.001,50,0);

% Function pptwopop
% This function solves Equations (26) and (27) in the phase plane
% by Heun's method and plots the solution on the plane X(z) versus V(z), or 
% W(z) versus U(z). The same plot shows also the equilibrium points (0,0) 
% and (1,0) and the intersection point of the solution with the 
% corresponding axis U(z)=0 or V(z)=0.
% INPUT ARGUMENTS:
% ** c, the wave speed, positive or negative, 0 < |c| < 2
% ** D, the relative diffusity, D > 0
% ** lambda, the relative proliferation rate, lambda > 0
% ** dz, the step size used to discretise the domain of z. The choice of dz
% must take into account the direction of the trajectory in the phase plane
% used when integrating Equations (26) and (27). Starting at the initial
% conditions, if the next value of the solution is in the positive 
% direction of z, then dz > 0 and z_begin < z_end, it the next value is in
% the negative direction of z, then dz < 0 and z_begin > z_end. 
% ** z_begin and z_end, the lower and upper limit of the domain of z, used
% to integrate numerically Equations (26) and (27) by Heun's method,
% such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
function pptwopop(c,D,lambda,dz,z_begin,z_end)
   
    % Using the values of c, D and lambda, setting the matrix of the
    % system of Equations (26) and (27) linearised around the
    % saddle point (1,0).
    % Finding the eigenvectors and the eigenvalues of the solutions. 
    if (dz > 0)
        DV = 1;
        lambdaV = 1;    
    else 
        DV = D;
        lambdaV = lambda;
    end
    % Value of V(z) around the saddle point
    Vs = 1;
    % Matrix of the linearised system of equations
    A = [0 1;(2*Vs-1)*lambdaV/DV -c/DV];
    % Multiplying the matrix by -1 if the integration follows the backward
    % direction of z
    if (dz < 0)
        A = -A;
    end
    % Finding the eigenvalues d and eigenvectors v of the linearised system
    [v,d]=eig(A);
    % Defining the initial conditions
    % Choosing the eigenvector where to start the solution whether the 
    % solution follows the stable or the unstable manifolds, in the forward
    % or the backward direction of z.
    IC = [0;0];
    if (d(1,1) > 0)
        IC = v(:,1);
    end
    if (d(2,2) > 0)
        IC = v(:,2);
    end
    if (c < 0 && dz > 0)
       IC(2,1) = -IC(2,1);
    end
    if (dz < 0 && (DV ~= 1 || lambda ~= 1))
        IC(1,1) = -IC(1,1);
        IC(2,1) = -IC(2,1);
    end
    % Setting the initial conditions
    IC1 = Vs+IC(1,1)*0.001;
    IC2 = IC(2,1)*0.001;

    % Solving Equations Equations (21) and (22) in the phase plane with
    % Heun's method
    [V, X] = heunSolver(c, DV, lambdaV, dz, z_begin, z_end, IC1, IC2);
     
    % Finding the intersection point with the corresponding axis U(z)=0 or 
    % V(z)=0
    intersect = 0;
    intersects = [];
    for ii = 1:length(V)-1
        if (V(ii,1)>= 0 && V(ii+1,1)<= 0 || V(ii,1)<= 0 && V(ii+1,1)>= 0 )
            var1 = [V(ii,1) V(ii+1,1)];
            var2 = [X(ii,1) X(ii+1,1)];
            intersects = [ intersects interp1(var1,var2,0,'linear')];
        end
    end
    if (dz > 0)
        intersect = min(intersects);
    else
        intersect = max(intersects);
    end

    % Setting and computing the field vectors of the solution
    y1 = linspace(-1.2,1.6,14);
    y2 = linspace(-1.2,1.6,14);
    [x,y] = meshgrid(y1,y2);
    dv = zeros(size(x));
    dx = zeros(size(x));
    if c < 0
        f1 = @(Y) [-Y(2);-(-c*Y(2)+lambda*(-Y(1)+Y(1).*Y(1)))/D;];
    else
        f1 = @(Y) [Y(2);(-c*Y(2)+lambda*(-Y(1)+Y(1).*Y(1)))/D;];
    end
    for i = 1:numel(x)
        Yprime = f1([x(i); y(i)]);
        dv(i) = Yprime(1);
        dx(i) = Yprime(2);
    end
    
    % Defining the colors - orange and green
    colors = [227/255 92/255 0; 0 117/255 94/255;];
    figure
    % Plotting the corresponding axis of the phase plane: U(z) and W(z), 
    % or V(z) and X(z)
    line([-1.5 1.5],[0 0],'Color','k','LineStyle','-','LineWidth',2);
    hold on
    line([0 0],[-1.5 1.5],'Color','k','LineStyle','-','LineWidth',2);
    % Plotting the solution 
    if (dz > 0)
        % W(z) and U(z)
         plot(V,X,'r:','LineWidth',2,'Color',colors(1,1:3));
    else
        % V(z) and X(z)
         plot(V,X,'r:','LineWidth',2,'Color',colors(2,1:3));
    end

    % Plotting the field vectors
    quiver(x,y,dv,dx,'b','LineWidth',1); 
    % Plotting the intersection point of the solution with the 
    % corresponding axis U(z) = 0 or V(z) = 0
    plot(0,intersect,'mo','MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',4);
    % Plotting the equilibrium points (0,0) and (1,0)
    plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    plot(1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);

    % Plotting the labels of the corresponding axis
    if (dz < 0)
        ylabel('$X(z)$','interpreter','latex');
        xlabel('$V(z)$','interpreter','latex');
    else
        ylabel('$W(z)$','interpreter','latex');
        xlabel('$U(z)$','interpreter','latex'); 
    end
    xlim([-1,1]);
    ylim([-1.2,1.2]);
    box on;
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 18);
end

% Function heunSolver
% This function solves Equations (21) and (22) by Heun's method 
% as described in the Supplementary material.
% INPUT ARGUMENTS:
% ** c, the wave speed, positive or negative, 0 < |c| < 2
% ** D, the relative diffusity, D > 0
% ** lambda, the relative proliferation rate, lambda > 0
% ** dz, the step size used to discretise the domain of z
% ** z_begin and z_end, the lower and upper limit of the numerical domain 
% of z such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
% ** V1, X1, the values of the initial conditions
% OUTPUT ARGUMENTS:
% ** Vout : The solution V(z)
% ** Xout : The solution X(z)
% The output size of Vout and Xout may not be equal to the original number
% of nodes correponding to (z_end-z_begin)/dz+1 in the discretised domain.
% If the solutions contains values superior to 3 or inferior to -3, 
% the arrays Vout and Xout are truncated such as all the values are between
% -3 and 3.
function [Vout, Xout] = heunSolver(c,D,lambda,dz,z_begin,z_end,V1,X1)

    % The discretised domain of z
    z = z_begin:dz:z_end;
    % Number of nodes in the domain
    sz = length(z);

    % Initialisation of the arrays X(z) and V(z) to zero
    X = zeros(sz,1);
    V = zeros(sz,1);
    % Initial conditions
    V(1) = V1;
    X(1) = X1;
   
    % Initialising the final number of nodes szout using the original
    % size of the domain
    szout = sz;
    
    % For all steps in the domain
    for i = 1:sz-1
        % Equations (26) and (27) solved with Heun's method
        % See Supplementary Material
        Vbar = dz * X(i) + V(i); 
        Xbar = dz * (-c*X(i)- lambda* V(i)*(1-V(i)))/D + X(i);
        V(i+1) = dz/2 * (X(i)+Xbar) + V(i); 
        X(i+1) = dz/2 * ((-c*X(i) - lambda * V(i)*(1-V(i))) + (-c*Xbar ...
            - lambda * Vbar *(1-Vbar)))/D + X(i); 
        % If the solution value is too large, stop the solver
        if (X(i+1) < -3 || X(i+1) > 3 || V(i+1) < -3 || V(i+1) > 3)
            % Modifying the number of nodes used, szout
            szout = i;
            break;
        end
    end
    
    % Truncating the output solutions using the number of nodes szout
    Vout = V(1:szout,1);
    Xout = X(1:szout,1);
end