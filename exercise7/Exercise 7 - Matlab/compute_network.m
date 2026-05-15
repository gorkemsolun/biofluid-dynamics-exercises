function [p,c] = compute_network(p,c,t,invivo)

nc = size(c,1); %number of connections
np = size(p,1); %number of nodes

% Set physical and geometrical parameters
a = c(:,8);         % capillary radius
mu = 1.5e-3;        % plasma dyn. viscosity
L = c(:,6);         % pipe length
u0 = 1e-3;          % inflow velocity
q0 = u0*pi*a(1)^2;  % mass flow rate at inlet


%% Compute transmissibilities

% Base transmissibility (plasma Poiseuille flow)
c(:,3) = a.^4*pi./(8*mu*L);

% Calculate the apparent viscosity by a formula from Pries et al. (2005)
mueffrel=compute_relative_effective_viscosity(2*c(:,8)*1e6,c(:,9),invivo);

% Update the transmissibility
c(:,3) = c(:,3)./(real(mueffrel));
%c(:,3) = c(:,3)./(1+KT*c(:,7));


%% Set pressure matrix
P = zeros(np); % initialize the pressure matrix with [np x np] zeros

% Loop through all vessels
for i = 1:nc

    ip0 = c(i,1); % p0 of current connection
    ip1 = c(i,2); % p1 of current connection
    T = c(i,3);   % Transmissibility

    P(ip0,ip0) = P(ip0,ip0) - T;
    P(ip0,ip1) = P(ip0,ip1) + T;

    P(ip1,ip0) = P(ip1,ip0) + T;
    P(ip1,ip1) = P(ip1,ip1) - T;
end


% Eliminate last node from system (pressure=0)
P = P(1:np-1,1:np-1);

%% Set right-hand side

rhs = -q0*eye(np-1,1); % eye = identity matrix

%% Solve for pressures

p(1:np-1,3) = P\rhs;

% Pressure outlet condition
p(np,3) = 0; 

%% Compute bulk flow velocities

% Loop through all vessels
for i = 1:nc
    
    ip0 = c(i,1); % p0 of current connection
    ip1 = c(i,2); % p1 of current connection
    T = c(i,3);   % Transmissibility

    c(i,5) = T*(p(ip0,3)-p(ip1,3))/(pi*a(i)^2); % T * deltaP / A

end

%% Compute average flow rates at nodes

% Inflow rates
p(:,5) = 0; % set all flow rates at the nodes to zero

% Loop through all vessels
for i = 1:nc
   if (c(i,5)<0)
       % Update flow at p0
       p(c(i,1),5) = p(c(i,1),5) + abs(c(i,5)*pi*a(i)^2);
   end
   if (c(i,5)>0)
       % Update flow at p1
       p(c(i,2),5) = p(c(i,2),5) + abs(c(i,5)*pi*a(i)^2);
   end
end

% Set the mass inflow condition at the inlet node
p(1,5) = p(1,5)+q0; 

end