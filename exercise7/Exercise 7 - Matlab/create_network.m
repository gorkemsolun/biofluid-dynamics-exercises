function [p,c] = create_network(n,Ledge)
%% define network parameters

% typical length of a capillary: Ledge

%% define grid
% x-axis:
gridx1 = zeros(1,n);
gridx2 = zeros(1,n);
for i = 2:length(gridx1)
    if mod(i,2)
        % if i is odd
        gridx1(i) = gridx1(i-1) + 1;
        gridx2(i) = gridx2(i-1) + 2;
    else
        gridx1(i) = gridx1(i-1) + 2;
        gridx2(i) = gridx2(i-1) + 1;
    end
end
gridx1 = gridx1*Ledge;
gridx2 = (gridx2+0.5)*Ledge;

% y-axis:
gridy = linspace(0,(2*n-2),(2*n-1));
gridy = gridy*Ledge*sqrt(3)./2; % sin(60)=sqrt(3)./2

%% define pressure nodes
np = n*(2*n-1)-1;
p = zeros(np,5);

for i = 1:n % loop over x-axis
    for j = 1:(2*n-1) % loop over y-axis
        
        ip = (j-1)*n+i; % point identification number
        
        % x location
        if mod(j,2)
            % if j is odd
            p(ip,1) = gridx1(i);
        else
            p(ip,1) = gridx2(i);
        end
        
        % y location
        p(ip,2) = gridy(j);
        
    end
end
% create outflow point above top right corner (if wanted)
% p(np+1,:) = [gridx2(end) gridy(end)+Ledge*sqrt(3)./2 0 0 0]
% plot to check network nodes
% figure
% plot(p(:,1),p(:,2),'*')
% axis equal

%% define connections
c = zeros(1,7);

% x connections
ic = 0;
for i = 1:n
    for j = 1:(2*n-1)
        
        if mod(j,2)
            % if j is odd
            if mod(i,2) == 0
                % if i is even
                ic = ic + 1; % connection identification number
                % left node
                c(ic,1) = (j-1)*n+i; % assign p0
                % right node
                c(ic,2) = (j-1)*n+i+1; % assign p1
                connected = 1;
            end
        else
            % if j is even
            if mod(i,2)
                % if i is odd
                ic = ic + 1; % connection identification number
                % left node
                c(ic,1) = (j-1)*n+i; % assign p0
                % right node
                c(ic,2) = (j-1)*n+i+1; % assign p1
                connected = 1;
            end
        end
    
        % update pipe-per-node count
        if ic && connected
            if c(ic,2) > np
                break
            end
            p(c(ic,1),4) = p(c(ic,1),4) + 1; % p0 has now one more node
            p(c(ic,2),4) = p(c(ic,2),4) + 1; % p1 has now one more node
        end
        connected = 0;
    end
end

% y connections
for i = 1:n
    for j = 1:(2*n-2)
        
        ic = ic + 1; % connection identification number
        % lower node
        c(ic,1) = (j-1)*n+i; % assign p0
        % upper node
        c(ic,2) = j*n+i; % assign p1
        
        if ic
            if c(ic,2) > np
                break
            end
            % update pipe-per-node count
            p(c(ic,1),4) = p(c(ic,1),4) + 1; % p0 has now one more node
            p(c(ic,2),4) = p(c(ic,2),4) + 1; % p1 has now one more node
        end
    end
end

% length of connections
c(:,6) = sqrt((p(c(:,1),1)-p(c(:,2),1)).^2 ...
    + (p(c(:,1),2)-p(c(:,2),2)).^2);

% pipe radii
c(:,8) = 2.5e-6;

% get rid of bottom left point (dead end)
% find node at the bottom left corner
to_kill = find(p(:,1) == min(p(:,1)) & p(:,2) == max(p(:,2)));
p(to_kill,:) = [];
c_to_kill = find(c(:,1) == to_kill); % search through all p0
c_to_kill2 = find(c(:,2) == to_kill); % search through all p1
c(c_to_kill,:) = []; % delete the connection entry
c(c_to_kill2,:) = []; % delete the connection entry
% reduce the count after the node that has been killed by one
for i = 1:size(c,1)
    if c(i,1) >= to_kill
        c(i,1) = c(i,1) - 1;
    end
    if c(i,2) >= to_kill
        c(i,2) = c(i,2) - 1;
    end
end
