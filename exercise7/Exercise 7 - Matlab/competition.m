close all
clear
clc
tic
if  not(isfolder('Figures'))
       mkdir('Figures')
end
%% Ask if you want to plot the pictures
selection = questdlg('Do you want to coninuosly plot the RBC movement? (takes longer to compute)','Figure',...
                        'Yes','No', 'No');
switch selection
    case 'Yes'
        imagesplot = 1;
    case 'No'
        imagesplot = 0;
end
%% Set simulation parameters
t = 0;          % simulation start time (dont change)
tout = 0;
dtout = 0.025;  % time interval (s) between two outputs
tend = 20;      % simulation stops at t=tend seconds

plotmode=7;     % hexagonal network
fontsize1 = 12; % define plot font size

% RBC volume in rat 55fl
Vrbc = 55e-18;

% Size of network (do NOT change)
N=11; % must be an odd number (hexagonal network)

% In vivo or in vitro calculation of the effects
invivo=0.0;

% Edge length
Ledge = 15e-06;

%% test Name
b= 'Simulation name';
a=1;
while a==1
    
name= inputdlg({b}, 'Information',[1 50]);
waitfor(name)
Path=cd;
if  contains(Path, '/')
    Folder=[Path '/Figures/' ];
else
    Folder=[Path '\Figures\'];
end

    if isfile([Folder char(name) '_Setting.mat']) % check if there is already a test with the same name
          answer = questdlg('A simulation with this name alredy exist. Do you want to overwrite the simulation?', ...
            'Overwrite test','Yes','No','No');
            waitfor(answer)
            % Handle response
            switch answer
                case 'Yes'
                    disp('The simulation will be overwritten.')
                    a=0;
                case 'No'
                    disp('Choose new name')
                    b= 'New simulation name';
                    a=1;
            end
    else 
        a=0;
    end
end


%% Define nodes and connections
disp('...calling Matlab function create_network.m');
[p,c] = create_network(N,Ledge);

%  vertices (points)
%       1      2   3  4       5     6      7
%       x      y   P  nc      f     Htavg  favg
% p = [ ...

% connections (between two points p0 and p1)
%        1   2   3   4     5       6        7   8  9   10
%        p0  p1  T   nrbc  uBulk   Ledge    HT  a  HD  vf
% c = [ ...

% red blood cells
%         1   2   3               4   5      6
%         x   y   connectionID    s   Tage   Tremain
% r = [ ...




% Number of nodes and connections
np = size(p,1);
nc = size(c,1);

%constriction = [];
%dilation = [];

% Connection distance of nodes from inlet
i=1;
gi=0;
gp=repmat(-1,np,1);
while ~isempty(i)
    gp(i)=gi;
    gi=gi+1;
    ic=any(ismember(c(:,1:2),i),2);
    i=setdiff(unique(c(ic,1:2)),find(gp>=0));
end
gc=min(gp(c(:,1:2)),[],2);
% indp=find(gp<=18);
% valc=find(any(ismember(c(:,1:2),indp),2));

% Select downstream vessels
xm=max(p(:,1));
ym=max(p(:,2));
l=cross([xm 0 1],[0 ym 1]);
n=l(1)*p(:,1)+l(2)*p(:,2)+l(3);
indp=find(n>=0);
valc=find(any(ismember(c(:,1:2),indp),2))';


%% Select vessels to change the diameter
% #########################################################################


% 1.) You can pick the capillaries manually by a GUI:

[constriction,dilation, done]=select_capillaries(c,p,valc,name,Folder);


% 2.) You can also select vessels directly by addressing their IDs:

% constriction = [205 206 208 209] ;
% dilation = [225 226 227];
% done=1; %Do NOT change => See select_capillaries function
% save([Folder char(name) '_Setting.mat'],'conTeststriction','dilation')


% 3.) You can pick the capillaries manually by a GUI starting with an existing list of IDs:

% constriction = [205 206 208 209] ;
% dilation = [225 226 227];
% [constriction,dilation,done]=select_capillaries(c,p,valc,name,Folder,constriction,dilation);


% #########################################################################
%% if not saved

if done==1  
else
    disp('Capillary network closed before it could save the user settings.')  
    return
end
%% Figures
% Get screen size
scrsz=get(0,'ScreenSize');

% Set up figure positions
f1 = figure(1);
set(f1,'OuterPosition',[20 50 (scrsz(3)/2-40) (2*(scrsz(4)/3))]);
f2 = figure(2);
set(f2,'OuterPosition',[(scrsz(3)/2+40) 50 (scrsz(3)/2-40) (2*(scrsz(4)/3))]);

%% Update diameters
% ensure that all selected vessels are valid options
dilation=intersect(dilation,valc);
constriction=intersect(constriction,valc);

% The initial vessel diameter (5 microns) is changed by 40% 
if constriction
    % Make them 40% smaller
    c(constriction,8) = c(constriction,8)*0.6;
end

if dilation 
    % Make them 40% larger
    c(dilation,8) = c(dilation,8)*1.4; 
end




%% RBC setup

% Law for RBC bifurcation
bifurcationlaw = 1; % RBCs follow the largest bulk flow velocity

% Markers
markerstart=0;
markerduration=1100;

% Average hematocrit in the network
H = 0.25;
% Total number of RBCs
nr = floor((H*nc*Ledge*pi*mean(c(:,8))^2)/Vrbc);
% Upper limit for hematocrit values
Htmax=1;

randcon=zeros(1,nr);
randleng=zeros(1,nr);

%define the random numbers but same if you restart the simulation with the
%same info Change to loaded random position.
% for i=1:nr
%       %rng(i);
%       randcon(i)=rand;
%       %rng(nr+1);
%       randleng(i)=rand;
%     
% end
load('randcon')
load('randleng')

% Loop through all RBCs
for i = 1:nr
    % Randomly assign connection IDs for all RBCs
  
    r(i,3) = floor(randcon(i)*nc)+1;

    % Increment the RBC count on the connection vessel
    c(r(i,3),4) = c(r(i,3),4) + 1;

    % Get the x-y coords of that vessel's start vertex
    x1 = p(c(r(i,3),1),1);
    y1 = p(c(r(i,3),1),2);
    
    % Get the x-y coords of that vessel's end vertex
    x2 = p(c(r(i,3),2),1);
    y2 = p(c(r(i,3),2),2);
    
    % Get the vessel length
    L = c(r(i,3),6);

    % Create a random location for the RBC along the vessel
    r(i,4) = randleng(i)*L;
    
    % Store the exact coordinates of the RBC
    r(i,1) = x1 + r(i,4)/L*(x2-x1);
    r(i,2) = y1 + r(i,4)/L*(y2-y1);

    % Assign an age to the RBC
    r(i,5) = -1;
end


%% Compute tube and discharge Hematocrit

% Update the tube hematocrit
c(:,7) = c(:,4)*Vrbc./(c(:,6).*pi.*c(:,8).^2);

% Update the discharge hematocrit
c(:,9) = compute_discharge_from_tube(c(:,7),2*c(:,8)*1e6,invivo);

% Initialize velocity correction
vf = c(:,9)./c(:,7);
for i=1:length(vf)
    if (vf(i) < 1)
        vf(i) = 1;
    end
end
c(:,10) = vf;

%% Compute network

disp('...calling Matlab function compute_network');
[p,c] = compute_network(p,c,0,invivo);


%% Initialize average hematocrit and flow velocity at nodes

p(:,6) = 0;
p(:,7) = 0;


%% Initialize history

Thistory = 0;            % store the current time
Hhistory = mean(c(:,7)); % store the current tube hematocrit ensemble average
Dhistory = std(c(:,7));  % store the current tube hematocrit standard deviation


%% Compute remaining time to next bifurcation
disp('...calling Matlab function bifurcation_time');
r(:,6) = bifurcation_time(r,p,c) + t;

% Time at next bifurcation
tnext = min(r(:,6));
irnext = find(r(:,6)==tnext,1);
icnext = r(irnext,3);


%% Integration in time
passagetime = zeros(100,2);

% Plot the network
disp('...calling Matlab function plot_network');
figure(1)
plot_network(r,p,c);

% Inform about user changes on the network
disp([sprintf('\nThe largest vessel has a diameter of ') num2str(max(c(:,8)*2)*1e6) ' microns.']); 
disp(['The smallest vessel has a diameter of ' num2str(min(c(:,8)*2)*1e6) ' microns.']);
numOfChanges = 0;
for vnr = 1:nc
    if (not (c(vnr,8) == 2.5e-06))
        numOfChanges = numOfChanges +1;
    end
end
disp(['There was a total of ' num2str(numOfChanges) sprintf(' changed vessel diameters.\n')]);

% Create the goal bin
binDelta = 0.2;
binStart = H-binDelta/2;
binEnd   = H+binDelta/2;
binCount = 0;
binCountAverage = 0;
avCount = 0;

% Count the hits in the bin
for i=1:nc
    if ( (c(i,7)>=binStart) && (c(i,7)<=binEnd) )
        binCount = binCount + 1;
    end
end

% Initialize a binCount history
binCountHistoryT = 0;
binCountHistoryH = binCount;

% Reset the binCount to zero
binCount = 0;
    
% Evolve in time
while (t < tend)

    dt = min(tout,tnext)-t;
    dt = max(1e-6,dt);

    %pause(1)
    t = t + dt;
    
    % Propagate all RBCs
    ic = r(:,3);

    % Calculate all RBC's vessel vertex locations
    x1 = p(c(ic,1),1);
    y1 = p(c(ic,1),2);
    x2 = p(c(ic,2),1);
    y2 = p(c(ic,2),2);
    
    % Calculate all RBC's vessel lengths
    L = c(ic,6);
    
    % Calculate all RBC's velocities
    v = c(ic,5);

    % Calculate velocity factor of fahraeus effect (HD/HT)
    vf = c(ic,9)./c(ic,7);
    for i=1:length(vf)
        if (vf(i) < 1)
            vf(i) = 1;
        end
    end
    c(ic,10) = vf;
    
    % Correct the RBC velocities with the Fahraeus effect
    v = v.*c(ic,10);

    % Get the next vessel location of the RBC,
    % but restrict its position to the length of the vessel
    r(:,4) = min(L,max(0,r(:,4) + dt*v));
    
    % Recalculate the exact RBC coordinates
    r(:,1) = x1 + r(:,4)./L.*(x2-x1);
    r(:,2) = y1 + r(:,4)./L.*(y2-y1);

    % Check for bifurcation incident
    if (t>=tnext)

        % Identify bifurcation node
        if (v(irnext)>0)
            bifurcationnode = c(icnext,2);
        else
            bifurcationnode = c(icnext,1);
        end

        % Exit condition (periodic)
        if (bifurcationnode==np)
            bifurcationnode = 1;

            if (r(irnext,5)>0)
                
                pt = round(t-r(irnext,5));
                pt = max(1,pt);
                
                % Red RBC exiting
                passagetime(pt,2) = passagetime(pt,2) + 1;
            end

            % Set new start time for RBC
            pt = min(floor(t+1),100);
            if (t>=markerstart && t<=markerstart+markerduration)
                
                % New blue RBC
                r(irnext,5) = t;
                passagetime(pt,1) = passagetime(pt,1) + 1;
            else
                r(irnext,5) = -1;
            end

        end

        % Bifurcation rule
        switch bifurcationlaw
            
            case 1
                %% RBC prefer largest bulk flow velocity
                vBulk = 0;
                newpipe = icnext;
                for j = 1:nc
                    % find largest bulk flow velocity at bifurcaiton node
                    if (c(j,1)==bifurcationnode)
                        vBulkPipe = -c(j,5); % Q in the node positive, in out negative
                        Ht = c(j,7);
                        if (vBulkPipe < vBulk && Ht < Htmax)
                            newpipe = j;
                            vBulk = vBulkPipe;
                            r(irnext,1) = p(c(newpipe,1),1);
                            r(irnext,2) = p(c(newpipe,1),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = 0;
                        end
                    end

                    if (c(j,2)==bifurcationnode)
                        vBulkPipe = c(j,5);
                        Ht = c(j,7);
                        %                if (rand<0.5)
                        if (vBulkPipe < vBulk && Ht < Htmax)
                            newpipe = j;
                            vBulk = vBulkPipe;
                            r(irnext,1) = p(c(newpipe,2),1);
                            r(irnext,2) = p(c(newpipe,2),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = c(newpipe,6);
                        end
                    end

                end

            case 2
                %% RBC prefer largest pressure gradient
                dpdxmin = 0;
                newpipe = icnext;
                for j = 1:nc

                    % find largest pressure gradient at bifurcaiton node
                    if (c(j,1)==bifurcationnode)
                        dpdx = -c(j,5)/c(j,8)^4 + (rand-0.5)*eps;
                        Ht = c(j,7);
                        if (dpdx<dpdxmin && Ht < Htmax)
                            %                if (rand<0.5)
                            newpipe = j;
                            dpdxmin = dpdx;
                            r(irnext,1) = p(c(newpipe,1),1);
                            r(irnext,2) = p(c(newpipe,1),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = 0;
                        end
                    end

                    if (c(j,2)==bifurcationnode)
                        dpdx = c(j,5)/c(j,8)^4 + (rand-0.5)*eps;
                        Ht = c(j,7);
                        %                if (rand<0.5)
                        if (dpdx<dpdxmin && Ht < Htmax)
                            newpipe = j;
                            dpdxmin = dpdx;
                            r(irnext,1) = p(c(newpipe,2),1);
                            r(irnext,2) = p(c(newpipe,2),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = c(newpipe,6);
                        end
                    end

                end

               
            case 3
                %% RBC bifurcate randomly into distal vessels
                dpdxmin = 0;
                newpipe = icnext;
                for j = 1:nc

                    % find largest pressure gradient at bifurcaiton node
                    if (c(j,1)==bifurcationnode)
                        dpdx = -rand*sign(c(j,5));
                        Ht = c(j,7);
                        if (dpdx<dpdxmin && Ht < Htmax)
                            %                if (rand<0.5)
                            newpipe = j;
                            dpdxmin = dpdx;
                            r(irnext,1) = p(c(newpipe,1),1);
                            r(irnext,2) = p(c(newpipe,1),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = 0;
                        end
                    end

                    if (c(j,2)==bifurcationnode)
                        dpdx = rand*sign(c(j,5));
                        Ht = c(j,7);
                        %                if (rand<0.5)
                        if (dpdx<dpdxmin && Ht < Htmax)
                            newpipe = j;
                            dpdxmin = dpdx;
                            r(irnext,1) = p(c(newpipe,2),1);
                            r(irnext,2) = p(c(newpipe,2),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = c(newpipe,6);
                        end
                    end

                end

               
            case 4
                %% RBC follow plasma flow as passive tracers

                rvalue = rand * p(bifurcationnode,5);

                for j = 1:nc

                    if (c(j,1)==bifurcationnode && c(j,5)>0)
                        rvalue = rvalue - c(j,5);
                        if (rvalue<eps)
                            rvalue = p(bifurcationnode,5);
                            newpipe = j;
                            r(irnext,1) = p(c(newpipe,1),1);
                            r(irnext,2) = p(c(newpipe,1),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = 0.01*c(newpipe,6);
                        end
                    end

                    if (c(j,2)==bifurcationnode && c(j,5)<0)
                        rvalue = rvalue + c(j,5);
                        if (rvalue<eps)
                            rvalue = p(bifurcationnode,5);
                            newpipe = j;
                            r(irnext,1) = p(c(newpipe,2),1);
                            r(irnext,2) = p(c(newpipe,2),2);
                            r(irnext,3) = newpipe;
                            r(irnext,4) = 0.99*c(newpipe,6);
                        end
                    end

                end

        end


        if (newpipe == icnext)
            
            % No new pipe found for RBC
            % --> re-set RBC to previous position
            r(irnext,4) = r(irnext,4) - 0.01*sign(c(icnext,5))*c(icnext,6);

        else
            
            % Book keeping for number of RBC in pipe
            c(icnext,4) = c(icnext,4) - 1;
            c(newpipe,4) = c(newpipe,4) + 1;

            % Update tube and discharge Hematocrit
            c(icnext,7) = c(icnext,4)*Vrbc./(c(icnext,6).*pi.*c(icnext,8).^2);
            c(newpipe,7) = c(newpipe,4)*Vrbc./(c(newpipe,6).*pi.*c(newpipe,8).^2);
            c(icnext,9) = compute_discharge_from_tube(c(icnext,7),2*c(icnext,8)*1e6,invivo);
            c(newpipe,9) = compute_discharge_from_tube(c(newpipe,7),2*c(newpipe,8)*1e6,invivo);
            
            % Recompute network
            %disp('...calling Matlab function compute_network');
            [p,c] = compute_network(p,c,t,invivo);

        end

        % Update bifurcation time
        r(:,6) = bifurcation_time(r,p,c) + t;

        % Time at next bifurcation
        tnext = min(r(:,6));
        irnext = find(r(:,6)==tnext,1);
        icnext = r(irnext,3);

    else

        tout = tout + dtout;
        %disp(['t = ' num2str(tout) 's']);

        % Compute time averaged hematocrit at nodes
        for i = 1:nc
           p(c(i,1),6) = p(c(i,1),6) + dtout*c(i,7)/p(c(i,1),4);
           p(c(i,2),6) = p(c(i,2),6) + dtout*c(i,7)/p(c(i,2),4);
        end

        % Compute average flow rate at nodes
        p(:,7) = p(:,7)+dtout*p(:,5);

        % Store the current time
        Thistory = [Thistory t];
        % Store the current tube hematocrit ensemble average
        Hhistory = [Hhistory mean(c(:,7))];
        % Store the current tube hematocrit standard deviation
        Dhistory = [Dhistory std(c(:,7))];
        
        % Count the hits in the bin
        for i=1:nc
            if ( (c(i,7)>=binStart) && (c(i,7)<=binEnd) )
                binCount = binCount + 1;
            end
        end

        % Initialize a binCount history
        binCountHistoryT = [binCountHistoryT t];
        binCountHistoryH = [binCountHistoryH binCount];

        % Reset the binCount to zero
        binCount = 0;
   if (t>10)
        avCount = avCount +1;
             binCountAverage = (binCountAverage*(avCount-1) ...
                +binCountHistoryH(end))/avCount;
   end
        %%  
        
  if imagesplot      
        % Plot the network
        figure(1);
        plot_network(r,p,c);
        
        f2=figure(2);
        figure(2);
        plot(binCountHistoryT,binCountHistoryH);
        xlim([0 20]);
        ylim([0 nc]);
        xlabel('time t [s]', 'FontSize',fontsize1);
        ylabel({'number of vessels with a hematocrit'; ['between ' ...
            num2str(H-binDelta/2) ' and ' num2str(H+binDelta/2)]}, ...
            'FontSize',fontsize1);
        
        % Start averaging the binCount after 10s
        if (t<=10)
            text(2, nc*0.8, 'reaching steady state', ...
                'FontSize',fontsize1);
        else
            avCount = avCount +1;
            binCountAverage = (binCountAverage*(avCount-1) ...
                +binCountHistoryH(end))/avCount;
            text(12, nc*0.8, 'temporal averaging active', ...
                'FontSize',fontsize1);
            text(12, nc*0.7, ['current bin count: ' ...
                num2str(binCountHistoryH(end))], 'FontSize',fontsize1);
            text(12, nc*0.6, ['average bin count: ' ...
                num2str(binCountAverage)], 'FontSize',fontsize1);
        end
  end

    end


end



if ~isfolder('Figures')
    mkdir('Figures')
end

figure(1);
f1=figure(1);
plot_network(r,p,c);
set(f1,'OuterPosition',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.8 scrsz(4)*0.8]);
for j=1:nc
    xvesselcenter = p(c(j,1),1) + (p(c(j,2),1) - p(c(j,1),1))/2;
    yvesselcenter = p(c(j,1),2) + (p(c(j,2),2) - p(c(j,1),2))/2;
    if (2*c(j,8)==5e-06)
        text(xvesselcenter,yvesselcenter,num2str(j),...
            'HorizontalAlignment','center',...
            'Color','black','BackgroundColor',[.9 .9 .9]);
    elseif (2*c(j,8)>5e-06)
        text(xvesselcenter,yvesselcenter,num2str(j),...
            'HorizontalAlignment','center',...
            'Color','black','BackgroundColor','green');
    else
        text(xvesselcenter,yvesselcenter,num2str(j),...
            'HorizontalAlignment','center',...
            'Color','black','BackgroundColor','red');
    end
end
%saveas(f1,[Folder char(name) '_Capillary_Network.png']);

figure(2);
 plot(binCountHistoryT,binCountHistoryH);
        xlim([0 20]);
        ylim([0 nc]);
        xlabel('time t [s]', 'FontSize',fontsize1);
        ylabel({'number of vessels with a hematocrit'; ['between ' ...
            num2str(H-binDelta/2) ' and ' num2str(H+binDelta/2)]}, ...
            'FontSize',fontsize1);
            text(12, nc*0.8, 'temporal averaging active', ...
                'FontSize',fontsize1);
            text(12, nc*0.7, ['current bin count: ' ...
                num2str(binCountHistoryH(end))], 'FontSize',fontsize1);
            text(12, nc*0.6, ['average bin count: ' ...
                num2str(binCountAverage)], 'FontSize',fontsize1);
saveas(f2,[ Folder char(name) '_Bin_Count.png']);


% Finalize the averaging of flow rates and hematocrit
p(:,6:7) = p(:,6:7)/t;


%% Create a plot that shows the time-averaged hematocrit
fprintf('\n');
disp('...plot the time-averaged hematocrit');

[phandle1] = plot_time_averaged_flow(c, p, 'Htavg');

% Display the constricted (red) and dilated (green) vessels
if (not (constriction==0)) 
    plot([p(c(constriction,1),1),p(c(constriction,2),1)]',[p(c(constriction,1),2) p(c(constriction,2),2)]','r-','linewidth',3);
end
if (not (dilation==0))
    plot([p(c(dilation,1),1),p(c(dilation,2),1)]',[p(c(dilation,1),2) p(c(dilation,2),2)]','g-','linewidth',3);
end
saveas(gcf,[Folder char(name) '_Hematocrit.png']);

%% Create a plot that shows the time-averaged flow
disp('...plot the time-averaged flow');

[phandle2] = plot_time_averaged_flow(c, p, 'favg');

% Display the constricted (red) and dilated (green) vessels
if (not (constriction==0)) 
    plot([p(c(constriction,1),1),p(c(constriction,2),1)]',[p(c(constriction,1),2) p(c(constriction,2),2)]','r-','linewidth',3);
end
if (not (dilation==0))
    plot([p(c(dilation,1),1),p(c(dilation,2),1)]',[p(c(dilation,1),2) p(c(dilation,2),2)]','g-','linewidth',3);
end
saveas(gcf,[Folder char(name) '_Flow.png']);

%% Plot the time history of the desired vessel
disp('...plot the history of the desired vessel');
f5 = figure(5);
set(f5,'OuterPosition',[1 scrsz(4)/2+1 scrsz(3)/3 scrsz(4)/2-20]);
plot(Thistory,Hhistory);
hold on
plot(Thistory,Hhistory+Dhistory,'r');
plot(Thistory,Hhistory-Dhistory,'r');
xlim([0 20]);
ylim([-0.2 1.2]);
xlabel('time [s]', 'FontSize',fontsize1);
ylabel('average tube hematocrit [-]', 'FontSize',fontsize1);
saveas(f5,[Folder char(name) '_Standard_Deviation.png']);

Percentage_bincount=binCountAverage/nc*100;
save([Folder char(name) '_Results.mat'],'Percentage_bincount')

toc

fprintf('\n');
disp(['Congrats! You achieved a desirable tube hematocrit in ' num2str(binCountAverage/nc*100) '% of your vessels!']);
fprintf('\n');
