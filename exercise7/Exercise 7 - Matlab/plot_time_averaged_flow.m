function [phandle] = plot_time_averaged_flow(c, p, varargin)
%  contour plots for mean flow and Ht 
% INPUT: c, p: as constructed by fork.m
% OUTPUT: phandle: Plot handle.
%varargin: 'favg' for plotting the flows or 'Htavg' for the Ht
    figure
    r = p(:,1:2);
    el = c(:,1:2);
    hold on;
    
    N = ceil(size(r,1)/1.8);
    l = max(c(:,6));
    xmax = max(r(:,1));
    ymax = max(r(:,2));
    l = min(c(:,6));
    if length(varargin) > 0
        [XI,YI] = meshgrid([0:xmax/N:xmax], [0:(ymax-l)/N:ymax]);
        x = r(:,1);
        y = r(:,2);
        if strcmp(varargin{1}, 'Htavg')
            z = p(:,6);
        elseif strcmp(varargin{1}, 'favg')
            z = p(:,7);
        end
        F = TriScatteredInterp(x,y,z);
        ZI = F(XI,YI);
        %ZI = griddata(x,y,z,XI,YI);
        ZI(isnan(ZI)) = 0;
        contourf(XI, YI, ZI);
        colorbar()
        if strcmp(varargin{1}, 'Htavg')
            title('Hematocrit');
        elseif strcmp(varargin{1}, 'favg')
            title('Flow');
        end
    end
	
	
	scatter(r(:, 1), r(:, 2));
	axis('equal');
	% Plot capillaries:
	ceids = 1:size(c,1);
	for j=1:length(ceids)
		i = ceids(j);
		line([r(el(i,1),1), r(el(i,2),1)],...
		     [r(el(i,1),2), r(el(i,2),2)], 'Color', 'black');
    end
	set(gca,'Ydir','reverse','xlim', [-l, xmax+l], 'ylim', [-l, ymax+l]);
    phandle = gca;
    
end


