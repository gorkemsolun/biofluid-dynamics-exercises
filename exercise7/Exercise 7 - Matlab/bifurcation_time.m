function [Tremain] = bifurcation_time(r,p,c)

nr = size(r,1);
nc = size(c,1);
np = size(p,1);

%% Compute time to next bifurcation

for i = 1:nr
    
    % number of pipe
    ic = r(i,3);

    % position in pipe
    s = r(i,4);
    
    % pipe length
    l = c(ic,6);

    % flow velocity in pipe
    v = c(ic,5);
    
    % Fahraeus effect
    v = v.*c(ic,10); % Correction

    % remaining time in pipe
    if (v>0)
        Tremain(i,1) = (l-s)/v;
    else
        Tremain(i,1) = -s/v;
    end

end

end