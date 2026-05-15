function [Percentage_bincount, valc, c, p, HtAvg] = competition_headless(constriction, dilation)
% Headless version of competition.m: no GUI, no plots, no saving.
% Inputs: constriction, dilation -- vectors of capillary IDs (top-left half).
% Outputs:
%   Percentage_bincount -- percentage of vessels in HT = 0.25 +- 0.1
%   valc -- list of valid (modifiable) capillary IDs
%   c    -- final capillary matrix
%   p    -- final pressure-node matrix
%   HtAvg -- time-averaged tube hematocrit per vessel (10s..20s)

%% Simulation parameters
t = 0;
tout = 0;
dtout = 0.025;
tend = 20;

Vrbc = 55e-18;
N = 11;
invivo = 0.0;
Ledge = 15e-06;

%% Build network
[p, c] = create_network(N, Ledge);

np = size(p,1);
nc = size(c,1);

% Connection distance from inlet (not strictly needed for evaluation, but keep parity)
i = 1;
gi = 0;
gp = repmat(-1, np, 1);
while ~isempty(i)
    gp(i) = gi;
    gi = gi + 1;
    ic = any(ismember(c(:,1:2), i), 2);
    i = setdiff(unique(c(ic,1:2)), find(gp >= 0));
end

% Select downstream / top-left-half vessels (variable capillaries)
xm = max(p(:,1));
ym = max(p(:,2));
l = cross([xm 0 1], [0 ym 1]);
n = l(1)*p(:,1) + l(2)*p(:,2) + l(3);
indp = find(n >= 0);
valc = find(any(ismember(c(:,1:2), indp), 2))';

%% Apply user changes (intersect with valc to enforce the rule)
dilation = intersect(dilation, valc);
constriction = intersect(constriction, valc);

if ~isempty(constriction)
    c(constriction, 8) = c(constriction, 8) * 0.6;
end
if ~isempty(dilation)
    c(dilation, 8) = c(dilation, 8) * 1.4;
end

%% RBC setup
bifurcationlaw = 1;

H = 0.25;
nr = floor((H * nc * Ledge * pi * mean(c(:,8))^2) / Vrbc);
Htmax = 1;

S = load('randcon');
randcon = S.randcon;
S = load('randleng');
randleng = S.randleng;

r = zeros(nr, 6);
for i = 1:nr
    r(i,3) = floor(randcon(i) * nc) + 1;
    c(r(i,3), 4) = c(r(i,3), 4) + 1;
    x1 = p(c(r(i,3),1), 1);
    y1 = p(c(r(i,3),1), 2);
    x2 = p(c(r(i,3),2), 1);
    y2 = p(c(r(i,3),2), 2);
    L = c(r(i,3), 6);
    r(i,4) = randleng(i) * L;
    r(i,1) = x1 + r(i,4)/L*(x2-x1);
    r(i,2) = y1 + r(i,4)/L*(y2-y1);
    r(i,5) = -1;
end

%% Initial tube & discharge hematocrit
c(:,7) = c(:,4) * Vrbc ./ (c(:,6) .* pi .* c(:,8).^2);
c(:,9) = compute_discharge_from_tube(c(:,7), 2*c(:,8)*1e6, invivo);

vf = c(:,9) ./ c(:,7);
vf(vf < 1) = 1;
c(:,10) = vf;

%% Solve network
[p, c] = compute_network(p, c, 0, invivo);

p(:,6) = 0;
p(:,7) = 0;

%% Time to next bifurcation
r(:,6) = bifurcation_time(r, p, c) + t;
tnext = min(r(:,6));
irnext = find(r(:,6) == tnext, 1);
icnext = r(irnext, 3);

%% Bin counting setup
binDelta = 0.2;
binStart = H - binDelta/2;
binEnd   = H + binDelta/2;
binCount = 0;
binCountAverage = 0;
avCount = 0;
binCountHistoryH = 0;

passagetime = zeros(100, 2);
markerstart = 0;
markerduration = 1100;

% Time-averaged tube hematocrit per vessel (window: t > 10s)
HtAcc = zeros(nc, 1);
HtAccCount = 0;

%% Main loop
while (t < tend)
    dt = min(tout, tnext) - t;
    dt = max(1e-6, dt);
    t = t + dt;

    ic = r(:,3);
    x1 = p(c(ic,1), 1);
    y1 = p(c(ic,1), 2);
    x2 = p(c(ic,2), 1);
    y2 = p(c(ic,2), 2);
    L = c(ic, 6);
    v = c(ic, 5);

    vf = c(ic, 9) ./ c(ic, 7);
    vf(vf < 1) = 1;
    c(ic, 10) = vf;
    v = v .* c(ic, 10);

    r(:,4) = min(L, max(0, r(:,4) + dt * v));
    r(:,1) = x1 + r(:,4) ./ L .* (x2 - x1);
    r(:,2) = y1 + r(:,4) ./ L .* (y2 - y1);

    if (t >= tnext)
        if (v(irnext) > 0)
            bifurcationnode = c(icnext, 2);
        else
            bifurcationnode = c(icnext, 1);
        end

        if (bifurcationnode == np)
            bifurcationnode = 1;
            if (r(irnext, 5) > 0)
                pt = round(t - r(irnext, 5));
                pt = max(1, pt);
                passagetime(pt, 2) = passagetime(pt, 2) + 1;
            end
            pt = min(floor(t + 1), 100);
            if (t >= markerstart && t <= markerstart + markerduration)
                r(irnext, 5) = t;
                passagetime(pt, 1) = passagetime(pt, 1) + 1;
            else
                r(irnext, 5) = -1;
            end
        end

        % Only case 1 (largest bulk flow velocity) is needed
        vBulk = 0;
        newpipe = icnext;
        for j = 1:nc
            if (c(j,1) == bifurcationnode)
                vBulkPipe = -c(j,5);
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
            if (c(j,2) == bifurcationnode)
                vBulkPipe = c(j,5);
                Ht = c(j,7);
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

        if (newpipe == icnext)
            r(irnext,4) = r(irnext,4) - 0.01 * sign(c(icnext,5)) * c(icnext,6);
        else
            c(icnext,4) = c(icnext,4) - 1;
            c(newpipe,4) = c(newpipe,4) + 1;
            c(icnext,7) = c(icnext,4) * Vrbc ./ (c(icnext,6) .* pi .* c(icnext,8).^2);
            c(newpipe,7) = c(newpipe,4) * Vrbc ./ (c(newpipe,6) .* pi .* c(newpipe,8).^2);
            c(icnext,9) = compute_discharge_from_tube(c(icnext,7), 2*c(icnext,8)*1e6, invivo);
            c(newpipe,9) = compute_discharge_from_tube(c(newpipe,7), 2*c(newpipe,8)*1e6, invivo);
            [p, c] = compute_network(p, c, t, invivo);
        end

        r(:,6) = bifurcation_time(r, p, c) + t;
        tnext = min(r(:,6));
        irnext = find(r(:,6) == tnext, 1);
        icnext = r(irnext, 3);

    else
        tout = tout + dtout;
        for i = 1:nc
            p(c(i,1),6) = p(c(i,1),6) + dtout*c(i,7)/p(c(i,1),4);
            p(c(i,2),6) = p(c(i,2),6) + dtout*c(i,7)/p(c(i,2),4);
        end
        p(:,7) = p(:,7) + dtout * p(:,5);

        binCount = 0;
        for i = 1:nc
            if (c(i,7) >= binStart) && (c(i,7) <= binEnd)
                binCount = binCount + 1;
            end
        end
        binCountHistoryH = binCount;

        if (t > 10)
            avCount = avCount + 1;
            binCountAverage = (binCountAverage*(avCount-1) + binCountHistoryH) / avCount;
            HtAcc = HtAcc + c(:,7);
            HtAccCount = HtAccCount + 1;
        end
    end
end

Percentage_bincount = binCountAverage / nc * 100;
if HtAccCount > 0
    HtAvg = HtAcc / HtAccCount;
else
    HtAvg = c(:,7);
end
end
