function run_strategy_v2()
% Refine around combo_K10 = 17.02%
D = load('diag.mat');
valc = D.valc(:)';
HtAvg = D.HtAvg;
ubulk = D.ubulk;

htv = HtAvg(valc);
[~, ord_hi] = sort(htv, 'descend');
valc_hi = valc(ord_hi);

ub = abs(ubulk(valc));
htl = htv;
htl(ub < 1e-8) = NaN;
[~, ord_lo] = sort(htl, 'ascend');
valc_lo_all = valc(ord_lo);
% Drop those that ended up NaN
mask = ~isnan(htl(ord_lo));
valc_lo = valc_lo_all(mask);

fid = fopen('headless_log.txt', 'a');
fprintf(fid, '--- v2 ---\n');
trials = {};

% Asymmetric K (constrict, dilate)
for kc = [8 9 10 11 12]
    for kd = [6 8 10 12 14]
        trials{end+1} = struct('name', sprintf('asym_c%d_d%d',kc,kd), ...
            'constriction', valc_hi(1:kc), 'dilation', valc_lo(1:kd));
    end
end

best_pct = -inf; best_name = ''; best_con = []; best_dil = [];

for k = 1:length(trials)
    tr = trials{k};
    t0 = tic;
    pct = competition_headless(tr.constriction, tr.dilation);
    el = toc(t0);
    line = sprintf('%-15s nc=%d nd=%d  pct=%.4f  (%.1fs)', tr.name, ...
        numel(tr.constriction), numel(tr.dilation), pct, el);
    fprintf('%s\n', line);
    fprintf(fid, '%s\n', line);
    if pct > best_pct
        best_pct = pct; best_name = tr.name;
        best_con = tr.constriction; best_dil = tr.dilation;
    end
end

fprintf('BEST v2: %s pct=%.4f\n', best_name, best_pct);
fprintf(fid, 'BEST v2: %s pct=%.4f\n', best_name, best_pct);
fclose(fid);
save('best_v2.mat', 'best_pct', 'best_name', 'best_con', 'best_dil');
end
