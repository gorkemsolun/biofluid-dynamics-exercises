function run_strategy_v1()
% Try several constrict-high / dilate-low strategies in one MATLAB session.

D = load('diag.mat');
valc = D.valc(:)';
HtAvg = D.HtAvg;
ubulk = D.ubulk;
nc = length(HtAvg);

% Order modifiable vessels by Ht (descending = high)
htv = HtAvg(valc);
[hi, ord_hi] = sort(htv, 'descend');
valc_hi = valc(ord_hi);

% Order modifiable vessels by Ht ascending (low) but require nonzero flow
ub = abs(ubulk(valc));
htl = htv;
htl(ub < 1e-8) = NaN; % drop vessels with no flow
[lo, ord_lo] = sort(htl, 'ascend');
valc_lo = valc(ord_lo);
valc_lo = valc_lo(~isnan(lo));

fid = fopen('headless_log.txt', 'a');
trials = {};
% Strategy A: constrict top K high-Ht vessels only
for K = [5 10 15 20 30 40]
    trials{end+1} = struct('name', sprintf('constrictHi_K%d',K), ...
                           'constriction', valc_hi(1:K), 'dilation', []);
end
% Strategy B: dilate top K low-Ht vessels only
for K = [5 10 15 20 30 40]
    trials{end+1} = struct('name', sprintf('dilateLo_K%d',K), ...
                           'constriction', [], 'dilation', valc_lo(1:K));
end
% Strategy C: combined
for K = [10 15 20 25 30 40]
    trials{end+1} = struct('name', sprintf('combo_K%d',K), ...
                           'constriction', valc_hi(1:K), 'dilation', valc_lo(1:K));
end

best_pct = -inf;
best_name = '';
best_con = [];
best_dil = [];

for k = 1:length(trials)
    tr = trials{k};
    t0 = tic;
    pct = competition_headless(tr.constriction, tr.dilation);
    el = toc(t0);
    line = sprintf('%-22s nc=%d nd=%d  pct=%.4f  (%.1fs)', tr.name, ...
        numel(tr.constriction), numel(tr.dilation), pct, el);
    fprintf('%s\n', line);
    fprintf(fid, '%s\n', line);
    if pct > best_pct
        best_pct = pct;
        best_name = tr.name;
        best_con = tr.constriction;
        best_dil = tr.dilation;
    end
end

fprintf('BEST: %s pct=%.4f\n', best_name, best_pct);
fprintf(fid, 'BEST: %s pct=%.4f\n', best_name, best_pct);
fclose(fid);
save('best_v1.mat', 'best_pct', 'best_name', 'best_con', 'best_dil');
end
