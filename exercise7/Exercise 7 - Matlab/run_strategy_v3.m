function run_strategy_v3()
% Refine the best from v2: run it, then add/swap candidates based on residual Ht
B = load('best_v2.mat');
con0 = B.best_con(:)';
dil0 = B.best_dil(:)';
fprintf('Starting from %s pct=%.4f (con=%d dil=%d)\n', B.best_name, B.best_pct, ...
    numel(con0), numel(dil0));

% Run once to get residual HtAvg
[base_pct, valc, c, p, HtAvg] = competition_headless(con0, dil0);
fprintf('Reconfirmed: pct=%.4f\n', base_pct);
fid = fopen('headless_log.txt', 'a');
fprintf(fid, '--- v3 ---\n');
fprintf(fid, 'Reconfirm c11_d14 pct=%.4f\n', base_pct);

ubulk = c(:,5);

% Identify modifiable vessels not yet touched, ranked by residual Ht
mods = unique([con0 dil0]);
free = setdiff(valc, mods);
ht_free = HtAvg(free);
ub_free = abs(ubulk(free));

% Sort residual high Ht (candidates for new constriction)
[~, oh] = sort(ht_free, 'descend');
cand_constrict = free(oh(1:min(15, numel(oh))));
% Sort residual low Ht with nonzero flow (candidates for new dilation)
ht_l = ht_free;
ht_l(ub_free < 1e-8) = NaN;
[lv, ol] = sort(ht_l, 'ascend');
mask = ~isnan(lv);
ol = ol(mask);
cand_dilate = free(ol(1:min(15, numel(ol))));

fprintf('Top residual high-Ht candidates to add (constrict): %s\n', mat2str(cand_constrict));
fprintf('Top residual low-Ht candidates to add (dilate):     %s\n', mat2str(cand_dilate));
fprintf(fid, 'cand_constrict=%s\n', mat2str(cand_constrict));
fprintf(fid, 'cand_dilate=%s\n', mat2str(cand_dilate));

trials = {};
% Add one new constriction at a time
for k = 1:numel(cand_constrict)
    new = [con0 cand_constrict(k)];
    trials{end+1} = struct('name', sprintf('addC_%d', cand_constrict(k)), ...
        'constriction', new, 'dilation', dil0);
end
% Add one new dilation at a time
for k = 1:numel(cand_dilate)
    new = [dil0 cand_dilate(k)];
    trials{end+1} = struct('name', sprintf('addD_%d', cand_dilate(k)), ...
        'constriction', con0, 'dilation', new);
end
% Remove one constriction at a time (drop test)
for k = 1:numel(con0)
    new = con0; new(k) = [];
    trials{end+1} = struct('name', sprintf('rmC_%d', con0(k)), ...
        'constriction', new, 'dilation', dil0);
end
% Remove one dilation at a time (drop test)
for k = 1:numel(dil0)
    new = dil0; new(k) = [];
    trials{end+1} = struct('name', sprintf('rmD_%d', dil0(k)), ...
        'constriction', con0, 'dilation', new);
end

best_pct = base_pct;
best_name = 'baseline_c11d14';
best_con = con0; best_dil = dil0;

for k = 1:length(trials)
    tr = trials{k};
    t0 = tic;
    pct = competition_headless(tr.constriction, tr.dilation);
    el = toc(t0);
    line = sprintf('%-20s nc=%d nd=%d  pct=%.4f  (%.1fs)', tr.name, ...
        numel(tr.constriction), numel(tr.dilation), pct, el);
    fprintf('%s\n', line);
    fprintf(fid, '%s\n', line);
    if pct > best_pct
        best_pct = pct; best_name = tr.name;
        best_con = tr.constriction; best_dil = tr.dilation;
    end
end

fprintf('BEST v3: %s pct=%.4f\n', best_name, best_pct);
fprintf(fid, 'BEST v3: %s pct=%.4f\n', best_name, best_pct);
fclose(fid);
save('best_v3.mat', 'best_pct', 'best_name', 'best_con', 'best_dil');
end
