function run_diag()
% Baseline diagnostics: per-vessel HtAvg, flow, position
t0 = tic;
[pct, valc, c, p, HtAvg] = competition_headless([], []);
fprintf('DIAG pct=%.4f elapsed=%.1fs\n', pct, toc(t0));

nc = size(c,1);
% Midpoint coordinates of each vessel
xmid = 0.5*(p(c(:,1),1) + p(c(:,2),1));
ymid = 0.5*(p(c(:,1),2) + p(c(:,2),2));
ubulk = c(:,5);
a = c(:,8);

T = table((1:nc)', xmid, ymid, a, ubulk, HtAvg, 'VariableNames', ...
    {'id','xmid','ymid','radius','ubulk','HtAvg'});
save('diag.mat', 'c', 'p', 'valc', 'HtAvg', 'xmid', 'ymid', 'ubulk');
writetable(T, 'diag.csv');
fprintf('Saved diag.mat and diag.csv\n');
end
