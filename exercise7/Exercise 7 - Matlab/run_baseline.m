function run_baseline()
% Baseline: no modifications
t0 = tic;
[pct, valc, c] = competition_headless([], []);
fprintf('BASELINE pct=%.4f  elapsed=%.1fs  nc=%d  |valc|=%d\n', pct, toc(t0), size(c,1), length(valc));
fid = fopen('headless_log.txt', 'a');
fprintf(fid, 'BASELINE pct=%.4f  elapsed=%.1fs  nc=%d  valc=%d\n', pct, toc(t0), size(c,1), length(valc));
fclose(fid);
save('valc.mat', 'valc');
end
