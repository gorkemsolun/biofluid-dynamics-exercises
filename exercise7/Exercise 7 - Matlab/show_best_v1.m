function show_best_v1()
B = load('best_v1.mat');
fprintf('Best: %s pct=%.4f\n', B.best_name, B.best_pct);
fprintf('Constriction (n=%d): %s\n', numel(B.best_con), mat2str(B.best_con(:)'));
fprintf('Dilation     (n=%d): %s\n', numel(B.best_dil), mat2str(B.best_dil(:)'));
end
