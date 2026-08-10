clear


EOS = ["DD2_hot_equal_pinned" "MPA1"];

n_eos = numel(EOS);

figure(3); hold on; box on;
for i=1:n_eos
    dir = strcat("/Users/horay/ptmp/GRASS/eos/",EOS(i),".dat");
    A = importdata(dir," ",1);
    e = A.data(:,1);
    p = A.data(:,2);
    h = A.data(:,3);
    n = A.data(:,4);
    plot( e, p, '-o' )
end

set(gca,'xscale','log'); set(gca,'yscale','log');
legend('Interpreter', 'latex', 'FontSize', 20,'Location','northwest')

xlabel('$e$', 'Interpreter', 'latex');
ylabel('$p$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 25, ...
    'LabelFontSizeMultiplier', 1.1, 'TickDir', 'in');
