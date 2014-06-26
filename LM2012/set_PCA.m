% Set PCA fig properties
fontsize = 8;
fontname = 'Arial';

ylim([-pi pi]);
set(gca, 'XTick', [0 1], 'YTick', [-pi 0 pi], 'TickDir', 'out');
set(gca, 'YTickLabel', ['D';'V'; 'D']);
ylabel([]);
xlabel([]);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', fontsize);
set(gca, 'FontName', fontname);

caxis(lims);
colormap(bluewhitered(256));

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 xSize ySize]);
set(gcf, 'Color', 'w');
