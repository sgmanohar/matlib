function export_fig_pdf(file)
% taken from the matlab website!
% sgm  2016
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-painters','-dpdf',file) ; 
% help page also suggests '-r0' flag