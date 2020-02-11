function export_fig_pdf(file, BLACK)
% export_fig_pdf( file, BLACK )
% taken from the matlab website!
% BLACK: if true, then turn off 'inverthardcopy' which by default renders 
%    black backgrounds as white. 
% sgm  2016
if ~exist('BLACK','var'), BLACK = false; end
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
if BLACK
  set(h,'inverthardcopy','off');
end
print(gcf,'-painters','-dpdf',file) ; 
% help page also suggests '-r0' flag