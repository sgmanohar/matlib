function gifti2obj( giftiStructure, outfile)
% write object text file (lightwave .obj format)
% sgm
if exist(outfile,'file')
  reply=input('Overwrite file Y/N ?','s')
  if ~(reply=='Y' || reply=='y') return; end
end

f=fopen(outfile, 'w');
fprintf(f,'# gifti from matlab\n');
for i=1:size(giftiStructure.vertices,1)
  x=giftiStructure.vertices(i,:);
  fprintf(f, 'v %g %g %g\n', x(1), x(2), x(3) );
end
for i=1:size(giftiStructure.faces,1)
  x=giftiStructure.faces(i,:);
  fprintf(f, 'f %g %g %g\n', x(1), x(2), x(3) );
end
fclose(f);
