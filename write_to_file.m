function write_to_file(matrix, L_1, L_2, L, H, step)
fid=fopen(sprintf('config_%d.txt',step),'w');
fprintf(fid,'ITEM: TIMESTEP\n');
aa=size(matrix);aa=aa(1);
fprintf(fid,'%.0f\n',step);
fprintf(fid,'ITEM: NUMBER OF ATOMS\n');
fprintf(fid,'%.0f\n',aa);
fprintf(fid,'ITEM: BOX BOUNDS pp pp pp\n');
% print box information to txt file
for x=1:1:3
    fprintf(fid,'%.0f %.0f\n',L,H);
end
fprintf(fid,'ITEM: ATOMS id type xu yu zu \n');
% print position information to txt file

type = 0;
for x=1:1:aa
    if (matrix(x,2) <= L_1)
        type = 1;
    elseif (matrix(x,2) < (L-L_2+1))
        type = 2;
    else
        type = 3;
    end
    fprintf(fid,'%.0f %.0f %.6f %.6f %.6f\n',aa,type,matrix(x,2),matrix(x,1),0);
end
fclose(fid);