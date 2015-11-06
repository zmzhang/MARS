function masstomsp(project_names,S,R,W)

[rows,cols]=size(S);
fid = fopen([project_names 'MS.msp'], 'wt'); % Open for writing
for i=1:rows
   fprintf(fid,'Name: scan %d \n',R(i));
   fprintf(fid,'Num peaks: %d \n',length(W));
for j=1:cols
   fprintf(fid,'%2.5f \t %d \n',W(j),S(i,j)); 
end
end
fclose(fid);
end
