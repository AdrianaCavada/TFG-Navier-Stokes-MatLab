function [] = print_matrix(U,V,p,errorU, errorV,n_step,time)
%PRINT_MATRIX Summary of this function goes here
%   Detailed explanation goes here

[rowsU colsU]=size(U);
% x=repmat('%10.6f\n',1,cols-1);
fileID= fopen('U.txt','wt');
for i=1:rowsU
    for j=1:colsU
        fprintf(fileID,'%f\t', U(i,j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

[rowsV colsV]=size(V);
% x=repmat('%10.6f\n',1,cols-1);
fileID= fopen('V.txt','wt');
for i=1:rowsV
    for j=1:colsV
        fprintf(fileID,'%f\t', V(i,j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

[rowsp colsp]=size(p);
% x=repmat('%10.6f\n',1,cols-1);
fileID= fopen('p.txt','wt');
for i=1:rowsp
    for j=1:colsp
        fprintf(fileID,'%f\t', p(i,j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);
fileID= fopen('errors.txt','wt');
fprintf(fileID,'%f\t', errorU);
fprintf(fileID, '\n');
fprintf(fileID,'%f\t', errorV);
fprintf(fileID, '\n');
fprintf(fileID,'%f\t', n_step);
fprintf(fileID, '\n');
fprintf(fileID,'%f\t', time);
fprintf(fileID, '\n');
fclose(fileID);



end

