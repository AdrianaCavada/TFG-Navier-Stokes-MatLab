function [M] = avg(m)
%AVG Summary of this function goes here
%   Detailed explanation goes here
[size_i,size_j]=size(m);

for i=1:size_i
    for j=1:size_j-1
        M(i,j)=(m(i,j)+m(i,j+1))/2;
    end
end

end

