function [ICC_ck,assign] = matrixICC(Matrix1,Matrix2,reorder)

addpath(genpath('~/Dropbox/Matlab/ICC'))
if size(Matrix1,2)~=size(Matrix2,2)
    warning('Matrices are not the same size, some columns will be ignored\n')
end

ICC_ck = zeros(size(Matrix1,2));
for i = 1:size(Matrix1,2)
    for j = 1:size(Matrix1,2)
        ICC_ck(i,j) = ICC([Matrix1(:,i) Matrix2(:,j)],'A-1');
    end
end

if reorder==1
    assign = munkres(-1*ICC_ck);
    Matrix2 = Matrix2(:,assign);
    for i = 1:size(Matrix1,2)
        for j = 1:size(Matrix1,2)
            ICC_ck(i,j) = ICC([Matrix1(:,i) Matrix2(:,j)],'A-1');
        end
    end
else
    assign = [];
end




