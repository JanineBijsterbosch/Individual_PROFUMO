function map = dscalar2double(dscalar_struct,remove_nans)

N = fieldnames(dscalar_struct);
ismap = strncmp('x',N,1);
map = zeros(size(dscalar_struct.pos,1),sum(ismap)); i = 1;
for n = 1:length(N)
    if ismap(n) == 1
        map(:,i) = dscalar_struct.(N{n});
        i = i+1;
    end
end
if remove_nans==1
    map(isnan(map(:,1)),:) = [];
end
