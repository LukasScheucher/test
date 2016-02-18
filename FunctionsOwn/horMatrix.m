function matrix = horMatrix(matrices,n,m)

matrix = zeros(n,m);
mi = 1;
for i = 1:length(matrices)
    miend = mi+size(matrices{i},2) - 1;
    matrix(:,mi:miend) = matrices{i};
    mi = miend + 1;
end

end