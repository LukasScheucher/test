function matrix = verMatrix(matrices,n,m)

matrix = zeros(n,m);
ni = 1;
for i = 1:length(matrices)
    niend = ni+size(matrices{i},1) - 1;
    matrix(ni:niend,:) = matrices{i};
    ni = niend + 1;
end

end