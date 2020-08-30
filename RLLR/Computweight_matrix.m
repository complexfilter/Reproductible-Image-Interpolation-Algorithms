function weight=Computweight_matrix(CB,X)
d = (CB-X).*(CB-X).*0.05;
dd = sqrt(sum(d));
weight = exp(-dd);
weight = weight./sum(weight);
weight = weight';