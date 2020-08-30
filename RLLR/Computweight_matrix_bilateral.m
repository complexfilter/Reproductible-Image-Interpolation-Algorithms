function weight=Computweight_matrix_bilateral(CB,X,location_diff)
d = (CB-X).*(CB-X).*0.05;
dd = sqrt(sum(d));
diff_location = sqrt(sum(location_diff .* location_diff.*0.004));
weight = exp(-dd).*exp(-diff_location);
weight = weight./sum(weight);
weight = weight';
end