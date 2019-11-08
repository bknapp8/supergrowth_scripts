function [S, V, L, W] = CENTERLINE_geom_props(x,y,stp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the centerline produced by CENTERLINE_pombe() and
% calculates surface area (S), volume (V), length (L), and median width (W)
% of a cell, described by the 

[geom, iner, cpmo] = polygeom(x,y);

center = [geom(2) geom(3)];

[M, W1, W2] = CENTERLINE_pombe(x,y,center, stp);


S = 0;
V = 0;
L = 0;
W = 0;

for k=1:length(M)-1
    m1 = M(k,:);
    m2 = M(k+1,:);
    w11 = W1(k,:);
    w12 = W1(k+1,:);
    w21 = W2(k,:);
    w22 = W2(k+1,:);

    midy1 = (norm(w11-m1) + norm(w12-m2))/2;
    midy2 = (norm(w21-m1) + norm(w22-m2))/2;

    dx = norm(m2-m1);

    % Make surface area
    S = S + (midy1+midy2)*pi*dx;

    % Make volume
    V = V + (midy1^2 + midy2^2)*pi/2*dx;

    % Make length
    L = L + dx;



end


% Make width
dw2 = (W2-W1).^2;
W = median(sqrt(sum(dw2')));
