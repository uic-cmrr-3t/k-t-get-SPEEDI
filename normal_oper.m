function UU = normal_oper(Uk_r, Vt_r, Rind, Cind, Mask)

N = size(Mask, 1);
M = size(Mask, 2);
r = size(Vt_r, 1);
G   = reshape(Uk_r, N, r)*Vt_r;
d = G(Mask);
clear G
UU = sparse(Rind, Cind, d, N, M)*Vt_r';
UU = reshape(UU, N*r, 1);
