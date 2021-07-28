function sc2=unit_convert(sc, from, to)
lmax = size(sc,1)-1;
llist = (0:lmax)';
tf1 = eigengrav(llist, from);
tf2 = eigengrav(llist, to);

rr = tf2./tf1;

rr(~isfinite(rr)) = 0;

sc2  = sc .* repmat(rr, 1, 2 * lmax + 1);
end