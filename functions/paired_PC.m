function [ipaired,inot]=paired_PC(PC,ithreshold,max_shift)
% {{check; pair; correlation}}
K = size(PC,2);
max_corr = zeros(K);
min_corr = zeros(K);
for ii = 1:K
    i1 = min([K,ii+2]);
    for ij = ii+1:i1
        [~,~,~,max_corr(ii,ij),min_corr(ii,ij)] = shift_corr(PC(:,ii),PC(:,ij),max_shift);
    end   
end

[a,b] = find(max_corr>=ithreshold & min_corr<=-ithreshold);
ipaired = unique([a;b]);

[~,inot] = array_contains(1:K,ipaired);

end