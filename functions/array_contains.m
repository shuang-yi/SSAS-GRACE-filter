function [a_share,a_only,b_only,b_share,ind_a,ind_b]=array_contains(ia,ib,varargin)
% check the share part between a and b
% 
% See also remove_by_seq
% {{array; contain; sequence}}

if nargin == 0
    help(mfilename);
    fun_ex1;
    return;
end

PAR0 = struct('misfit',0);
PAR = var_initial(PAR0, varargin);

% --------------------end of head-----------
ia = ia(:);
ib = ib(:);
dd = abs(ia-ib');
ind_b = min(dd,[],1) <= PAR.misfit;
b_share = ib(ind_b);
b_only = ib(~ind_b);

ind_a = min(dd,[],2) <= PAR.misfit;
a_share = ia(ind_a);
a_only = ia(~ind_a);

end

%% example

function fun_ex1()
ia = 1:6
ib = [3,5,6]+0.1

[b_share,b_only,a_only,a_share]=array_contains(ia,ib,'misfit',0)
end