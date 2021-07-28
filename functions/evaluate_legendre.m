function Plm = evaluate_legendre(Nmax, x, fnpl)
Plm = NaN(length(x),addmup(Nmax)); % ( Nlat, (Nmax+1)*(Nmax+2)/2 )
if Nmax>200
    h=waitbar(0,'Evaluating all Legendre polynomials');
end

in2=0;
for in=0:Nmax
    in1 = in2+1;
    in2 = in2+in+1;
    Plm(:,in1:in2)=(legendre(in,x(:)','sch')*sqrt(2*in+1))';
    
    if Nmax>200
        waitbar((in+1)/(Nmax+1),h)
    end
end
if Nmax>200
    delete(h)
end
if length(x) > 90
    save(fnpl,'Plm')
end
end