%subfunction to search for lindblad
function L = lindbladsearch(k,natom,v,e,nevaltruc)
sZ = [1 0; 0 -1];
z = ceil(k/8);

beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi;
betainv = 2.6*1e9;

sZ_1 = kron(sZ, speye(2^(natom-1)));
sZ_2 = kron(kron(speye(2),sZ),speye(2^(natom-2)));
sZ_3 = kron(kron(speye(2^2),sZ),speye(2^(natom-3)));
sZ_4 = kron(kron(speye(2^3),sZ),speye(2^(natom-4)));
sZ_5 = kron(kron(speye(2^4),sZ),speye(2^(natom-5)));
sZ_6 = kron(kron(speye(2^5),sZ),speye(2^(natom-6)));
sZ_7 = kron(kron(speye(2^6),sZ),speye(2^(natom-7)));
sZ_8 = kron(speye(2^(natom-1)),sZ);

A_1 = sZ_1;
A_2 = sZ_2;
A_3 = sZ_3;
A_4 = sZ_4;
A_5 = sZ_5;
A_6 = sZ_6;
A_7 = sZ_7;
A_8 = sZ_8;


X = bsxfun(@minus, e.', e);
%[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol

if z == 1
   column = ceil(length(sortedOutput)/2);   
elseif logical(mod(z,2))
   column = ceil(length(sortedOutput)/2) - (z-1)/2;
else
   column = ceil(length(sortedOutput)/2) + z/2;
end



w = sortedOutput(column);
gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
if isnan(gamma) || isinf(gamma)
    %gamma = 2.*g.^2.*pi.*beta.^(-1);
    gamma = gsq2pi*betainv;
end
[b, a] = ind2sub(size(X), find(abs(X - w)< 0.1));% AbsTol
L = zeros(nevaltruc,nevaltruc); 
if mod(k,8) == 1
    A = A_1;
elseif mod(k,8) == 2
    A = A_2;
elseif mod(k,8) == 3
    A = A_3;
elseif mod(k,8) == 4
    A = A_4;
elseif mod(k,8) == 5
    A = A_5;
elseif mod(k,8) == 6
    A = A_6;
elseif mod(k,8) == 7
    A = A_7;
else
    A = A_8;
end

for s = 1:length(b)
  matrixelement = v(:,a(s))'*A*v(:,b(s));
  Lcomponent = matrixelement*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
  L = L + Lcomponent;
end
L = sqrt(gamma)*L;
