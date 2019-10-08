%function drhodt = lindblad(~, rho, Hsd, natom, g, beta, wc,A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,v,e)

%function needed for ode intergration: calculate the lindblad terms
function drhodt = lindblad(~, rho, Hsd, natom, gsq2pi, beta, betainv, wc,A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,v,e)
neval = 2^natom;
nevaltruc = 18;
rhomcb = sparse(reshape(rho, [nevaltruc,nevaltruc]));

drhodt = -1i*(Hsd*rhomcb - rhomcb*Hsd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = bsxfun(@minus, e.', e);

%sorted jump operators
[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol
length(sortedOutput(sortedOutput>0));
count = 0;
%Lindblad for w>0
for w = sortedOutput(sortedOutput>0)
    gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
    if isnan(gamma) || isinf(gamma)
        gamma = gsq2pi*betainv;
    end
    [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
    count = count+length(b);

    %use sparse matrix
    Lpcomponents1 = sparse(nevaltruc,nevaltruc); 
    Lpcomponents2 = sparse(nevaltruc,nevaltruc);
    Lpcomponents3 = sparse(nevaltruc,nevaltruc);
    Lpcomponents4 = sparse(nevaltruc,nevaltruc);
    Lpcomponents5 = sparse(nevaltruc,nevaltruc);
    Lpcomponents6 = sparse(nevaltruc,nevaltruc);
    Lpcomponents7 = sparse(nevaltruc,nevaltruc);
    Lpcomponents8 = sparse(nevaltruc,nevaltruc);
    for s = 1:length(b)
      %matrixelement = v(:,i(s))'*A*v(:,j(s));
      matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));     %j<->a, i<->b in paper
      matrixelement2 = v(:,a(s))'*A_2*v(:,b(s));
      matrixelement3 = v(:,a(s))'*A_3*v(:,b(s));
      matrixelement4 = v(:,a(s))'*A_4*v(:,b(s));
      matrixelement5 = v(:,a(s))'*A_5*v(:,b(s));
      matrixelement6 = v(:,a(s))'*A_6*v(:,b(s));
      matrixelement7 = v(:,a(s))'*A_7*v(:,b(s));
      matrixelement8 = v(:,a(s))'*A_8*v(:,b(s));

      Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent2 = matrixelement2*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent3 = matrixelement3*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent4 = matrixelement4*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent5 = matrixelement5*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent6 = matrixelement6*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent7 = matrixelement7*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
      Lpcomponent8 = matrixelement8*sparse(a(s),b(s),1,nevaltruc,nevaltruc);

      Lpcomponents1 = Lpcomponents1 + Lpcomponent1;
      Lpcomponents2 = Lpcomponents2 + Lpcomponent2;
      Lpcomponents3 = Lpcomponents3 + Lpcomponent3;
      Lpcomponents4 = Lpcomponents4 + Lpcomponent4;
      Lpcomponents5 = Lpcomponents5 + Lpcomponent5;
      Lpcomponents6 = Lpcomponents6 + Lpcomponent6;
      Lpcomponents7 = Lpcomponents7 + Lpcomponent7;
      Lpcomponents8 = Lpcomponents8 + Lpcomponent8;
    end

    Lncomponents1 = Lpcomponents1';
    Lncomponents2 = Lpcomponents2';
    Lncomponents3 = Lpcomponents3';
    Lncomponents4 = Lpcomponents4';
    Lncomponents5 = Lpcomponents5';
    Lncomponents6 = Lpcomponents6';
    Lncomponents7 = Lpcomponents7';
    Lncomponents8 = Lpcomponents8';
    
    Lp1 = sqrt(gamma)*Lpcomponents1;
    Lp2 = sqrt(gamma)*Lpcomponents2;
    Lp3 = sqrt(gamma)*Lpcomponents3;
    Lp4 = sqrt(gamma)*Lpcomponents4;
    Lp5 = sqrt(gamma)*Lpcomponents5;
    Lp6 = sqrt(gamma)*Lpcomponents6;
    Lp7 = sqrt(gamma)*Lpcomponents7;
    Lp8 = sqrt(gamma)*Lpcomponents8;

    Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;
    Ln2 = sqrt(gamma*exp(-beta*w))*Lncomponents2;
    Ln3 = sqrt(gamma*exp(-beta*w))*Lncomponents3;
    Ln4 = sqrt(gamma*exp(-beta*w))*Lncomponents4;
    Ln5 = sqrt(gamma*exp(-beta*w))*Lncomponents5;
    Ln6 = sqrt(gamma*exp(-beta*w))*Lncomponents6;
    Ln7 = sqrt(gamma*exp(-beta*w))*Lncomponents7;
    Ln8 = sqrt(gamma*exp(-beta*w))*Lncomponents8;   
    
    


%%%%%%%%%%Superoperator2%%%%%%%%%  
drhodt = drhodt + (Lp1*rhomcb*Lp1'-0.5*(Lp1)'*Lp1*rhomcb - 0.5*rhomcb*(Lp1)'*Lp1) + (Ln1*rhomcb*Ln1'-0.5*(Ln1)'*Ln1*rhomcb - 0.5*rhomcb*(Ln1)'*Ln1)...
    + (Lp2*rhomcb*Lp2'-0.5*(Lp2)'*Lp2*rhomcb - 0.5*rhomcb*(Lp2)'*Lp2) + (Ln2*rhomcb*Ln2'-0.5*(Ln2)'*Ln2*rhomcb - 0.5*rhomcb*(Ln2)'*Ln2) ...
    + (Lp3*rhomcb*Lp3'-0.5*(Lp3)'*Lp3*rhomcb - 0.5*rhomcb*(Lp3)'*Lp3) + (Ln3*rhomcb*Ln3'-0.5*(Ln3)'*Ln3*rhomcb - 0.5*rhomcb*(Ln3)'*Ln3)...
    + (Lp4*rhomcb*Lp4'-0.5*(Lp4)'*Lp4*rhomcb - 0.5*rhomcb*(Lp4)'*Lp4) + (Ln4*rhomcb*Ln4'-0.5*(Ln4)'*Ln4*rhomcb - 0.5*rhomcb*(Ln4)'*Ln4) ...
    + (Lp5*rhomcb*Lp5'-0.5*(Lp5)'*Lp5*rhomcb - 0.5*rhomcb*(Lp5)'*Lp5) + (Ln5*rhomcb*Ln5'-0.5*(Ln5)'*Ln5*rhomcb - 0.5*rhomcb*(Ln5)'*Ln5) ...
    + (Lp6*rhomcb*Lp6'-0.5*(Lp6)'*Lp6*rhomcb - 0.5*rhomcb*(Lp6)'*Lp6) + (Ln6*rhomcb*Ln6'-0.5*(Ln6)'*Ln6*rhomcb - 0.5*rhomcb*(Ln6)'*Ln6) ...
    + (Lp7*rhomcb*Lp7'-0.5*(Lp7)'*Lp7*rhomcb - 0.5*rhomcb*(Lp7)'*Lp7) + (Ln7*rhomcb*Ln7'-0.5*(Ln7)'*Ln7*rhomcb - 0.5*rhomcb*(Ln7)'*Ln7) ...
    + (Lp8*rhomcb*Lp8'-0.5*(Lp8)'*Lp8*rhomcb - 0.5*rhomcb*(Lp8)'*Lp8) + (Ln8*rhomcb*Ln8'-0.5*(Ln8)'*Ln8*rhomcb - 0.5*rhomcb*(Ln8)'*Ln8);
end

gamma0 = gsq2pi*betainv;
%[b0, a0] = ind2sub(size(X), find(X == 0));
[b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol

L01 = sparse(nevaltruc,nevaltruc);
L02 = sparse(nevaltruc,nevaltruc);
L03 = sparse(nevaltruc,nevaltruc);
L04 = sparse(nevaltruc,nevaltruc);
L05 = sparse(nevaltruc,nevaltruc);
L06 = sparse(nevaltruc,nevaltruc);
L07 = sparse(nevaltruc,nevaltruc);
L08 = sparse(nevaltruc,nevaltruc);
for s = 1:length(b0)
    matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));
    matrixelement2 = v(:,a0(s))'*A_2*v(:,b0(s));
    matrixelement3 = v(:,a0(s))'*A_3*v(:,b0(s));
    matrixelement4 = v(:,a0(s))'*A_4*v(:,b0(s));
    matrixelement5 = v(:,a0(s))'*A_5*v(:,b0(s));
    matrixelement6 = v(:,a0(s))'*A_6*v(:,b0(s));
    matrixelement7 = v(:,a0(s))'*A_7*v(:,b0(s));
    matrixelement8 = v(:,a0(s))'*A_8*v(:,b0(s));
    L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component2 = matrixelement2*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component3 = matrixelement3*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component4 = matrixelement4*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component5 = matrixelement5*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component6 = matrixelement6*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component7 = matrixelement7*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L0component8 = matrixelement8*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
    L01 = L01 + L0component1;
    L02 = L02 + L0component2;
    L03 = L03 + L0component3;
    L04 = L04 + L0component4;
    L05 = L05 + L0component5;
    L06 = L06 + L0component6;
    L07 = L07 + L0component7;
    L08 = L08 + L0component8;
end

L01 = sqrt(gamma0)*L01;
L02 = sqrt(gamma0)*L02;
L03 = sqrt(gamma0)*L03;
L04 = sqrt(gamma0)*L04;
L05 = sqrt(gamma0)*L05;
L06 = sqrt(gamma0)*L06;
L07 = sqrt(gamma0)*L07;
L08 = sqrt(gamma0)*L08;


%%%%%%%%%Superoperator2%%%%%%%%%    
drhodt = drhodt + (L01*rhomcb*(L01)'-0.5*(L01)'*L01*rhomcb - 0.5*rhomcb*(L01)'*L01)...
    + (L02*rhomcb*(L02)'-0.5*(L02)'*L02*rhomcb - 0.5*rhomcb*(L02)'*L02)...
    + (L03*rhomcb*(L03)'-0.5*(L03)'*L03*rhomcb - 0.5*rhomcb*(L03)'*L03)...
    + (L04*rhomcb*(L04)'-0.5*(L04)'*L04*rhomcb - 0.5*rhomcb*(L04)'*L04)...
    + (L05*rhomcb*(L05)'-0.5*(L05)'*L05*rhomcb - 0.5*rhomcb*(L05)'*L05)...
    + (L06*rhomcb*(L06)'-0.5*(L06)'*L06*rhomcb - 0.5*rhomcb*(L06)'*L06)...
    + (L07*rhomcb*(L07)'-0.5*(L07)'*L07*rhomcb - 0.5*rhomcb*(L07)'*L07)...
    + (L08*rhomcb*(L08)'-0.5*(L08)'*L08*rhomcb - 0.5*rhomcb*(L08)'*L08);


%%%%%%%%%Superoperator2%%%%%%%%%
% drhodt = L*rhovcb;
%%%%%%%%%Superoperator2%%%%%%%%%
%whos
drhodt = drhodt(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
