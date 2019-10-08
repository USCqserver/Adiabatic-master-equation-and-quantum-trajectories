%kawayip@usc.edu
%8-qubit_chain
%function ame use the function lindblad on line 189 to run. 
function ame
%warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym')
% % Define Pauli matrices
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
unit = speye(2);
natom = 8;

%nevel means number of eigenvalues
neval = 2^natom;
%nevaltruc means number of truncations
nevaltruc = 18;

% Define tensor product of Pauli matrices
sX_1 = kron(sX, speye(2^(natom-1)));
sX_2 = kron(kron(speye(2),sX),speye(2^(natom-2)));
sX_3 = kron(kron(speye(2^2),sX),speye(2^(natom-3)));
sX_4 = kron(kron(speye(2^3),sX),speye(2^(natom-4)));
sX_5 = kron(kron(speye(2^4),sX),speye(2^(natom-5)));
sX_6 = kron(kron(speye(2^5),sX),speye(2^(natom-6)));
sX_7 = kron(kron(speye(2^6),sX),speye(2^(natom-7)));
sX_8 = kron(speye(2^(natom-1)),sX);

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

sZsZIIIIII = kron(kron(kron(kron(kron(kron(kron(sZ,sZ),unit),unit),unit),unit),unit),unit);
IsZsZIIIII = kron(kron(kron(kron(kron(kron(kron(unit,sZ),sZ),unit),unit),unit),unit),unit);
IIsZsZIIII = kron(kron(kron(kron(kron(kron(kron(unit,unit),sZ),sZ),unit),unit),unit),unit);
IIIsZsZIII = kron(kron(kron(kron(kron(kron(kron(unit,unit),unit),sZ),sZ),unit),unit),unit);
IIIIsZsZII = kron(kron(kron(kron(kron(kron(kron(unit,unit),unit),unit),sZ),sZ),unit),unit);
IIIIIsZsZI = kron(kron(kron(kron(kron(kron(kron(unit,unit),unit),unit),unit),sZ),sZ),unit);
IIIIIIsZsZ = kron(kron(kron(kron(kron(kron(kron(unit,unit),unit),unit),unit),unit),sZ),sZ);


plus = [1/sqrt(2); 1/sqrt(2)];

%these lines mean to read from file
dlm = dlmread('DW1_parameters.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).';
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);


%Temperature and coupling parameters
beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.27)*1e-4*2*pi; %g^2*pi
betainv = 2.6*1e9; %1/beta


%Gibbs state initialization
Hs = -1e9.*A_sp1(0).*(sX_1+sX_2+sX_3+sX_4+sX_5+sX_6+sX_7+sX_8) + 1e9.*B_sp1(0).*((-1).*((1/4).*sZ_1)...
 + ((-1).*sZsZIIIIII + (-1).*IsZsZIIIII + (-1).*IIsZsZIIII + (-1).*IIIsZsZIII + (-1).*IIIIsZsZII + (-1).*IIIIIsZsZI...
 + (-1).*IIIIIIsZsZ));
[V,D] = eig(full(Hs));
if ~issorted(diag(D))
    [V,D] = eig(full(Hs));
    [D,I] = sort(diag(D));
    D = diag(D);
    V = V(:, I);
    fprintf('sorted');
end
%e, v are eigenvalus and eigenvectors
e = zeros(1,nevaltruc);
v = sparse(V(:,1:nevaltruc));
for i = 1:nevaltruc
    e(i) = sparse(D(i,i));
end
Z = 0;
for i = 1:nevaltruc
    Z = Z + exp(-beta*e(i));
end

gibbs = 0;
for i = 1:nevaltruc
    gibbs = gibbs + (exp(-beta*e(i))/Z)*v(:,i)*v(:,i)';
end

% %Pure state initialization
% psi0 = sparse(kron(kron(kron(kron(kron(kron(kron(plus,plus),plus),plus),plus),plus),plus),plus));
% rho0 = psi0*psi0';
% %rho = rho0(:);
rho = gibbs(:);

tf = 10e-6;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timestep between each ode executions
dt_me = tf/1000;
tstep_me = 0:dt_me:tf;
fidelitylist_me = zeros(1, numel(tstep_me));

%a for loop inside each step use ode solver
for index = 1:numel(tstep_me)
    Hs = -1e9.*A_sp1(tstep_me(index)./tf).*(sX_1+sX_2+sX_3+sX_4+sX_5+sX_6+sX_7+sX_8) + 1e9.*B_sp1(tstep_me(index)./tf).*((-1).*((1/4).*sZ_1)...
     + ((-1).*sZsZIIIIII + (-1).*IsZsZIIIII + (-1).*IIsZsZIIII + (-1).*IIIsZsZIII + (-1).*IIIIsZsZII + (-1).*IIIIIsZsZI...
     + (-1).*IIIIIIsZsZ));

    [V,D] = eig(full(Hs));
    %make sure eigenmatrix is sorted
    if ~issorted(diag(D))
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        fprintf('sorted');
    end
    %e, v eigenvalues and eigenvectors.
    v = sparse(neval,nevaltruc);
    e = sparse(1,nevaltruc);
    
    %make them sparse
    for ii = 1:nevaltruc
        v(:,ii) = sparse(V(:,ii));
        e(ii) = sparse(D(ii,ii));
    end
    
    %diagonalized H and rho
    Hsd = v'*Hs*v;
    v0 = v(:,1);
    rhom = reshape(rho,[neval,neval]);
    rhomcb = sparse(v'*rhom*v);
    
    %ground state population
    fidelity = rhomcb(1,1);
    fidelitylist_me(1, index) = fidelity;
    
    %ode solver
    alindblad  = @(t, rho)lindblad(t, rho, Hsd, natom, gsq2pi, beta, betainv, wc, A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,v,e);
    t    = [tstep_me(index), tstep_me(index) + dt_me];
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    %[~, RHO] = ode45(alindblad, t, rhomcb(:));
    [~, RHO] = ode45(alindblad, t, rhomcb(:),options);
    % Final value of rho is initial value for next step:
    rho = RHO(end, :); 
    rho = v*reshape(rho,[nevaltruc,nevaltruc])*v';   %in computational basis
    rho = rho(:);
end

%time taken to run
eptime = toc

%plot fidelity
figure(1)
plot(tstep_me, fidelitylist_me,'-b','LineWidth',2);
xlabel('$t$','Interpreter','latex')
xlim([0 tf])
ylabel('$fidelity$','Interpreter','latex')
title(['tf: ' num2str(tf)])

figure(2)
plot(tstep_me./tf, fidelitylist_me,'-b','LineWidth',2); 
xlabel('$s$','Interpreter','latex')
ylabel('$fidelity$','Interpreter','latex')
title(['tf: ' num2str(tf)])

%store fidelity
txt1 = sprintf('ame%d.txt',tf);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d\n',[tstep_me;fidelitylist_me]);
fclose(fid1);


