% Subfunction for single-trajectory of deterministic evolution ODE 
% Fast: Constant Hamiltonian
function unpsidot = dunpsifast(unpsi,hbar,H_eff)
% No jump
unpsidot = -1i*hbar*H_eff*unpsi;
