%--------------------------------------------------------------------------
% Calculate energy values and check for unrealisic values.
% Use trapazoidal integration to calulate internal, external, and kinetic 
% energy.
%--------------------------------------------------------------------------
function [Wint_n,Wext_n,WKE,Wvdamp_n,energy_value, max_energy] = check_energy_explicit(FEM,BC,GLOBAL,...
       disp_n, disp_prev,fi_n,fi_prev,fe_n,fe_prev,fvd_n,fvd_prev,Wint_n,Wext_n,Wvdamp_n,time)

%internal work
sum = (disp_n(BC.hostdof) - disp_prev(BC.hostdof))' * (fi_n(BC.hostdof) + fi_prev(BC.hostdof));
Wint_n = Wint_n + 0.5 * sum;

%external work
sum = (disp_n(BC.hostdof) - disp_prev(BC.hostdof))' * (fe_n(BC.hostdof) + fe_prev(BC.hostdof));
Wext_n = Wext_n + 0.5 * sum;

% kinetic energy
WKE = 0;
for i=1:FEM(1).mesh.n_dofs
    WKE = WKE + GLOBAL.M_KE(i,i) *( GLOBAL.velocities(i) * GLOBAL.velocities(i));
end
WKE = WKE * 0.5;

%energy from viscous damping
sum = (disp_n(BC.hostdof) - disp_prev(BC.hostdof))' * (fvd_n(BC.hostdof) + fvd_prev(BC.hostdof));
Wvdamp_n = Wvdamp_n + 0.5 * sum;

%calculate the total energy
energy_value = abs(WKE + Wint_n - Wext_n);
numbers = [WKE, Wint_n, Wext_n];
max_energy = max(numbers);
energy_tolerance = 0.01;
if(energy_value > (energy_tolerance * max_energy))
    disp('Energy Violation!')
end



