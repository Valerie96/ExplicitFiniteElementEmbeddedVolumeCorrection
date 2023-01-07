%--------------------------------------------------------------------------
%  Update coodinates of displacement (Dirichlet) boundary conditions.
%--------------------------------------------------------------------------
function [x, v]       = update_prescribed_displacements_explicit(dofprescribed,x0,x,...
                    v, presc_displacement,time_n, total_time)

%Apply the displacements as a linear ramp over the full simulation time
AppliedDisp = presc_displacement(dofprescribed);
ramp = time_n * (AppliedDisp / total_time);

x(dofprescribed) = x0(dofprescribed) + ramp;
v(dofprescribed) = (AppliedDisp / total_time);