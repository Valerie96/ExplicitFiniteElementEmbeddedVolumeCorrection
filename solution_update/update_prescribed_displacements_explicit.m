%--------------------------------------------------------------------------
%  Update coodinates of displacement (Dirichlet) boundary conditions.
%--------------------------------------------------------------------------
function [x, v]       = update_prescribed_displacements_explicit(dofprescribed,x0,x,...
                    v, presc_displacement,dt,amp)

%Apply the displacements as a linear ramp over the full simulation time
AppliedDisp = presc_displacement(dofprescribed);
xnew = x0(dofprescribed) + AppliedDisp*amp;


v(dofprescribed) = (xnew-x(dofprescribed))/(dt);
x(dofprescribed) = x0(dofprescribed) + AppliedDisp*amp;
