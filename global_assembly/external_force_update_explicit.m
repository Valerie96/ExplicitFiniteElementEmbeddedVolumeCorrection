%--------------------------------------------------------------------------
% Update nodal forces (excluding pressure) and gravity.
%--------------------------------------------------------------------------
function [Residual,external_force] = external_force_update_explicit(total_load,...
                                     Residual,external_force,dlamb)
% Residual = Residual - dlamb*total_load;
external_force = dlamb*total_load;


