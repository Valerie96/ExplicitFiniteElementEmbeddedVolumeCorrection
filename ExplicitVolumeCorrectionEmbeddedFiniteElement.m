%To run, either fill out the intro to this or run command
%ExplicitVolumeCorrectionEmbededFiniteElement(inputfile,outputfeq,vtuOn);

clear; clc; close all; 
basedir_fem='C:/Users/Valerie/Documents/GitHub/ExplicitFiniteElementEmbeddedVolumeCorrection/';
inputfile = 'Cube_8h_4t.dat';
inputfile = 'SmallTension_Speed.dat';
% inputfile = 'RussellTensile-Half_5000Fibers7_discritized.dat'
inputfile ='VolRedGuidelines_DyneemaCube_0fibers.dat';

outputfreq = 20;
vtuOn = 0;

ansmlv ='y'; 
prefactor = 0.7;

d1=digits(64);
CON.OUTPUT.incout = outputfreq;
CON.OUTPUT.nwant = 0;
CON.OUTPUT.iwant = 1;
% parpool('local');
tic
%% Input_data_and_initilaization.m
[PRO,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CONSTANT,GLOBAL,...
   PLAST,KINEMATICS,INITIAL_KINEMATICS,STRESS,simtime,Wint,Wext,Wvdamp,...
   EmbedElt,VolumeCorrect] = input_data_and_initialisation(basedir_fem,...
   ansmlv,inputfile,CON,vtuOn);

%% ExplicitDynamics_algorithm
%--------------------------------------------------------------------------
% Explicit central diff. algorithm 
% Based on Beylthereohetg Box xxx
%--------------------------------------------------------------------------

%Artifical bulk viscosity 
DAMPING.b1 = 0.04; %Linear bulk viscosity damping
DAMPING.b2 = 1.2; %Quadratic bulk viscosity damping

%Step 1 - Initialisation
%       - this is done in the intialisation.m file, line 68
velocities_half = zeros(FEM(1).mesh.n_dofs,1);
disp_n = zeros(FEM(1).mesh.n_dofs,1);
disp_prev = zeros(FEM(1).mesh.n_dofs,1);
CON.xlamb = 1;
CON.incrm = 0; 

output(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS,0,0);
CON.incrm = CON.incrm + 1; 

%Step 2 - Get Force
[GLOBAL,updated_PLAST,GEOM.Jn_1,GEOM.VolRate,f_damp] = getForce_explicit(CON.xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,INITIAL_KINEMATICS,BC,DAMPING,STRESS,EmbedElt,VolumeCorrect,1);      
     
%Step 3 - Compute accelerations.
GLOBAL.accelerations(BC.hostdof(:,1)) = inv(GLOBAL.M)*(GLOBAL.external_load(BC.hostdof(:,1)) - GLOBAL.T_int(BC.hostdof(:,1)));

%Step 4 - Time update/iterations
Time = 0; 
tMax = simtime; % in seconds
GLOBAL.tMax = tMax;
dt = prefactor * CalculateTimeStep(FEM,GEOM,MAT,DAMPING); % in seconds
time_step_counter = 0;

% Start explicit loop
while(Time<tMax)
    %Calculate time incremetation
    t_n       = Time;
    t_np1     = Time + dt;
    Time      = t_np1; % update the time by adding full time step
    dt_nphalf = dt; % equ 6.2.1
    t_nphalf  = 0.5 *(t_np1 + t_n); %equ 6.2.1
    
    %If on the final step, adjust dt to end the time at exactly the maximum
    %sepecified time
    if Time>=tMax
        fprintf('%d\n',Time);
        dt        = tMax - t_n;
        Time      = tMax;
        t_np1     = Time;
        dt_nphalf = dt;
        t_nphalf  = 0.5 *(t_np1 + t_n);
        fprintf('%d\n',Time);
    end        
    
% Step 5 - Update velocities
    velocities_half = GLOBAL.velocities + (t_nphalf - t_n) * GLOBAL.accelerations;
 
% Step 7 Update nodal displacments 
    % store old displacements for energy computation
    disp_prev = disp_n;
    % update nodal displacements 
    disp_n(BC.freedof) = disp_n(BC.freedof) + dt_nphalf *velocities_half(BC.freedof);
    
  %----------------------------------------------------------------
  % Update stored coodinates.
  displ = disp_n-disp_prev; 
  GEOM.x = update_geometry_explicit(GEOM.x,GEOM.x0,1,disp_n(BC.freedof),BC.freedof);
  %----------------------------------------------------------------

  % save external force, to be used in energy computation
  fe_prev = GLOBAL.external_load + GLOBAL.Reactions;
  
% Step 6 - enforce displacement BCs 
  %--------------------------------------------------------------------
  % Update nodal forces (excluding pressure) and gravity. 
  %--------------------------------------------------------------------
   CON.dlamb  = t_np1/tMax;
   [GLOBAL.Residual,GLOBAL.external_load] = external_force_update_explicit ...
       (GLOBAL.nominal_external_load,...
        GLOBAL.Residual,GLOBAL.external_load,CON.dlamb);

    GLOBAL.external_load_effective = GLOBAL.external_load + GLOBAL.Reactions;
    
  %--------------------------------------------------------------------
  % Update nodal forces and stiffness matrix due to external pressure 
  % boundary face (line) contributions. 
  %--------------------------------------------------------------------      
  if LOAD.n_pressure_loads      
     GLOBAL = pressure_load_and_stiffness_assembly(GEOM,MAT,FEM,...
              GLOBAL,LOAD,QUADRATURE.boundary,CON.dlamb);    
  end
  %--------------------------------------------------------------------
  % Update applied displacements (incrementation based on a smooth ramp 
  % over the total sim time 
  %--------------------------------------------------------------------
  if  BC.n_prescribed_displacements > 0
      [GEOM.x ,velocities_half]  = update_prescribed_displacements_explicit(BC.dofprescribed,...
               GEOM.x0,GEOM.x,velocities_half,BC.presc_displacement,t_np1,tMax); 
       disp_n(BC.fixdof) = GEOM.x(BC.fixdof) - GEOM.x0(BC.fixdof);  
  end
  %----------------------------------------------------------------
  % Update nodal forces due to pressure. 
  %----------------------------------------------------------------
  if LOAD.n_pressure_loads
      GLOBAL = pressure_load_and_stiffness_assembly(GEOM,MAT,FEM,...
               GLOBAL,LOAD,QUADRATURE.boundary,CON.xlamb);
  end
  
  %--------------------------------------------------------------------       
  % Update coodinates of embedded nodes (if there are any)  
  %--------------------------------------------------------------------
      if EmbedElt == 1
          [GEOM.x,velocities_half, GLOBAL.accelerations ] = update_embedded_displacements_explicit(BC.tiedof, BC.tienodes,...
                FEM,GEOM, velocities_half, GLOBAL.accelerations); 
          disp_n(BC.tiedof) = GEOM.x(BC.tiedof) - GEOM.x0(BC.tiedof);
      end
  %----------------------------------------------------------------   

  %Save internal force, to be used in energy computation
  fi_prev = GLOBAL.T_int; 
  f_damp_prev = f_damp;
  
%Step 8 - Get Force
  [GLOBAL,updated_PLAST,GEOM.Jn_1,GEOM.VolRate,f_damp] = getForce_explicit(CON.xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,INITIAL_KINEMATICS,BC,DAMPING,STRESS,EmbedElt,VolumeCorrect,dt);

  GLOBAL.external_load_effective = GLOBAL.external_load + GLOBAL.Reactions;
  
  % updated stable time increment based on current deformation     
  dt_old = dt;
  dt = prefactor * CalculateTimeStep(FEM,GEOM,MAT,DAMPING);

%Step 9 - Compute accelerations.       
  AccOld = GLOBAL.accelerations;
  GLOBAL.accelerations(BC.hostdof(:,1)) = inv(GLOBAL.M)*(GLOBAL.external_load_effective(BC.hostdof(:,1)) - GLOBAL.T_int(BC.hostdof(:,1)));
  GLOBAL.NetForce = (GLOBAL.external_load_effective - GLOBAL.T_int);
  
% Step 10 Second partial update of nodal velocities
  VelOld = GLOBAL.velocities;
  GLOBAL.velocities = velocities_half + (t_np1 - t_nphalf) * GLOBAL.accelerations;  
  %Update velocity and acceleration of embedded nodes (if there are any)  
      if EmbedElt == 1
          [GEOM.x,GLOBAL.velocities, GLOBAL.accelerations ] = update_embedded_displacements_explicit(BC.tiedof, BC.tienodes,...
                FEM,GEOM, GLOBAL.velocities, GLOBAL.accelerations); 
           disp_n(BC.tiedof) = GEOM.x(BC.tiedof) - GEOM.x0(BC.tiedof);
      end
  
 %--------------------------------------------------------------
% Step 11 Check energy
  [Wint,Wext,WKE,Wvdamp,energy_value, max_energy] = check_energy_explicit(FEM,BC,GLOBAL,...
      disp_n, disp_prev,GLOBAL.T_int,fi_prev,...
      GLOBAL.external_load + GLOBAL.Reactions,fe_prev,f_damp,f_damp_prev,Wint,Wext,Wvdamp,Time);
  %--------------------------------------------------------------
  %Plot every # steps
  if( mod(time_step_counter,outputfreq) == 0 || Time == tMax)
      output(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS,Time,dt);
      write_energy_output(PRO,CON,Wint,Wext,WKE,Wvdamp,Time);
      if vtuOn
        output_vtu(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS);
      end
%       PLAST = save_output(updated_PLAST,PRO,FEM,GEOM,QUADRATURE,BC,...
%                 MAT,LOAD,CON,CONSTANT,GLOBAL,PLAST,KINEMATICS);  
        disp(['step = ',sprintf('%d', time_step_counter),'     time = ',... 
          sprintf('%.2e', t_n), ' sec.     dt = ', sprintf('%.2e', dt_old) ,...
          ' sec.'])
  end
  %--------------------------------------------------------------
  time_step_counter = time_step_counter + 1;  
  
  % this is set just to get the removal of old vtu files in output.m
  CON.incrm =  CON.incrm + 1;

end % end on while loop

fprintf(' Normal end of PROGRAM. \n');

toc

delete *RESTART* 
digits(d1);

% delete(gcp)