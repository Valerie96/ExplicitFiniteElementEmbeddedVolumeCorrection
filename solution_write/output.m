%--------------------------------------------------------------------------
% TXT file save of converged solution every CON.OUTPUT.incout increments.
% - Coordinates, element connectivity and stress.
% - For node CON.OUTPUT.nwant and dof CON.OUTPUT.iwant output displacement
%   and corresponding force (file name '...FLAGOUT.TXT').
%--------------------------------------------------------------------------
function output(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS,time,dt)
%--------------------------------------------------------------------------
% Restart or write from sratch.
%--------------------------------------------------------------------------
string='a';
if (~PRO.rest && CON.incrm==0)
   string='w';
end
%--------------------------------------------------------------------------
% Open file for single node output.
%--------------------------------------------------------------------------
if (~isempty (CON.OUTPUT.nwant) && CON.OUTPUT.nwant~=0)
   fid_flagout = fopen(PRO.outputfile_name_flagout,string);
end
%--------------------------------------------------------------------------
% Open file for output.
%--------------------------------------------------------------------------
fid = fopen(PRO.outputfile_name,string);
%--------------------------------------------------------------------------
% Print title, load increment number and load factor.
%--------------------------------------------------------------------------
space = '       ';
output_title = [PRO.title space 'at increment:' space  ...
               num2str(CON.incrm) ', ' space 'load:  ' num2str(CON.xlamb) space 'time:  ' num2str(time) space 'dt:  ' num2str(dt)];
fprintf(fid,'%c',output_title);
fprintf(fid,'\n');

%--------------------------------------------------------------------------
% Print number of nodes.
%--------------------------------------------------------------------------
fprintf(fid,'%d',GEOM.npoin);
fprintf(fid,'\n');
fprintf(fid,'Node#, BC, Coordinates, Total Force\n');
fprintf(fid,['     ' repmat('     %s     ',1,GEOM.ndime)],'x', 'y', 'z');
fprintf(fid,['    ' repmat('    %s     ',1,GEOM.ndime) '\n'],'Fx', 'Fy', 'Fz');
%--------------------------------------------------------------------------
% Print boundary codes, coordinates, reactions and external loads.
%--------------------------------------------------------------------------
info                      =  zeros(GEOM.npoin,2 + 2*GEOM.ndime);
info(:,1)                 =  (1:GEOM.npoin)';
info(:,2)                 =  BC.icode;
aux                       =  zeros(FEM(1).mesh.n_dofs,1);
aux(BC.fixdof)            =  GLOBAL.Reactions(BC.fixdof);
aux(BC.freedof)           =  GLOBAL.external_load(BC.freedof);
% aux                       = GLOBAL.external_load - GLOBAL.T_int; %Output Internal Element forces
aux                       =  reshape(aux,GEOM.ndime,[]);
info(:,3:end)             =  [GEOM.x'  aux'];
format                    =  ['%d %d ' repmat('% -1.4E ',1,2*GEOM.ndime) '\n'];
fprintf(fid,format,info');
fprintf(fid,'\n');

fprintf(fid,'Node#, Displacement, Velocity, Acceleration\n');
fprintf(fid,[' ' repmat('     %s     ',1,GEOM.ndime)],'dx', 'dy', 'dz');
fprintf(fid,['    ' repmat('     %s     ',1,GEOM.ndime)],'vx', 'vy', 'vz');
fprintf(fid,['    ' repmat('     %s     ',1,GEOM.ndime) '\n'],'ax', 'ay', 'az');
%--------------------------------------------------------------------------
% Print nodal displacement velocity and acceleration
%--------------------------------------------------------------------------
info                      =  zeros(GEOM.npoin,1 + GEOM.ndime*3);
info(:,1)                 =  (1:GEOM.npoin)';
aux                       =  GLOBAL.velocities;
aux_v                     =  reshape(aux,GEOM.ndime,[]);
aux                       =  GLOBAL.accelerations;
aux_a                     =  reshape(aux,GEOM.ndime,[]);
info(:,2:end)             =  [GEOM.x'-GEOM.x0' aux_v' aux_a'];
format                    =  ['%d ' repmat('% -1.4E ',1,GEOM.ndime) ' || ' repmat('% -1.4E ',1,GEOM.ndime) ' || ' repmat('% -1.4E ',1,GEOM.ndime) '\n'];
fprintf(fid,format,info');

fprintf(fid,'\n');
fprintf(fid,'Element Types: %d \n',FEM(1).n_elet_type);

for i=1:FEM(1).n_elet_type
    %--------------------------------------------------------------------------
    % Print element type. 
    %--------------------------------------------------------------------------
    fprintf(fid,'\n');
    fprintf(fid,'%c',FEM(i).mesh.element_type);
    fprintf(fid,'\n');
    %--------------------------------------------------------------------------
    % Print material type and connectivities.
    % How about don't print connectivities, that's a waste of time
    %--------------------------------------------------------------------------
    fprintf(fid,'Elements: %d',FEM(i).mesh.nelem);
    fprintf(fid,'\n');
    info                      =  zeros(FEM(i).mesh.nelem,2+FEM(i).mesh.n_nodes_elem);
    info(:,1)                 =  (1:FEM(i).mesh.nelem)';
    info(:,2)                 =  MAT(i).matno;
    info(:,3:end)             =  FEM(i).mesh.connectivity';
    switch FEM(i).mesh.element_type
        case 'truss2'
            fprintf(fid,'     %s           %s     \n','Sxx','LExx');    
        case 'hexa8'
        fprintf(fid,[repmat('    %s     ',1,6)], 'Sxx','Sxy','Sxz','Syy','Syz','Szz');
        fprintf(fid,'    ');
        fprintf(fid,[repmat('   %s     ',1,6) '\n'], 'LExx','LExy','LExz','LEyy','LEyz','LEzz');
    end
    %--------------------------------------------------------------------------
    % Print stresses.
    %--------------------------------------------------------------------------
    for ielement=1:FEM(i).mesh.nelem  
        %----------------------------------------------------------------------
        % Temporary variables associated with a particular element
        % ready for stress outpue calculation.
        %----------------------------------------------------------------------
        global_nodes    = FEM(i).mesh.connectivity(:,ielement); 
        material_number = MAT(i).matno(ielement);               
        matyp           = MAT(i).matyp(material_number);        
        properties      = MAT(i).props(:,material_number);      
        xlocal          = GEOM.x(:,global_nodes);            
        x0local         = GEOM.x0(:,global_nodes);               
        Ve              = GEOM.Ve(ielement,i);                 
        %----------------------------------------------------------------------
        % Select internal variables within the element (PLAST).
        %----------------------------------------------------------------------
        PLAST_element = selecting_internal_variables_element(PLAST(i),matyp,ielement);    
        %----------------------------------------------------------------------
        % Compute stresses
        %----------------------------------------------------------------------
        [Stress, LE] = stress_output_from_mem(GEOM.ndime,ielement,matyp,xlocal,x0local,...
                               properties,QUADRATURE(i).element,GEOM,STRESS);
        %----------------------------------------------------------------------
        % Print stress.
        %----------------------------------------------------------------------    
        StressStrain = [Stress; LE];
        format = [repmat('% -1.4E ',1,size(Stress,1)) ' || ' repmat('% -1.4E ',1,size(LE,1)) '\n'];
        fprintf(fid,format',StressStrain);
        
        switch FEM(i).mesh.element_type    
        case 'hexa8'
            xx = Stress(1); xy = Stress(2); xz = Stress(3);
            yy = Stress(4); yz = Stress(5); zz = Stress(6);
            VM = sqrt(0.5*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2) + 3*(xy^2+xz^2+yz^2));
            fprintf(fid,'VM: % -1.4E \n', VM);
        end
    end 
end
%--------------------------------------------------------------------------
% - For node CON.OUTPUT.nwant and dof CON.OUTPUT.iwant output displacement
%   and corresponding force (file name '...FLAGOUT.TXT').
%--------------------------------------------------------------------------
if (~isempty (CON.OUTPUT.nwant) && CON.OUTPUT.nwant~=0)
    increment  = CON.incrm;   
    coordinate = GEOM.x(CON.OUTPUT.iwant,CON.OUTPUT.nwant);
    Force      = GLOBAL.external_load(FEM(1).mesh.dof_nodes(CON.OUTPUT.iwant,CON.OUTPUT.nwant));
    xlamb      = CON.xlamb;
    radius     = CON.ARCLEN.arcln;
    format     = [repmat('% -1.4E ',1,5) '\n'];
    fprintf(fid_flagout,format,[increment,coordinate,Force,xlamb,radius]);
end
fprintf(fid,['\n' repmat('-',1,length(output_title)*2)]);
fprintf(fid,['\n' repmat('-',1,length(output_title)*2) '\n\n']);
fclose(fid);
if (~isempty (CON.OUTPUT.nwant) && CON.OUTPUT.nwant~=0)
   fclose(fid_flagout);
end



