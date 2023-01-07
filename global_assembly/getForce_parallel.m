%--------------------------------------------------------------------------
% Computes and assemble residual force vector and global tangent stiffness
% matrix except surface (line) element pressure contributions.
%--------------------------------------------------------------------------
function [GLOBAL,updated_PLAST,Jn_1_vec,VolRate_vec,globef_damp,STRESS] = getForce_parallel(xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,INITIAL_KINEMATICS,BC,DAMPING,STRESS,EmbedElt,VolumeCorrect,dt)    
      
Step_globalT_int = zeros(size(GLOBAL.T_int,1),1); 
Step_globalvdamp = zeros(size(GLOBAL.T_int,1),1); 


%If embedded elements are being used, calculated the stresses of the
%embedded elements first so they can be used in the internal force
%calculation of the host elements
if EmbedElt==1
    for ielet=1:FEM(2).mesh.nelem
        global_nodes    = FEM(2).mesh.connectivity(:,ielet);   
        material_number = MAT(2).matno(ielet);
        properties      = MAT(2).props(:,material_number); 
        xlocal          = GEOM.x(:,global_nodes);                     
        x0local         = GEOM.x0(:,global_nodes);  
        [T_internal,~,~,~,Cauchy,LE, ~] = element_force_truss(...
          properties,xlocal,x0local,PLAST(1),GEOM,DAMPING,dt);
       STRESS(2).Cauchy(ielet) = Cauchy;
       STRESS(2).LE(ielet) = LE;
       STRESS(2).InternalForce(ielet,:)=T_internal';
    end
end

% for k = [1:FEM(1).n_elet_type-e]
k=1;


    %--------------------------------------------------------------------------
    % Initialisation of the updated value of the internal variables.
    %--------------------------------------------------------------------------
    updated_PLAST = PLAST(k);
    %--------------------------------------------------------------------------
    % Initialises the total external load vector except pressure contributions.
    %--------------------------------------------------------------------------
    % GLOBAL.external_load = xlamb*GLOBAL.nominal_external_load;

    Jn_1_vec=ones(FEM(1).mesh.nelem,1);
    VolRate_vec=zeros(FEM(1).mesh.nelem,QUADRATURE(1).element.ngauss);                                  

    % Main element loop.
%     ticBytes(gcp);
    %--------------------------------------------------------------------------
    %OG: 144.122 s
    %Parfor 443.804 s
    %Parfor with less broadcast variables 229.719

    connectivity    = FEM(k).mesh.connectivity;   
    matno_list      = MAT(k).matno;     
    matyp_list      = MAT(k).matyp;        
    props_list      = MAT(k).props; 
    xgeom           = GEOM.x;                     
    x0geom          = GEOM.x0;                       
    Ve_list         = GEOM.Ve;
    Tint            = GLOBAL.T_int;
    elt_type        = FEM(k).mesh.element_type;
    DN_Xtemp        = INITIAL_KINEMATICS(1).DN_X;

    temp_Cauchy=zeros(FEM(1).mesh.nelem,6,QUADRATURE(1).element.ngauss);
    temp_LE=zeros(FEM(1).mesh.nelem,6,QUADRATURE(1).element.ngauss);
    QUAD=QUADRATURE(1).element;


    parfor ielet=1:FEM(k).mesh.nelem
    %for ielet=1:FEM(k).mesh.nelem
        temp_globalT_int = zeros(size(Tint,1),1); 
        temp_globalvdamp = zeros(size(Tint,1),1); 
    
            %----------------------------------------------------------------------
            % GATHER Temporary variables associated with a particular element.
            %----------------------------------------------------------------------
            global_nodes    = connectivity(:,ielet);   
            material_number = matno_list(ielet);     
            matyp           = matyp_list(material_number);        
            properties      = props_list(:,material_number); 
            xlocal          = xgeom(:,global_nodes);                     
            x0local         = x0geom(:,global_nodes);                       
            Ve              = Ve_list(ielet);      
%             v_e             = GLOBAL.velocities(global_nodes);
%             a_e             = GLOBAL.accelerations(global_nodes);

            %----------------------------------------------------------------------
            % Select internal variables within the element (plasticity).
            %----------------------------------------------------------------------
            PLAST_element = selecting_internal_variables_element(PLAST(k),matyp,ielet);    
            %----------------------------------------------------------------------
            % Compute internal force and stiffness matrix for an element.
            %----------------------------------------------------------------------    
            switch elt_type
              case 'truss2'
               [T_internal,~,~,~,~,~, ~] = element_force_truss(...
                  properties,xlocal,x0local,PLAST(1),GEOM,DAMPING,dt);
              f_damp=T_internal*0;
              
%                 Jn_1_vec(ielement) = 1;
%                 VolRate_vec(ielement) = 1;
                otherwise
                 
                 DN_X=DN_Xtemp{ielet,1};

                [T_internal,Jn_1,VolRate,f_damp,elt_STRESS] = ...
                InternalForce_explicit(ielet,FEM,xlocal,x0local,...
                Ve,QUAD,properties,CONSTANT,GEOM,PLAST_element,matyp,...
                KINEMATICS,DN_X,MAT,DAMPING,STRESS,VolumeCorrect,dt);


                Jn_1_vec(ielet) = Jn_1;
                VolRate_vec(ielet) = VolRate;
                temp_Cauchy(ielet,:,:) = elt_STRESS.Cauchy(1,:,:);
                temp_LE(ielet,:,:) = elt_STRESS.LE(1,:,:);

            end
%             tocBytes(gcp)
            %----------------------------------------------------------------------
            % Assemble element contribution into global internal force vector.   
            %----------------------------------------------------------------------
%             Step_globalT_int = force_vectors_assembly(T_internal,global_nodes,...
%                            Step_globalT_int,FEM(k).mesh.dof_nodes);
%             Step_globalvdamp = force_vectors_assembly(f_damp,global_nodes,...
%                            Step_globalvdamp,FEM(k).mesh.dof_nodes); 

            global_dofs = FEM(k).mesh.dof_nodes(:,global_nodes);
            temp_globalT_int(global_dofs,1) = T_internal;
            temp_globalvdamp(global_dofs,1) = f_damp;
            
            Step_globalT_int=Step_globalT_int+temp_globalT_int;
            Step_globalvdamp=Step_globalvdamp+temp_globalvdamp;

            %----------------------------------------------------------------------
            % Storage of updated value of the internal variables. 
            %----------------------------------------------------------------------    
%             updated_PLAST = plasticity_storage(PLAST_element,updated_PLAST,matyp,...
%                                                ielet);    

     end
        GLOBAL.T_int = Step_globalT_int;
        STRESS(1).Cauchy = temp_Cauchy;
        STRESS(1).LE     = temp_LE;
%         globef_damp = Step_globalvdamp;

% end
    GLOBAL.T_int = Step_globalT_int;
    globef_damp = Step_globalvdamp;

  if  BC.n_prescribed_displacements > 0
    GLOBAL.external_load_effective(BC.fixdof) = GLOBAL.T_int(BC.fixdof);
    
  end


% algorithm in box 6.1 gives f_n = f_ext - f_int
GLOBAL.Residual              = GLOBAL.external_load - GLOBAL.T_int;
GLOBAL.Reactions(BC.fixdof)  = GLOBAL.T_int(BC.fixdof);

    
end


 


 
