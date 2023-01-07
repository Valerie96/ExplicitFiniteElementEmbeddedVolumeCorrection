%--------------------------------------------------------------------------
% Computes and assemble lumped mass matrix. Mass Correct gives the option
% to remove redundant mass from the system when using embedded elements.
%--------------------------------------------------------------------------
function [GLOBAL] = effective_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE,VolumeCorrect)   

%number of dimensions
ndims = GEOM.ndime; 
massSize = ndims*GEOM.npoin;
LumpedMass_eff = zeros(massSize,massSize); %Embedded element mass included in host element mass
LumpedMass_act = zeros(massSize,massSize); %Host elements without embedded element mass


M_total = 0;
effective_M_total = 0;
%Loop over all (host) elements
for ii=1:FEM(1).mesh.nelem
    ielement=ii;
    %----------------------------------------------------------------------
    % Temporary variables associated with a particular element.
    %----------------------------------------------------------------------
    n_nodes_elem    = FEM(1).mesh.n_nodes_elem;
    global_nodes    = FEM(1).mesh.connectivity(:,ielement);   
    material_number = MAT(1).matno(ielement);     
    matyp           = MAT(1).matyp(material_number);        
    properties      = MAT(1).props(:,material_number); % needs density included.
    Ve              = GEOM.Ve(ielement,1);
   
   
    % assign density
    rho = properties(1);
       
    % quadrature locations
    QUADRATURE(1).element.Chi;
    % quadrature weights
    QUADRATURE(1).element.W;
   

    %b: calculate element mass
    me = rho * Ve;
    
    %c: Loop over embedded elements
    mf_tot = 0; mc =0; 
    node_flag = [0 0];
    search = 0;
        for jj=1:FEM(2).mesh.nelem
           jelement = jj; 
           if GEOM.embedded.ElementHost(jelement,1) == ielement 
               node_flag(1) = 1;
               search = search + 1;
           end
           if GEOM.embedded.ElementHost(jelement,2) == ielement 
               node_flag(2) = 1;
               search = search + 1;
           end
           if node_flag(1) + node_flag(2) >= 1
               %i: get element vol and other properties
                    material_numberf = MAT(2).matno(jelement);     
                    properties_f      = MAT(2).props(:,material_numberf); % needs density included.
                    Vf               = GEOM.Ve(jelement,2);
                    rho_f            = properties_f(1);
                    
                    %If only one of the embedded element nodes is in this
                    %host, only part of the total embedded element volume
                    %is associated with that node. Persentage that goes to
                    %mass correction is caluclated by inverse_mapping.m
                    vol_node1 = Vf*node_flag(1)*GEOM.embedded.ElementHost(jelement,3);
                    vol_node2 = Vf*node_flag(2)*GEOM.embedded.ElementHost(jelement,4);
               
                %ii: Calculate embedded element mass
                    mf = rho_f * vol_node1 + rho_f * vol_node2;
                    mf_tot = mf_tot + mf;
                %iii: Calculate correction mass 
                    mc = mc + rho * vol_node1 + rho * vol_node2;
                    
              if search >= GEOM.embedded.HostTotals(ielement,1)
                  break;
              end
                    
           end
            
        end
    
        
        %d: calculate effective mass
        if VolumeCorrect
            meff = me + mf_tot - mc;
            mact = me - mc;
        else
            meff = me + mf_tot;
            mact = me;
        end
        effective_M_total = effective_M_total + meff;
        M_total = M_total + mact + mf_tot;
        
        %e: divide mass equally among nodes
        Meff = (meff/n_nodes_elem);
        Mact = (mact/n_nodes_elem);
        
        %f: form diagonal mass matrix
        Me_eff = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Meff;
        Me_act = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Mact;

        
        %g: transform and scatter scatter mass matrix to global 
        global_dof = FEM(1).mesh.dof_nodes(:,global_nodes);
        LumpedMass_eff(global_dof,global_dof) = LumpedMass_eff(global_dof,global_dof) + Me_eff;
        LumpedMass_act(global_dof,global_dof) = LumpedMass_act(global_dof,global_dof) + Me_act;
        
end % loop on elements


%Loop over all (embedded) elements to add them to the mass matrix
%They need to be in a mass matrix to calculate the total kinetic energy
for jj=1:FEM(2).mesh.nelem
    jelement = jj; 
           %i: get element vol and other properties
            n_nodes_elem_f    = FEM(2).mesh.n_nodes_elem;
            global_nodes_f    = FEM(2).mesh.connectivity(:,jelement);   
            material_numberf = MAT(2).matno(jelement);     
            matyp_f           = MAT(2).matyp(material_numberf);        
            properties_f      = MAT(2).props(:,material_numberf); % needs density included.
            Vf               = GEOM.Ve(jelement,2);
            rho_f            = properties_f(1);
            
            %ii: Calculate embedded element mass
            mf = rho_f * Vf;
            effective_M_total = effective_M_total + mf;
            
                    %e: divide mass equally among embedded nodes
                    Meff = (mf/n_nodes_elem_f);

                    %f: form diagonal mass matrix
                    Identity = eye(n_nodes_elem_f* ndims,n_nodes_elem_f* ndims);
                    Me = Identity * Meff;
                    
                    %g: scatter mass matrix to global 
                    global_dof = FEM(2).mesh.dof_nodes(:,global_nodes_f);
%                     LumpedMass_eff(global_dof,global_dof) = LumpedMass_eff(global_dof,global_dof) + Me;
                    LumpedMass_act(global_dof,global_dof) = LumpedMass_act(global_dof,global_dof) + Me;
end

GLOBAL.M    = LumpedMass_eff;
GLOBAL.M_KE = LumpedMass_act;

M_total = 0 ;
for i = 1:massSize
    M_total = M_total + LumpedMass_act(i,i);
end

end

 


 
