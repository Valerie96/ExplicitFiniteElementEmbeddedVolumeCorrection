function plot_Eulerian_strain(GEOM,FEM,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,fid3)

space = '   ';


if GEOM.ndime == 2
    fprintf(fid3,'%s%s%s%s<DataArray type="Float32" Name="e (Eulerian)" ',space,space,space,space);
    fprintf(fid3,'NumberOfComponents="4" ');
    fprintf(fid3,'ComponentName0="xx" '); 
    fprintf(fid3,'ComponentName1="xy" ');
    fprintf(fid3,'ComponentName2="yx" ');
    fprintf(fid3,'ComponentName3="yy" ');
    fprintf(fid3,'format="ascii">\n');
elseif GEOM.ndime == 3
    fprintf(fid3,'%s%s%s%s<DataArray type="Float32" Name="e (Eulerian)" ',space,space,space,space);
    fprintf(fid3,'NumberOfComponents="9" ');
    fprintf(fid3,'ComponentName0="xx" ');
    fprintf(fid3,'ComponentName1="xy" ');
    fprintf(fid3,'ComponentName2="xz" ');
    fprintf(fid3,'ComponentName3="yx" ');
    fprintf(fid3,'ComponentName4="yy" ');
    fprintf(fid3,'ComponentName5="yz" ');
    fprintf(fid3,'ComponentName6="zx" ');
    fprintf(fid3,'ComponentName7="zy" ');
    fprintf(fid3,'ComponentName8="zz" ');    
    fprintf(fid3,'format="ascii">\n');
end

for nt = 1:FEM(1).n_elet_type 
for ielement=1:FEM(nt).mesh.nelem  

    %----------------------------------------------------------------------
    % Temporary variables associated with a particular element
    % ready for stress outpue calculation.
    %----------------------------------------------------------------------
    global_nodes    = FEM(nt).mesh.connectivity(:,ielement); 
    material_number = MAT(nt).matno(ielement);               
    matyp           = MAT(nt).matyp(material_number);        
    properties      = MAT(nt).props(:,material_number);      
    xlocal          = GEOM.x(:,global_nodes);            
    x0local         = GEOM.x0(:,global_nodes);               
    Ve              = GEOM.Ve(ielement,nt);                    
    %----------------------------------------------------------------------

    if strcmp(FEM(nt).mesh.element_type,'truss2')
        L       = norm(x0local(:,2) - x0local(:,1));  
        dx      = xlocal(:,2) - xlocal(:,1);        
        l       = norm(dx);                            
        AlStrain = (l^2-L^2)/(2*l^2);
        
        if GEOM.ndime == 2
            fprintf(fid3,'%s%s%s%s%s%.10e %.10e ',space,space,space,space,space,...
                   AlStrain,0);
           fprintf(fid3,'%.10e %.10e\n',0,0);

        elseif GEOM.ndime == 3
           fprintf(fid3,'%s%s%s%s%s%.5e %.5e %.5e ',space,space,space,space,space,...
                   AlStrain,0,0);
           fprintf(fid3,'%.5e %.5e %.5e ', 0,0,0);
           fprintf(fid3,'%.5e %.5e %.5e\n',0,0,0);
        end
        fprintf(fid3,'%s%s%s%s</DataArray>\n',space,space,space,space); 
       return
    end
     
    KINEMATICS(nt) = gradients(xlocal,x0local,FEM(nt).interpolation.element.DN_chi,...
                          QUADRATURE(nt).element,KINEMATICS(nt))  ;     
    %KINEMATICS.F
    F_avg_over_gauss_pts=zeros(GEOM.ndime,GEOM.ndime);
    for igauss=1:QUADRATURE(nt).element.ngauss 
         kinematics_gauss = kinematics_gauss_point(KINEMATICS(nt),igauss);
         F_avg_over_gauss_pts= F_avg_over_gauss_pts+ kinematics_gauss.F;
%          KINEMATICS.F
        
    end
    
    %KINEMATICS.b;
    b_avg = F_avg_over_gauss_pts * F_avg_over_gauss_pts';
    
    e=(1./2.)*(eye(GEOM.ndime)-inv(b_avg));
   
    %----------------------------------------------------------------------
    % Print Eulerian Strain.
    %----------------------------------------------------------------------   
    switch FEM(nt).mesh.element_type
        case'quad4'
           if QUADRATURE(nt).element.ngauss == 4
               fprintf(fid3,'%s%s%s%s%s%.10e %.10e %.10e %.10e \n',space,space,space,space,space,...
                       e(1,1),e(1,2), e(2,1),e(2,2));
           end
        case 'hexa8'
           if QUADRATURE(nt).element.ngauss == 8
                fprintf(fid3,'%s%s%s%s%s%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e  \n',space,space,space,space,space,...
                       e(1,1),e(1,2),e(1,3),e(2,1),e(2,2),e(2,3),e(3,1),e(3,2),e(3,3)   );          
           end
    end
    
 
end
end 
fprintf(fid3,'%s%s%s%s</DataArray>\n',space,space,space,space); 



end

