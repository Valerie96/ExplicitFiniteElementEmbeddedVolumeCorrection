function plot_lnV(GEOM,FEM,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS,fid3)

space = '   ';


if GEOM.ndime == 2
    fprintf(fid3,'%s%s%s%s<DataArray type="Float32" Name="LE (lnV)" ',space,space,space,space);
    fprintf(fid3,'NumberOfComponents="4" ');
    fprintf(fid3,'ComponentName0="xx" '); 
    fprintf(fid3,'ComponentName1="xy" ');
    fprintf(fid3,'ComponentName2="yx" ');
    fprintf(fid3,'ComponentName3="yy" ');
    fprintf(fid3,'format="ascii">\n');
elseif GEOM.ndime == 3
    fprintf(fid3,'%s%s%s%s<DataArray type="Float32" Name="LE (lnV)" ',space,space,space,space);
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
        LE      = log(l/L);
        
        if GEOM.ndime == 2
            fprintf(fid3,'%s%s%s%s%s%.10e %.10e ',space,space,space,space,space,...
                   LE,0);
           fprintf(fid3,'%.10e %.10e\n',0,0);

        elseif GEOM.ndime == 3
           fprintf(fid3,'%s%s%s%s%s%.5e %.5e %.5e ',space,space,space,space,space,...
                   LE,0,0);
           fprintf(fid3,'%.5e %.5e %.5e ', 0,0,0);
           fprintf(fid3,'%.5e %.5e %.5e\n',0,0,0);
        end
        continue;
end
        

%     DN_Xtemp=INITIAL_KINEMATICS(nt).DN_X;
%     DN_X=DN_Xtemp{ielement,1};
%     KINEMATICS(nt) = gradients(xlocal,x0local,FEM(nt).interpolation.element.DN_chi,...
%                           QUADRATURE(nt).element,KINEMATICS(nt),DN_X)  ;     
%     %KINEMATICS.F
%     F_avg_over_gauss_pts=zeros(GEOM.ndime,GEOM.ndime);
%     for igauss=1:QUADRATURE(nt).element.ngauss 
%          kinematics_gauss = kinematics_gauss_point(KINEMATICS(nt),igauss);
%          F_avg_over_gauss_pts= F_avg_over_gauss_pts+ kinematics_gauss.F;
%          % KINEMATICS.F
%     end
%     % compute average F (over gauss points)
%     F_avg_over_gauss_pts= F_avg_over_gauss_pts/QUADRATURE(nt).element.ngauss;
%   
%     %KINEMATICS.b;
%     b_avg = F_avg_over_gauss_pts * F_avg_over_gauss_pts';
%    
%     % can be used to check calculation
%     %logm(sqrtm(b_avg));
%     
%     [b_e_vectors,b_e_values] = eig(b_avg);
% 
%     
%     % take care of ln(0) = -Inf 

% switch FEM(nt).mesh.element_type 
%        case 'quad4'
%        if QUADRATURE(nt).element.ngauss == 4 % the 2 is b/c it is the largest e-value from matlab 
%            V=  sqrt(b_e_values(1,1))*b_e_vectors(:,1)*b_e_vectors(:,1)'+ ...
%                sqrt(b_e_values(2,2))*b_e_vectors(:,2)*b_e_vectors(:,2)';
%            
%            lnV=log(sqrt(b_e_values(1,1)))*b_e_vectors(:,1)*b_e_vectors(:,1)'+ ...
%                log(sqrt(b_e_values(2,2)))*b_e_vectors(:,2)*b_e_vectors(:,2)';
%            
%              
%            %lnV=log(V_avg);
%          
% %            
% %            if isinf(lnV(1,2)) == 1
% %                lnV(1,2)=0;
% %            end
% %            if isinf(lnV(2,1)) == 1
% %                lnV(2,1)=0;
% %            end
%        end
%         case 'hexa8'
%        if QUADRATURE(nt).element.ngauss == 8 % the 3 is b/c it is the largest e-value from matlab 
%           V=  sqrt(b_e_values(1,1))*b_e_vectors(:,1)*b_e_vectors(:,1)'+ ...
%               sqrt(b_e_values(2,2))*b_e_vectors(:,2)*b_e_vectors(:,2)'+ ...
%               sqrt(b_e_values(3,3))*b_e_vectors(:,3)*b_e_vectors(:,3)';
%            
%            lnV=log(sqrt(b_e_values(1,1)))*b_e_vectors(:,1)*b_e_vectors(:,1)'+ ...
%                log(sqrt(b_e_values(2,2)))*b_e_vectors(:,2)*b_e_vectors(:,2)'+ ...
%                log(sqrt(b_e_values(3,3)))*b_e_vectors(:,3)*b_e_vectors(:,3)';
%            
% %            if isinf(lnV(1,2)) == 1
% %                lnV(1,2)=0;
% %            end
% %            if isinf(lnV(2,1)) == 1
% %                lnV(2,1)=0;
% %            end
% %            if isinf(lnV(1,3)) == 1
% %                lnV(1,3)=0;
% %            end
% %            if isinf(lnV(3,1)) == 1
% %                lnV(3,1)=0;
% %            end
% %            if isinf(lnV(2,3)) == 1
% %                lnV(2,3)=0;
% %            end
% %            if isinf(lnV(3,2)) == 1
% %                lnV(3,2)=0;
% %            end
%        end

        [~, LE] = stress_output_from_mem(GEOM.ndime,ielement,matyp,xlocal,x0local,...
                       properties,QUADRATURE(nt).element,GEOM,STRESS);
        LE_ave=[mean(LE(1,:)),mean(LE(2,:)),mean(LE(3,:)),mean(LE(4,:)),mean(LE(5,:)),mean(LE(6,:))];
        lnV=[LE_ave(1),LE_ave(2),LE_ave(3);LE_ave(2),LE_ave(4),LE_ave(5);LE(3),LE_ave(5),LE_ave(6)];
% end
    %lnV
%     Abaqus_NE= V-eye(GEOM.ndime);


       
    %----------------------------------------------------------------------
    % Print Logrithmic Strain.
    %----------------------------------------------------------------------   
    switch FEM(nt).mesh.element_type 
        case 'truss2'
            %Assuming this is 3D
        fprintf(fid3,'%s%s%s%s%s%.10e %.10e %.10e %.10e  %.10e %.10e %.10e  %.10e %.10e\n',space,space,space,space,space,...
                   lnVt,0,0,0,0,0);  
          
        case'quad4'
       if QUADRATURE(nt).element.ngauss == 4
           fprintf(fid3,'%s%s%s%s%s%.10e %.10e %.10e %.10e \n',space,space,space,space,space,...
                   lnV(1,1),lnV(1,2), lnV(2,1),lnV(2,2));
           %fprintf(fid3,'%.10e %.10e\n',Green_Strain_Avg(2,1),Green_Strain_Avg(2,2));
       end
        case 'hexa8'
       if QUADRATURE(nt).element.ngauss == 8
           fprintf(fid3,'%s%s%s%s%s%.10e %.10e %.10e %.10e  %.10e %.10e %.10e  %.10e %.10e\n',space,space,space,space,space,...
                   lnV(1,1),lnV(1,2), lnV(1,3),lnV(2,1),lnV(2,2), lnV(2,3),lnV(3,1), lnV(3,2),lnV(3,3)); 

       end
    end

end
end
fprintf(fid3,'%s%s%s%s</DataArray>\n',space,space,space,space); 








end

