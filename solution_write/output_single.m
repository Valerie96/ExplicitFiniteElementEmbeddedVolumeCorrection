%Single element and node output
function output_single(NodeDisp,NodeLoad,ElementStress,ElementType,PRO,CON,time,dt,fid)
GEOM.ndime = 3;
%--------------------------------------------------------------------------
% Print title, load increment number and load factor.
%--------------------------------------------------------------------------
space = ' ';
output_title = [PRO.title space 'increment:' space  ...
               num2str(CON.incrm) ', ' space 'load:  ' num2str(CON.xlamb) space 'time:  ' num2str(time) space 'dt:  ' num2str(dt)];
% fprintf(fid,'%c',output_title);
fprintf(fid,['increment: ',num2str(CON.incrm) '\n']);
fprintf(fid,['load: ',num2str(CON.xlamb) '\n']);
fprintf(fid,['time: ',num2str(time) '\n']);
fprintf(fid,['dt: ',num2str(dt) '\n']);
fprintf(fid,'\n');

%--------------------------------------------------------------------------
% Print nodes output
%--------------------------------------------------------------------------

% fprintf(fid,'Node#, BC, Coordinates, Total Force\n');
fprintf(fid,'Node Label: %u \n',CON.OUTPUT.nodeout);
format                    =  [' ' repmat(' % -1.4E ',1,GEOM.ndime) '\n'];
fprintf(fid,['  ' repmat('     %s      ',1,GEOM.ndime) '\n'],'x ', 'y ', 'z ');
fprintf(fid,format,NodeLoad(3:5)');

% Print nodal displacement velocity and acceleration
% fprintf(fid,'Node#, Displacement, Velocity, Acceleration\n');
format =  [' ' repmat(' % -1.4E ',1,GEOM.ndime) '\n'];
fprintf(fid,['  ' repmat('     %s      ',1,GEOM.ndime) '\n'],'dx', 'dy', 'dz');
fprintf(fid,format,NodeDisp(2:4)'); 
fprintf(fid,['  ' repmat('     %s      ',1,GEOM.ndime) '\n'],'vx', 'vy', 'vz');
fprintf(fid,format,NodeDisp(5:7)'); 
fprintf(fid,['  ' repmat('     %s      ',1,GEOM.ndime) '\n'],'ax', 'ay', 'az');
fprintf(fid,format,NodeDisp(8:10)');

fprintf(fid,['  ' repmat('     %s      ',1,GEOM.ndime) '\n'],'Fx', 'Fy', 'Fz');
fprintf(fid,format,NodeLoad(6:end)');

%--------------------------------------------------------------------------
% Print element type. 
%--------------------------------------------------------------------------
fprintf(fid,'\n');
fprintf(fid,'Element Type: %s \n',ElementType);
fprintf(fid,'Element Label: %u \n',CON.OUTPUT.eltout(2));
%--------------------------------------------------------------------------
% Print stresses.
%--------------------------------------------------------------------------   
AveStress=mean(ElementStress,2);
switch ElementType
    case 'truss2'
        fprintf(fid,'     %s           \n','Sxx'); 
        format = [repmat('% -1.4E ',1,1) '\n'];
        fprintf(fid,format',AveStress(1)); fprintf(fid,'\n');
        fprintf(fid,'     %s           \n','LExx'); 
        fprintf(fid,format',AveStress(2));fprintf(fid,'\n');
    case 'hexa8'
        fprintf(fid,[repmat('    %s     ',1,6) '\n'], 'Sxx','Sxy','Sxz','Syy','Syz','Szz');
        format = [repmat('% -1.4E ',1,6) '\n'];
        fprintf(fid,format',AveStress(1:6)); fprintf(fid,'\n');
        fprintf(fid,[repmat('   %s     ',1,6) '\n'], 'LExx','LExy','LExz','LEyy','LEyz','LEzz');
        fprintf(fid,format',AveStress(7:12)); fprintf(fid,'\n');

        xx = AveStress(1); xy = AveStress(2); xz = AveStress(3);
        yy = AveStress(4); yz = AveStress(5); zz = AveStress(6);
        VM = sqrt(0.5*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2) + 3*(xy^2+xz^2+yz^2));
        fprintf(fid,'VM: % -1.4E \n', VM);
end

fprintf(fid,['\n' repmat('-',1,75)]);
fprintf(fid,['\n' repmat('-',1,75) '\n\n']);

end
