%Plot a curved shape from Abaqus
%but with a lot of assumptions
%
%6/29/22

%To use, create a curved plate model as you want it in Abaqus with 
%everything you want except for the embedded elements. Dimentions should 
%be in meters

%The plate should be made by creating a rectangle with two curved sides and
%then extruding that shape to make a plate

%A node set of the outside nodes of the part must be definied in the
%assembly.

%An element set called HOST of all the elements in the part must be defined 
%in the assembly

%If you want to define general contact for all interior surfaces (say for a
%model that will include element deletion), create a element surface called
%ALLCONT in the assembly. This should include both interior and exterior
%faces. DO NOT CREATE A GENERAL CONTACT STEP IN ABAQUS In the Matlab user 
%input toggle on ALLSELFCONTACT.

%Write the inp file for the model. Enter the name of the part to be
%embedded with truss elements and the other info in the user input section,
%I think it's pretty self explanitory for now 


%LIMITATIONS
%PART must be oriented so that z+ is the through thickness direction
%Creates alternating layers of 0/90 fibers in the xy plane, stacking 
%direction z+

clear; clc; close all;
tic

OriginalINP="INPfiles/CurvedTest4.inp";
PartName="CurvedTest";
NewINP="CurvedTest4_Fibers.inp";
smallest_elt_length = 0.01;

%Plot part in MATLAB for diganostics
PLOT=true;

%Fiber Material (must be a material definied in Abaqus cae)
Fmat = "FIBERS";

%Create an all with self contact definition
ALLSELFCONTACT=false;
InteractionProperty="IntProp-1";

%Truss Mesh Density (nodes/length)
meshseed = 0.01;
meshseed = smallest_elt_length*1.25;

%Set number of points per row
nptsx = 30;
nptsy = 30;

%Option to enter actual fiber size and number of fibers represented per
%truss and have matlab calculate the rest
%Actual Fiber size
Fact=17E-6;

% Number of fibers per truss element
FpT = 8.6505e+04;

%Calculate FpT from a desired truss diameter
%FpT = (trussD/Fact)^2;
%FpT = (trussD/17E-6)^2;

    %Calculate Truss area
    TArea = FpT*0.25*pi()*Fact^2;
    %Calculate dist between fibers as fiber diameter plus a bit
    TDia = 2*sqrt(TArea/pi())
    DBT = TDia;% + TDia/6;
    DBL = DBT;
PlyThickness=60E-6;
LayersPerPly = PlyThickness/TDia;
%End User Input
%% Code Stuff

% %Read Data from Abaqus INP file
% [Elements,Nodes,Outside]=ReadINP(OriginalINP,PartName,"Outside");
% 
% %Create Surfaces from the part mesh
% [PartSurf]=CreatePartSurface(Elements,Outside);

%Read Data from Abaqus INP file
[BC,CONSTRAINT,ELSET,INSTANCE,NSET,PART,STEP] = ReadAbaqusInputFile(OriginalINP);

[Outside,Nodes,Elements,INSTANCE]=CollectInfo(ELSET,INSTANCE,NSET,PART,STEP,"Outside");

%%
%Find the corners and top point of the curved plate
[Ztop,CornerNodes]=FindBoundsAndCorners(Outside);

%Create Surfaces from the part mesh
[PartSurf]=CreatePartSurface(Elements,Outside);

%Plot curved plate shape, along with identified corner points
figure, hold on, view(3), grid on
patch('Vertices',Nodes(:,:), 'Faces', PartSurf ,'FaceColor','g','FaceAlpha',0.2);
xlabel("x");
ylabel("y");
zlabel("z");
plot3(Ztop(:,2),Ztop(:,3),Ztop(:,4),'kx')
plot3(CornerNodes(:,2),CornerNodes(:,3),CornerNodes(:,4),'bx')
hold off;
%%
%Find the plate arc lengths. This will be the length of the plate in flat
%space Use Corner Nodes 5 and 6, Top node 1
a=CornerNodes(5,2)-CornerNodes(1,2); b=CornerNodes(5,4)-CornerNodes(1,4);
h=sqrt(a^2+b^2);

%The top an bottom of the plate have different arc lengths, so the
%calculations are repeated for each
CordLength_top=CornerNodes(6,2)-CornerNodes(5,2);
CordHight_top=Ztop(1,4)-CornerNodes(6,4);
ArcRadius_top=CordLength_top^2/(8*CordHight_top)+CordHight_top/2;
ArcLength_top=2*ArcRadius_top*asin(CordLength_top/(2*ArcRadius_top));
L=ArcLength_top/2;

CordLength_bottom=CornerNodes(2,2)-CornerNodes(1,2);
CordHight_bottom=Ztop(1,4)-h-CornerNodes(2,4);
ArcRadius_bottom=CordLength_bottom^2/(8*CordHight_bottom)+CordHight_bottom/2;
ArcLength_bottom=2*ArcRadius_bottom*asin(CordLength_bottom/(2*ArcRadius_bottom));
L2=ArcLength_bottom/2;

ArcCenter_top=[Ztop(1,2),Ztop(1,3),Ztop(1,4)-ArcRadius_top];
ArcCenter_bottom=[Ztop(1,2),Ztop(1,3),Ztop(1,4)-h-ArcRadius_bottom];
ArcCenter=ArcCenter_top;

%The "Nodes" used to map from flat to curved space are the four corners of
%the flat plate
Xe=CornerNodes(:,2:4);

%Create the plate in Flat Space
%Set coordinates of flat rectangle node points
w=(CornerNodes(3,3)-CornerNodes(2,3))/2;
x1=[-L,-w,-h/2]; x5=[-L,-w, h/2];
x2=[ L,-w,-h/2]; x6=[ L,-w, h/2];
x3=[ L, w,-h/2]; x7=[ L, w, h/2];
x4=[-L, w,-h/2]; x8=[-L, w, h/2];
xe=[x1;x2;x3;x4;x5;x6;x7;x8];

Surf=[1,2,6,5;2,3,7,6;3,4,8,7;4,1,5,8;4,3,2,1;5,6,7,8];
    FlatSurf=(zeros(size(Surf,1),3));

    for i=1:size(Surf,1)
        j=1+(i-1)*2;
       FlatSurf(j,:)=Surf(i,1:3);
       FlatSurf(j+1,:)=[Surf(i,3), Surf(i,4), Surf(i,1)];
    end


%Add fibers in flat space using a point cloud created in Flat Space
Xbound=[-L,L]; Ybound=[-w,w]; Zbound=[-h/2,h/2];
[PointCloudX,PointCloudY,distx,disty] = CreatePointCloud(Xbound,Ybound,Zbound,DBT,nptsx,nptsy);

%Plot the mesh shape and point cloud in Flat Space
PlotMeshShape(xe,Outside,FlatSurf,PointCloudX,PointCloudY,PLOT);

%Determine Point Cloud Points that are outside the part
    XIN=inpolyhedron(FlatSurf,xe,PointCloudX);
    YIN=inpolyhedron(FlatSurf,xe,PointCloudY);

% Identify End points of nodes
[EndpointsX,EndpointsY]=FindTrussEndpoints(PointCloudX,PointCloudY,XIN,YIN,distx,disty);

% Plot Truss Points (in and out) and Endpoints, Plot Truss lines inside Part
PlotTrussPoints(xe,FlatSurf,EndpointsX,EndpointsY,PLOT)

% Mesh the truss elements
[NodesX,NodesY,totalNodes,numXnodes,numYnodes]=MeshTrussElements(Xbound,Ybound,EndpointsX,EndpointsY,meshseed);

%% Map Flat Space Truss Nodes to Curved Space
CurvedNodesX=NodesX;
CurvedNodesY=NodesY;


for i=1:size(CurvedNodesX,1)
    for j=1:3:size(CurvedNodesX,2)
        N=ShapeFunctions(NodesX(i,j:j+2),L,w,h);
        CurvedNodesX(i,j:j+2)=N*Xe;

        rz=ArcRadius_top*(1/2+NodesX(i,j+2)/h)-ArcRadius_bottom*(NodesX(i,j+2)/h-1/2);
        z1=sqrt(rz^2-(CurvedNodesX(i,j)-ArcCenter(1))^2)+ArcCenter(3);
        CurvedNodesX(i,j+2)=z1;
    end
end

for i=1:size(CurvedNodesY,1)
    for j=1:3:size(CurvedNodesY,2)
        N=ShapeFunctions(NodesY(i,j:j+2),L,w,h);
        CurvedNodesY(i,j:j+2)=N*Xe;

        rz=ArcRadius_top*(1/2+NodesY(i,j+2)/h)-ArcRadius_bottom*(NodesY(i,j+2)/h-1/2);
        z1=sqrt(rz^2-(CurvedNodesY(i,j)-ArcCenter(1))^2)+ArcCenter(3);
        CurvedNodesY(i,j+2)=z1;
    end
end

lx=size(CurvedNodesX,2);
ly=size(CurvedNodesY,2);

% Plot mesh for visulization of density
figure, hold on, view(3), grid on
lightGrey1   = [0.85 0.85 0.85];
patch('Vertices',Nodes(:,:), 'Faces', PartSurf ,'FaceColor',lightGrey1,'FaceAlpha',0.2);
plot3(CurvedNodesX(:,[1:3:lx]),CurvedNodesX(:,[2:3:lx]),CurvedNodesX(:,[3:3:lx]),'bo');
plot3(CurvedNodesY(:,[1:3:ly]),CurvedNodesY(:,[2:3:ly]),CurvedNodesY(:,[3:3:ly]),'go');
for i=1:size(CurvedNodesX,1)
    for j=1:3:(size(CurvedNodesX,2)-3)
        plot3([CurvedNodesX(i,j) CurvedNodesX(i,j+3)]', [CurvedNodesX(i,j+1) CurvedNodesX(i,j+4)]', [CurvedNodesX(i,j+2) CurvedNodesX(i,j+5)]','b');
    end
end
for i=1:size(CurvedNodesY,1)
    for j=1:3:(size(CurvedNodesY,2)-3)
        plot3([CurvedNodesY(i,j) CurvedNodesY(i,j+3)]', [CurvedNodesY(i,j+1) CurvedNodesY(i,j+4)]', [CurvedNodesY(i,j+2) CurvedNodesY(i,j+5)]','g');
    end
end
% zlim([-0.04,0.6]);
title("Truss Element Mesh Density");
xlabel("x");
ylabel("y");
zlabel("z");

% Write New INP File
[numXelts,numYelts]=WriteInputFile(OriginalINP,NewINP,Fmat,TArea,ALLSELFCONTACT,CurvedNodesX,CurvedNodesY,numXnodes,numYnodes);

toc
fprintf("Total number of Truss nodes: %u\n", totalNodes);
fprintf("Total number of Truss elements: %u\n", numXelts+numYelts);
fprintf("Total number of nodes: %u\n", totalNodes+size(Nodes,1));
fprintf("Total number of elements: %u\n\n", numXelts+numYelts+size(Elements,1));

%% Functions
function [Elements,Nodes,Outside]=ReadINP(OriginalINP,PartName,OutsideNSetName)
%Open text file and read data  
%Check if file exists and opens sucessfully
   if ~isfile(OriginalINP)
       fprintf("%s Not Found\n", OriginalINP);
   end
    fid = fopen(OriginalINP,'r');
        if fid < 0
            %disp('Error')
            fprintf("%s Not Found\n", OriginalINP);
            
        else
            
        %Look for part mesh
        tline = fgetl(fid);
            while ~contains(tline, strcat("*Part, name=",PartName))
                tline = fgetl(fid);
            end
            tline = fgetl(fid);

        %Read in node coordinates
            Nodes = fscanf(fid, '%f %*c %f %*c %f %*c %f ',[4, inf]);
            Nodes = Nodes';
            
        %Read in Element Type
            tline = fgetl(fid);
            EltType = tline(1,16:end);

        %Read in element definition 
            Elements = fscanf(fid, '%u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u',[9, inf]);
            Elements = Elements';     
  
        %Find Instance Definition
            while ~contains(tline, strcat("*Instance, name=",PartName))
                tline = fgetl(fid);
            end
            tline = fgetl(fid);
                if contains(tline, "*End Instance")
                    translate = [0 0 0];
                    rotvec = [0 0 0; 0 0 0];
                    rotangle = 0;
                else
                    translate = sscanf(tline, '%f %*c %f %*c %f ',[3, 1])';
                    tline = fgetl(fid);
                    if contains(tline, "*End Instance")
                        rotvec = [0 0 0; 0 0 0];
                        rotangle = 0;
                    else
                        rotvec = sscanf(tline, '%f %*c %f %*c %f %*c %f %*c %f %*c %f %*c %*f',[3, 2])';
                        rotangle = sscanf(tline,'%*f %*c %*f %*c %*f %*c %*f %*c %*f %*c %*f %*c %f');
                        tline = fgetl(fid);
                    end
                end
                        
        
        %Translate and Rotate the instance to global position
        if rotangle == 0
            R = eye(3);
        else
            u = rotvec(2,:) - rotvec(1,:);
            unit = u/abs(u);
            c = cosd(rotangle); s = sind(rotangle);
            R = [c+u(1)^2*(1-c), u(1)*u(2)*(1-c)-u(3)*s, u(1)*u(3)*(1-c)+u(2)*s; u(1)*u(2)*(1-c)+u(3)*s, c+u(2)^2*(1-c), u(2)*u(3)*(1-c)-u(1)*s; u(1)*u(3)*(1-c)-u(2)*s, u(2)*u(3)*(1-c)+u(1)*s, c+u(3)^2*(1-c)];
        end

        for j=1:size(Nodes,1)
            Nodes(j,2:4) = Nodes(j,2:4) + translate;
            x = Nodes(j,2:4) - rotvec(1,:);
            x = R*x';
            Nodes(j,2:4) = x' + rotvec(1,:);
        end
            
            %Find Outside Surface Node Set
            while ~contains(tline, OutsideNSetName)
                tline = fgetl(fid);
            end

         %Read in Outside Nodes
         
             %This is tricky because the delimiting character in the inp is a
             %space and a comma instead of a line break. But there could be
             %several line breaks in the information anyaway

             %Read in all of the node numbers into the character array Nums.
             %Read number characters only
                tline = fgetl(fid);
                Nums=NaN;
                while ~contains(tline, "*")
                   Nums = [Nums regexp(tline, '\d*', 'Match')];
                   tline = fgetl(fid);
                end

             %Cast the character array Nums as an integer array Outside   
                for ii = 2:length(Nums)
                    if ~isempty(Nums{ii})
                        Outside(ii-1,1) = str2double(Nums{ii});
                    else
                        Outside(ii-1,1) = NaN;
                    end
                end

                %If the outside node set is defined by an itteration
                if length(Outside)<5 && Outside(end,1)==1
                    Outside = [Outside(1,1):Outside(2,1)]';
                end

                %Add node coordinates to the Outside Node set
                Outside = [Outside zeros(length(Outside), 3)];
                for jj = 1:length(Outside)
                    Outside(jj,:)=Nodes(Outside(jj,1),:);
                end
            %Close the inp file
            fclose(fid);

        end
end

function [Outside,GlobalNodes,GlobalElements,INSTANCE]=CollectInfo(ELSET,INSTANCE,NSET,PART,STEP,OutsideNSet)
    GlobalNodes=zeros(1,3); GlobalElements=zeros(1,8);
    %Translate Instances to Global Positions
    n_Instances = length(INSTANCE);
    
    %count total nodes and total elements
    %assuming there are no more than 2 element types and they are C3D8 or T3D2 
    n_nodes = 0;
    n_elements = zeros(2,1);
    n_Parts = length(PART);
    for i=1:n_Instances
        for j=1:n_Parts
            if strcmp(INSTANCE(i).part, PART(j).name)
                %Assign global node numbers to the instance
                n_nodes_new = n_nodes + size(PART(j).nodes,1);
                global_nodes = [n_nodes+1:n_nodes_new]; 
                INSTANCE(i).global_nodes = [global_nodes' PART(j).nodes];
                INSTANCE(i).global_elements = PART(j).elements + ones(size(PART(j).elements,1),size(PART(j).elements,2))*n_nodes(end);
                n_nodes_old=n_nodes;
                n_nodes = n_nodes_new;
    
                INSTANCE(i).eltType = PART(j).eltType;
                if contains(PART(j).eltType, "C3D8")
                    n_old = n_elements(1);
                    n_elements(1) = n_elements(1) + size(PART(j).elements,1);
                    GlobalElements(n_old+1:n_elements(1),:)=INSTANCE(i).global_elements;
                elseif contains(PART(j).eltType,"T3D2")
                    n_elements(2) = n_elements(2) + size(PART(j).elements,1);
                end
                break;
            end
        end
    
            %Translate and Rotate the instance to global position
            if INSTANCE(i).rotangle == 0
                R = eye(3);
            else
                u = INSTANCE(i).rotvec(2,:) - INSTANCE(i).rotvec(1,:);
                c = cosd(INSTANCE(i).rotangle); s = sind(INSTANCE(i).rotangle);
                R = [c+u(1)^2*(1-c), u(1)*u(2)*(1-c)-u(3)*s, u(1)*u(3)*(1-c)+u(2)*s; u(1)*u(2)*(1-c)+u(3)*s, c+u(2)^2*(1-c), u(2)*u(3)*(1-c)-u(1)*s; u(1)*u(3)*(1-c)-u(2)*s, u(2)*u(3)*(1-c)+u(1)*s, c+u(3)^2*(1-c)];
            end
    
            for j=1:size(INSTANCE(i).global_nodes,1)
                INSTANCE(i).global_nodes(j,2:4) = INSTANCE(i).global_nodes(j,2:4) + INSTANCE(i).translate;
                x = INSTANCE(i).global_nodes(j,2:4) - INSTANCE(i).rotvec(1,:);
                x = R*x';
                INSTANCE(i).global_nodes(j,2:4) = x' + INSTANCE(i).rotvec(1,:);
            end

            GlobalNodes(n_nodes_old+1:n_nodes_new,:) = INSTANCE(i).global_nodes(:,2:4);
    end

    %Find Outside set
    Outside=zeros(1,4); row=0;
    for i=1:length(NSET)
        if strcmp(NSET(i).name, OutsideNSet)
            for j=1:n_Instances
                if strcmp(INSTANCE(j).name, NSET(i).instance)
                    len=row+size(NSET(i).nodes,2);
                    Outside(row+1:len,:) = INSTANCE(j).global_nodes(NSET(i).nodes,:);
                    row=len;
                end
            end
        end
        
    end
    



end

function [Ztop,CornerNodes]=FindBoundsAndCorners(Outside)

%Find the 8 corners of the plate
%first find part bounds
    Xbound = [0, 0];
    Ybound = [0, 0];
    Zbound = [0, 0];

    for i=1:length(Outside(:,1))
        if Outside(i,2)>Xbound(2)
            Xbound(2) = Outside(i,2);
        elseif Outside(i,2)<Xbound(1)
            Xbound(1) = Outside(i,2);
        end

        if Outside(i,3)>Ybound(2)
            Ybound(2) = Outside(i,3);
        elseif Outside(i,3)<Ybound(1)
            Ybound(1) = Outside(i,3);
        end

        if Outside(i,4)>Zbound(2)
            Zbound(2) = Outside(i,4);
        elseif Outside(i,4)<Zbound(1)
            Zbound(1) = Outside(i,4);
        end
    end
    
    %Find the nodes that touch these bounds
    %Right now this will only work for this one specific orientation of the
    %part with respect to the global axis, but this is proof of concept so
    %it's ok for now
    Ztop=zeros(2,4);
    CornerNodes=zeros(8,4);
    xave=(Xbound(1)+Xbound(2))/2;
    for i=1:length(Outside(:,1))

        if Outside(i,2)==Xbound(2)
            if Outside(i,3)==Ybound(2)
                CornerNodes(7,:)=Outside(i,:);
            elseif Outside(i,3)==Ybound(1)
                CornerNodes(6,:)=Outside(i,:);
            end

        elseif Outside(i,2)==Xbound(1)
            if Outside(i,3)==Ybound(2)
                CornerNodes(8,:)=Outside(i,:);
            elseif Outside(i,3)==Ybound(1)
                CornerNodes(5,:)=Outside(i,:);
            end
        end

        if Outside(i,4)==Zbound(2)
            if Outside(i,3)==Ybound(2)
                Ztop(1,:)=Outside(i,:);
            elseif Outside(i,3)==Ybound(1)
                Ztop(2,:)=Outside(i,:);
            end

        elseif Outside(i,4)==Zbound(1)
            if Outside(i,3)==Ybound(2)
                if Outside(i,2)>xave
                    CornerNodes(3,:)=Outside(i,:);
                else
                    CornerNodes(4,:)=Outside(i,:);
                end
            elseif Outside(i,3)==Ybound(1)
                if Outside(i,2)>xave
                    CornerNodes(2,:)=Outside(i,:);
                else
                    CornerNodes(1,:)=Outside(i,:);
                end
            end
        end
    end
end

function [Face3]=CreatePartSurface(Elements,Outside)
% Rearrange node connectivies to create faces instead of elements 
%For C3D8, each element has 6 surfaces, but not all are unique
%But who cares we're only interested in the outside surfaces, which we know
%because of the outside node set
 
    Surfaces = zeros(size(Elements,1)*6,4);
    for i=1:size(Elements,1)
        j=1+(i-1)*6;
        Surfaces(j,:)=[Elements(i,4), Elements(i,3), Elements(i,2), Elements(i,1)];
        Surfaces(j+1,:)=[Elements(i,1), Elements(i,5), Elements(i,8), Elements(i,4)];
        Surfaces(j+2,:)=[Elements(i,1), Elements(i,2), Elements(i,6), Elements(i,5)];
        Surfaces(j+3,:)=Elements(i,5:8);
        Surfaces(j+4,:)=[Elements(i,2), Elements(i,3), Elements(i,7), Elements(i,6)];
        Surfaces(j+5,:)=[Elements(i,3), Elements(i,4), Elements(i,8), Elements(i,7)];
    end

    %Find the surfaces that contain only outside nodes and add them to
    %OutsideSurfaces
    
    %Loop through Surfaces. For each suface (i) check if every element (j)
    %is also in Outside. If one element is not part of Outside, the surface
    %is not an outside face
        OutSurface=logical(zeros(size(Surfaces,1),1));
        for i=1:size(Surfaces,1)
           for j=1:4
               Out=false;
               for k=1:size(Outside,1)
                  if Surfaces(i,j)==Outside(k,1)
                      Out=true;
                  end               
               end

               if Out==false
                   break;
               end           
           end
           OutSurface(i,1)=Out;
        end
    
    %Create a smaller matrix of only outside surfaces (faces)
        Face4=Surfaces(OutSurface(:,1)',:);

    %Divide Faces into triangles
        Face3=(zeros(size(Face4,1),3));
    
        for i=1:size(Face4,1)
            j=1+(i-1)*2;
           Face3(j,:)=Face4(i,1:3);
           Face3(j+1,:)=[Face4(i,3), Face4(i,4), Face4(i,1)];
        end

end

function [N]=ShapeFunctions(x,l,w,h)
    t=h/2; s=1/(8*l*w*t); 
    N1 = -s*(x(1)-l)*(x(2)-w)*(x(3)-t);
    N2 =  s*(x(1)+l)*(x(2)-w)*(x(3)-t);
    N3 = -s*(x(1)+l)*(x(2)+w)*(x(3)-t);
    N4 =  s*(x(1)-l)*(x(2)+w)*(x(3)-t);
    N5 =  s*(x(1)-l)*(x(2)-w)*(x(3)+t);
    N6 = -s*(x(1)+l)*(x(2)-w)*(x(3)+t);
    N7 =  s*(x(1)+l)*(x(2)+w)*(x(3)+t);
    N8 = -s*(x(1)-l)*(x(2)+w)*(x(3)+t);
    N=[N1 N2 N3 N4 N5 N6 N7 N8];
end

function [TrussXX,TrussYY,distx,disty] = CreatePointCloud(Xbound,Ybound,Zbound,DBT,nptsx,nptsy)
   
%Create Point Cloud for Fiber Layers based on the X and Y boundaries

%Trusses with their axies aligned perpedicular to the x axis (ie each truss
%can be defined by it's x coordinate) are denoted by X
%Trusses with their axies aligned perpedicular to the y axis (ie each truss
%can be defined by it's y coordinate) are denoted by Y
  
%Calculate the number of trusses in a direction
    ntrussx = (Xbound(2)-Xbound(1))/DBT;
    ntrussy = (Ybound(2)-Ybound(1))/DBT;
    ntrussz = (Zbound(2)-Zbound(1))/DBT;
    
%Round to nearest int
    ntrussx=uint32(ntrussx);
    ntrussy=uint32(ntrussy);
    ntrussz=uint32(ntrussz);

%Create Matrix of truss nodes
    TrussX = zeros(ntrussx*nptsx,3);
    TrussY = zeros(ntrussy*nptsy,3);
    
%Define the distance between points on a fiber
    disty = (Ybound(2)-Ybound(1))/(nptsy-1);
    distx = (Xbound(2)-Xbound(1))/(nptsx-1);
%First trusses are not quite located on the edge of the bounds. First layer
%is not quite at the lowest z value
    xt = Xbound(1)+DBT/4;
    yt = Ybound(1)+DBT/4;    
    z1 = Zbound(1)+(DBT/4);
    for i=1:nptsx:size(TrussX,1)
        ypt=Ybound(1);
        for j=0:nptsx-1
            TrussX(i+j,:) = [xt, ypt, z1];
            ypt=ypt+disty;
        end
        xt=xt+DBT;
    end
    for i=1:nptsy:size(TrussY,1)
        xpt=Xbound(1);
            for j=0:nptsy-1
                TrussY(i+j,:) = [xpt, yt, z1];
                xpt=xpt+distx;
            end
        yt=yt+DBT;
    end
    
%Copy Truss Matricies into 3D matricies with different z coordinates
%Create matricies for each truss layer
    ntrussz2 = uint32(ntrussz/2+1);
    TrussXX = zeros(uint32(size(TrussX,1))*ntrussz2,3);
    TrussYY = zeros(uint32(size(TrussY,1))*ntrussz2,3);

    %For large ntrussz, TrussXX gets too large for int16 indexes to handle
    %(matrix dimentions are represented as int16 which only contains numbers
    %from -32768 to 32767 (2^15))
    
%Loop through larger truss matricies, copy the x y values from the smaller
%matricies and then add z coordinates

    for i=1:ntrussz2
        kx=1+(i-1)*size(TrussX,1);
        ky=1+(i-1)*size(TrussY,1);
        TrusXXCheck=TrussXX(kx:kx+size(TrussX,1)-1,:);
        TrussXX(kx:kx+size(TrussX,1)-1,:)=TrussX;
        TrussYY(ky:ky+size(TrussY,1)-1,:)=TrussY;

        %Change z values
        TrussXX(kx:kx+size(TrussX,1)-1,3)=ones(1,size(TrussX,1))*(z1+double(i-1)*2*DBT);
        TrussYY(ky:ky+size(TrussY,1)-1,3)=ones(1,size(TrussY,1))*(z1+double(1+(i-1)*2)*DBT);
    end

end

function PlotMeshShape(Nodes,Outside,PartSurf,PointCloudX,PointCloudY,PLOT)
    if PLOT
        figure, hold on, view(3), grid on
        patch('Vertices',Nodes, 'Faces', PartSurf ,'FaceColor','g','FaceAlpha',0.2);
        for i=1:size(PointCloudX,1)
            plot3(PointCloudX(i,1), PointCloudX(i,2), PointCloudX(i,3), 'b.');
        end     
        for i=1:size(PointCloudY,1)
            plot3(PointCloudY(i,1), PointCloudY(i,2), PointCloudY(i,3), 'mo');
        end
        title("Point Cloud");
        xlabel("x");
        ylabel("y");
        zlabel("z");
%         zlim([0 0.005]);
        hold off

%         figure, hold on, view(3), grid on
%             for i=1:size(Outside,1)
%                 plot3(Outside(i,2), Outside(i,3), Outside(i,4), 'k.');
%             end
%         title("Outside Suface Nodes");
%         hold off
    end
end

function [EndpointsX,EndpointsY]=FindTrussEndpoints(PointCloudX,PointCloudY,XIN,YIN,distx,disty)
    TrussXin=PointCloudX(XIN,:);
    TrussYin=PointCloudY(YIN,:);

%Create a matrix for the trussX endpoints. Endpoint 1,2 [ x1 y1 z1 x2 y2 z2]
    EndpointsX=zeros(size(TrussXin,1),6);
    nEndpt=1;
    EndpointsX(nEndpt,1:3)=TrussXin(1,1:3);
    for i=1:size(TrussXin,1)-1
       %Check if we're still on the same X coordinate
       if TrussXin(i,1)~=TrussXin(i+1,1)
           EndpointsX(nEndpt,4:6)=TrussXin(i,1:3); %End of this truss
           EndpointsX(nEndpt+1,1:3)=TrussXin(i+1,1:3); %Next pt is start of next truss
           nEndpt=nEndpt+1;
       else
           %If we are still on the same X coordinate, check if the y coordinates
           %are more than one disty away from each other
           if TrussXin(i+1,2) > (TrussXin(i,2)+disty)
           EndpointsX(nEndpt,4:6)=TrussXin(i,1:3);
           EndpointsX(nEndpt+1,1:3)=TrussXin(i+1,1:3);
               nEndpt=nEndpt+1;
           end
       end
    end
    EndpointsX(nEndpt,4:6)=TrussXin(end,1:3);
    EndpointsX=EndpointsX(1:nEndpt,:);

%Create a matrix for the trussY endpoints. Endpoint 1,2 [ x1 y1 z1 x2 y2 z2]
    EndpointsY=zeros(size(TrussYin,1),6);
    nEndpt=1;
    EndpointsY(nEndpt,1:3)=TrussYin(1,1:3);
    for i=1:size(TrussYin,1)-1
       %Check if we're still on the same Y coordinate
       if TrussYin(i,2)~=TrussYin(i+1,2)
           EndpointsY(nEndpt,4:6)=TrussYin(i,1:3); %End of this truss
           EndpointsY(nEndpt+1,1:3)=TrussYin(i+1,1:3); %Next pt is start of next truss
           nEndpt=nEndpt+1;
       else
           %If we are still on the same Y coordinate, check if the x coordinates
           %are more than one distx away from each other
           if TrussYin(i+1,1) > (TrussYin(i,1)+distx)
           EndpointsY(nEndpt,4:6)=TrussYin(i,1:3);
           EndpointsY(nEndpt+1,1:3)=TrussYin(i+1,1:3);
               nEndpt=nEndpt+1;
           end
       end
    end
    EndpointsY(nEndpt,4:6)=TrussYin(end,1:3);
    EndpointsY=EndpointsY(1:nEndpt,:);
end

function PlotTrussPoints(Nodes,PartSurf,EndpointsX,EndpointsY,PLOT)
    if PLOT
       figure, hold on, view(3), grid on
       patch('Vertices',Nodes, 'Faces', PartSurf ,'FaceColor','k','FaceAlpha',0.1);
            plot3(EndpointsX(:,1), EndpointsX(:,2), EndpointsX(:,3), 'bs');
            plot3(EndpointsX(:,4), EndpointsX(:,5), EndpointsX(:,6), 'bs');
            plot3(EndpointsY(:,1), EndpointsY(:,2), EndpointsY(:,3), 'gs');
            plot3(EndpointsY(:,4), EndpointsY(:,5), EndpointsY(:,6), 'gs');
            title("Fiber Endpoints");
            xlabel("x");
            ylabel("y");
            zlabel("z");
        hold off;
    
       figure, hold on, view(3), grid on
            plot3([EndpointsX(:,1) EndpointsX(:,4)]', [EndpointsX(:,2) EndpointsX(:,5)]', [EndpointsX(:,3) EndpointsX(:,6)]','b','LineWidth',3);
            plot3([EndpointsY(:,1) EndpointsY(:,4)]', [EndpointsY(:,2) EndpointsY(:,5)]', [EndpointsY(:,3) EndpointsY(:,6)]','g','LineWidth',3);
            title("Fibers")
            xlabel("x");
            ylabel("y");
            zlabel("z");
    %         zlim([-20E-3 20E-3]);
    %         zlim([0 0.05]);
    %         ylim([0 0.05]);
    %         xlim([0 0.05]);
        hold off;
    %     
    end
end

function [NodesX,NodesY,totalNodes,numXnodes,numYnodes]=MeshTrussElements(Xbound,Ybound,EndpointsX,EndpointsY,meshseed)
    %Create matricies to hold the truss mesh nodes. Each row represents a truss
    
    %The maximum number of nodes in any truss element will be the X (or Y)
    %bounds divided by the meshseed size -1
        Nxmax=int16(((Ybound(2)-Ybound(1))/meshseed))+2;
        Nymax=int16(((Xbound(2)-Xbound(1))/meshseed))+2;
    
        NodesX=zeros(size(EndpointsX,1),3*Nxmax);
        NodesY=zeros(size(EndpointsY,1),3*Nymax);
    
    
    %Loop through the list of endpoints and add them and the other mesh nodes
    %to NodesX
        numXnodes=0;
        for i=1:size(EndpointsX,1)
           tleng=EndpointsX(i,5)-EndpointsX(i,2);
           midnodes=int16(tleng/meshseed)-1;
           nodespace=tleng/double(midnodes+1);
    
           NodesX(i,1:3)=EndpointsX(i,1:3);
           numXnodes=numXnodes+1;
           ycord = EndpointsX(i,2);
    
           j=4;
           while j<(3+midnodes*3)
              ycord=ycord+nodespace; 
              NodesX(i,j:j+2)=[EndpointsX(i,1) ycord  EndpointsX(i,3)];
              j=j+3;
              numXnodes=numXnodes+1;
           end
    
           NodesX(i,j:j+2)=EndpointsX(i,4:6);
           numXnodes=numXnodes+1;
           %Assign any remaining places to NaN
           len_end = size(NodesX(i,j+3:end));
           NodesX(i,j+3:end)=NaN(len_end);
        end
    
    %Loop through the list of endpoints and add them and the other mesh nodes
    %to NodesY
        numYnodes=0;
        for i=1:size(EndpointsY,1)
           tleng=EndpointsY(i,4)-EndpointsY(i,1);
           midnodes=int16(tleng/meshseed)-1;
           nodespace=tleng/double(midnodes+1);
    
           NodesY(i,1:3)=EndpointsY(i,1:3);
           numYnodes=numYnodes+1;
           xcord = EndpointsY(i,1);
    
           j=4;
           while j<(3+midnodes*3)
              xcord=xcord+nodespace; 
              NodesY(i,j:j+2)=[xcord EndpointsY(i,2) EndpointsY(i,3)];
              numYnodes=numYnodes+1;
              j=j+3;
           end
    
           NodesY(i,j:j+2)=EndpointsY(i,4:6);
           numYnodes=numYnodes+1;
           len_end = size(NodesY(i,j+3:end));
           NodesY(i,j+3:end)=NaN(len_end);
        end
    
        totalNodes=numXnodes+numYnodes;

%         NodesX = NodesX(:,all(~isnan(NodesX)));
        lx=size(NodesX,2);
%         NodesY = NodesY(:,all(~isnan(NodesY)));
        ly=size(NodesY,2);

    % Plot mesh for visulization of density

    figure, hold on, view(3), grid on
    plot3([EndpointsX(:,1) EndpointsX(:,4)]', [EndpointsX(:,2) EndpointsX(:,5)]', [EndpointsX(:,3) EndpointsX(:,6)]','b');
    plot3([EndpointsY(:,1) EndpointsY(:,4)]', [EndpointsY(:,2) EndpointsY(:,5)]', [EndpointsY(:,3) EndpointsY(:,6)]','g');
    plot3(NodesX(:,[1:3:lx]),NodesX(:,[2:3:lx]),NodesX(:,[3:3:lx]),'bo');
    plot3(NodesY(:,[1:3:ly]),NodesY(:,[2:3:ly]),NodesY(:,[3:3:ly]),'go');
    title("Truss Element Mesh Density");
    xlabel("x");
    ylabel("y");
    zlabel("z");
%     zlim([-20E-3 20E-3]);
end

function [numXelts,numYelts]=WriteInputFile(OriginalINP,NewINP,Fmat,TArea,ALLSELFCONTACT,NodesX,NodesY,numXnodes,numYnodes)
    %Open the origianl inp for reading, create a new inp for writing
   if isfile(OriginalINP)
    fidOld = fopen(OriginalINP,'r');
        if fidOld < 0
            disp('Error')
        else
    fidNew = fopen(NewINP, 'w');
    
    %Replicate the original inp down to the end of the first part
    tline = fgets(fidOld);
    while ~contains(tline,"*End Part")
        
        fwrite(fidNew,tline);
        tline = fgets(fidOld);
    end
        fwrite(fidNew,tline);
    
    %Insert Truss elements as a new part
        fprintf(fidNew, "**\n");
        fprintf(fidNew, "*Part, name=Fibers0\n");
        fprintf(fidNew, "*Node\n");
        nodenum=0;
        for i = 1:size(NodesX,1)
            j=1; 
            while ~isnan(NodesX(i,j))
            nodenum=nodenum+1;    
            fprintf(fidNew, "      %u, %13f, %13f, %13f\n",nodenum, NodesX(i,j), NodesX(i,j+1), NodesX(i,j+2));
            j=j+3;
            end
        end

        fprintf(fidNew, "*Element, type=T3D2\n");
        eltnum=1;
        nodenum=1;
        for i = 1:size(NodesX,1)
           j=4;
           while ~isnan(NodesX(i,j))
               fprintf(fidNew,"%u, %u, %u\n", eltnum, nodenum, nodenum+1);
               eltnum=eltnum+1;
               nodenum=nodenum+1;
               j=j+3;
           end
           nodenum=nodenum+1;
        end
        numXelts=eltnum;
        
        fprintf(fidNew,"*Nset, nset=Fib0, internal, generate\n");
        fprintf(fidNew,"1, %u, 1\n", numXnodes);
        fprintf(fidNew,"*Elset, elset=Fib0, internal, generate\n");
        fprintf(fidNew,"1, %u, 1\n", numXelts);
        fprintf(fidNew, "**Section: TrussSection\n");
        fprintf(fidNew,"*Solid Section, elset=Fib0, material=%s\n%E,\n", Fmat,TArea);
        fprintf(fidNew, "*End Part\n");
        
 
        fprintf(fidNew, "**\n");
        fprintf(fidNew, "*Part, name=Fibers90\n");
        fprintf(fidNew, "*Node\n");
        nodenum=0;
        for i = 1:size(NodesY,1)
            j=1; 
            while ~isnan(NodesY(i,j))
            nodenum=nodenum+1;    
            fprintf(fidNew, "      %u, %13f, %13f, %13f\n",nodenum, NodesY(i,j), NodesY(i,j+1), NodesY(i,j+2));
            j=j+3;
            end
        end

        fprintf(fidNew, "*Element, type=T3D2\n");
        eltnum=1;
        nodenum=1;
        for i = 1:size(NodesY,1)
           j=4;
           while ~isnan(NodesY(i,j))
               fprintf(fidNew,"%u, %u, %u\n", eltnum, nodenum, nodenum+1);
               eltnum=eltnum+1;
               nodenum=nodenum+1;
           j=j+3;
           end
           nodenum=nodenum+1;
        end
        numYelts=eltnum;
        
        fprintf(fidNew,"*Nset, nset=Fib90, internal, generate\n");
        fprintf(fidNew,"1, %u, 1\n", numYnodes);
        fprintf(fidNew,"*Elset, elset=Fib90, internal, generate\n");
        fprintf(fidNew,"1, %u, 1\n", numYelts);
        fprintf(fidNew, "**Section: TrussSection\n");
        fprintf(fidNew,"*Solid Section, elset=Fib90, material=%s\n%E,\n", Fmat,TArea);
        fprintf(fidNew, "*End Part\n");
        
    %Continue Replicating the Original inp until the end of the first instance  
    tline = fgets(fidOld);
    while ~contains(tline,"*End Instance")
        
        fwrite(fidNew,tline);
        tline = fgets(fidOld);
    end
        fwrite(fidNew,tline);
        
    %Insert Instances of the Truss Parts
    fprintf(fidNew, "**\n*Instance, name=Fibers0, part=Fibers0\n");
    fprintf(fidNew,"%E, %E, %E\n", 0, 0 ,0);
    fprintf(fidNew,"*End Instance\n");
    
    fprintf(fidNew, "**\n*Instance, name=Fibers90, part=Fibers90\n");
    fprintf(fidNew,"%E, %E, %E\n", 0, 0 ,0);
    fprintf(fidNew,"*End Instance\n");
    
    
   %Create Node and Element set of all fiber elements for the embedded constraint
   fprintf(fidNew, "*Nset, nset=ALLFIBERS, instance=FIBERS0, generate\n");
   fprintf(fidNew,"1, %u, 1\n", numXnodes);
   fprintf(fidNew, "*Nset, nset=ALLFIBERS, instance=FIBERS90, generate\n");
   fprintf(fidNew,"1, %u, 1\n", numYnodes);
   fprintf(fidNew, "*Elset, elset=ALLFIBERS, instance=FIBERS0, generate\n");
   fprintf(fidNew,"1, %u, 1\n", numXelts);
   fprintf(fidNew, "*Elset, elset=ALLFIBERS, instance=FIBERS90, generate\n");
   fprintf(fidNew,"1, %u, 1\n", numYelts);
   
   %Finish Replicating the rest of the assembly from the orriginal inp
    tline = fgets(fidOld);
    while ~contains(tline,"*End Assembly") 
        fwrite(fidNew,tline);
        tline = fgets(fidOld);
    end
    
   %Create ALLSURF for All with self contact
    if ALLSELFCONTACT
       fprintf(fidNew, "*Surface, type ELEMENT, name=ALLSURF\n");
       fprintf(fidNew, " ALLCONT\n");
       fprintf(fidNew, " ALLFIBERS\n");
       fprintf(fidNew, "**\n");
    end
    
  %Add the Embedded Element Constraint
    fprintf(fidNew,"** Constraint: Embed-1\n");
    fprintf(fidNew,"*Embedded Element, host elset=HOST\nALLFIBERS\n");
 

   %Finish Replicating the rest of the assembly from the orriginal inp
   if ALLSELFCONTACT
    while ~contains(tline,"INTERACTIONS")     
        fwrite(fidNew,tline);
        tline = fgets(fidOld);
    end
  
    %Create General Contact Interaction
       fwrite(fidNew,tline);
       fprintf(fidNew, "**\n** Interaction: general_contact\n");
       fprintf(fidNew, "*Contact, op=NEW\n");
       fprintf(fidNew, "*Contact Inclusions\n");
       fprintf(fidNew, "ALLSURF , \n");
       fprintf(fidNew, "*Contact Property Assignment\n");
       fprintf(fidNew, " ,  , INTPROP-1\n");
       tline = fgets(fidOld);
    end
    
  %Finish Replicating the rest of the orriginal inp
    while tline~=-1
        fwrite(fidNew,tline);
        tline = fgets(fidOld);
    end
       
    fclose(fidOld);
    fclose(fidNew);
        end
        
   end
end
