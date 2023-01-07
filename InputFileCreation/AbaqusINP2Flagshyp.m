%Convert Abaqus Input file to Flaghsyp format
%11/29/21

%assuming there are no more than 2 element types and they are C3D8 or T3D2 
%and assuming an explicit model

%Someday update to read the correct material properties from Abaqus
%Also need to update to read Abaqus files without embedded elements

%Part structure holds part name, nodes, mesh
%% User Input

clear; clc; close all; 

OriginalINP = "DyneemaCube_10000fibers.inp";
NewINP = "DyneemaCube_10000fibers.dat";   

% OriginalINP = "RussellTensile-1-5_5000Fibers7_discritized.inp";
% NewINP = "RussellTensile-1-5_5000Fibers7_discritized.dat";
% 
% OriginalINP = "AttwoodCompression-1_1000Fibers7_discritized.inp";
% NewINP = "AttwoodCompression-1_1000Fibers7_discritized.dat";
% 
% OriginalINP = "RussellTensile-Half_5000Fibers7_discritized.inp";
% NewINP = "RussellTensile-Half_5000Fibers7_discritized.dat";

[BC,CONSTRAINT,ELSET,INSTANCE,NSET,PART,STEP] = ReadAbaqusInputFile(OriginalINP);

%% Calculations to convert to Flagshyp format
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
            n_nodes = n_nodes_new;

            INSTANCE(i).eltType = PART(j).eltType;
            if contains(PART(j).eltType, "C3D8")
                n_elements(1) = n_elements(1) + size(PART(j).elements,1);
            elseif contains(PART(j).eltType,"T3D2")
                n_elements(2) = n_elements(2) + size(PART(j).elements,1);
            end

            break;
        end
    end
    if n_elements(1) == 0 || n_elements(2) == 0
        n_elt_type = 1;
    else
        n_elt_type = 2;
    end



        %Translate and Rotate the instance to global position
        if INSTANCE(i).rotangle == 0
            R = eye(3);
        else
            u = INSTANCE(i).rotvec(2,:) - INSTANCE(i).rotvec(1,:);
            unit = u/abs(u);
            c = cosd(INSTANCE(i).rotangle); s = sind(INSTANCE(i).rotangle);
            R = [c+u(1)^2*(1-c), u(1)*u(2)*(1-c)-u(3)*s, u(1)*u(3)*(1-c)+u(2)*s; u(1)*u(2)*(1-c)+u(3)*s, c+u(2)^2*(1-c), u(2)*u(3)*(1-c)-u(1)*s; u(1)*u(3)*(1-c)-u(2)*s, u(2)*u(3)*(1-c)+u(1)*s, c+u(3)^2*(1-c)];
        end

        for j=1:size(INSTANCE(i).global_nodes,1)
            INSTANCE(i).global_nodes(j,2:4) = INSTANCE(i).global_nodes(j,2:4) + INSTANCE(i).translate;
            x = INSTANCE(i).global_nodes(j,2:4) - INSTANCE(i).rotvec(1,:);
            x = R*x';
            INSTANCE(i).global_nodes(j,2:4) = x' + INSTANCE(i).rotvec(1,:);
        end
end



%Determine Boundary Conditions of Each Global Node
BCmatrix = zeros(n_nodes,4);
for i=1:length(BC)
    switch BC(i).dof(2)
        case 1
            multiply=100;
        case 2
            multiply=10;
        case 3
            multiply=1;
    end

    for j=1:length(NSET)
        if strcmp(BC(i).nset, NSET(j).name)
            for k = 1:length(INSTANCE)
                if strcmp(NSET(j).instance, INSTANCE(k).name)
                    global_nodes = INSTANCE(k).global_nodes(NSET(j).nodes,1);
                    BCmatrix(global_nodes,1) = BCmatrix(global_nodes,1) + ones(length(global_nodes),1)*multiply;
                    if BC(i).disp ~= 0
                        BCmatrix(global_nodes,BC(i).dof(2)+1) = ones(length(global_nodes),1) * BC(i).disp;
%                         BCmatrix(global_nodes,3) = ones(length(global_nodes),1) * BC(i).dof(2);
                    end
                end
            end
            break;
        end
    end
end

perscribed_dof = zeros(n_nodes,3);
perscribed_x = zeros(n_nodes,3);
perscribed_y = zeros(n_nodes,3);
perscribed_z = zeros(n_nodes,3);
for i=1:n_nodes
    switch BCmatrix(i,1)
        case 0
            BCmatrix(i,1) = 0; 
        case 100
            BCmatrix(i,1) = 1; 
        case 10
            BCmatrix(i,1) = 2; 
        case 110
            BCmatrix(i,1) = 3; 
        case 1
            BCmatrix(i,1) = 4; 
        case 101
            BCmatrix(i,1) = 5; 
        case 11
            BCmatrix(i,1) = 6; 
        case 111
            BCmatrix(i,1) = 7; 
    end

    
    if BCmatrix(i,2)~=0
        perscribed_x(i,1) = i;
        perscribed_x(i,2) = 1;
        perscribed_x(i,3) = BCmatrix(i,2);
    end
    if BCmatrix(i,3)~=0
        perscribed_y(i,1) = i;
        perscribed_y(i,2) = 2;
        perscribed_y(i,3) = BCmatrix(i,3);
    end
    if BCmatrix(i,4)~=0
        perscribed_z(i,1) = i;
        perscribed_z(i,2) = 3;
        perscribed_z(i,3) = BCmatrix(i,4);
    end
    
end
perscribed_x = perscribed_x(perscribed_x(:,3)~=0,:);
perscribed_y = perscribed_y(perscribed_y(:,3)~=0,:);
perscribed_z = perscribed_z(perscribed_z(:,3)~=0,:);

perscribed_dof = [perscribed_x; perscribed_y; perscribed_z];


%Flag elements as embedded(1) or host(0) elements 
embed_flag = 0;
if ~isempty(CONSTRAINT)
    embed_flag = 1;
    for j=1:length(ELSET)
        if strcmpi(CONSTRAINT.Embed.HostSet, ELSET(j).name)
            for k = 1:length(INSTANCE)
                if strcmpi(ELSET(j).instance, INSTANCE(k).name)
                    INSTANCE(k).elt_ishost = zeros(size(INSTANCE(k).global_elements,1),1);
                end
            end
        end
        if strcmpi(CONSTRAINT.Embed.EmbedSet, ELSET(j).name)
            for k = 1:length(INSTANCE)
                if strcmpi(ELSET(j).instance, INSTANCE(k).name)
                    INSTANCE(k).elt_ishost = ones(size(INSTANCE(k).global_elements,1),1);
                end
            end
        end
    end

    %Add boundary conditions to embedded nodes
    for j=1:length(NSET)
        if strcmpi(CONSTRAINT.Embed.EmbedSet, NSET(j).name)
            for k = 1:length(INSTANCE)
                if strcmpi(NSET(j).instance, INSTANCE(k).name)
                    BCmatrix(INSTANCE(k).global_nodes(:,1),1) = 8;
                end
            end
        end
    end
else
    for k = 1:length(INSTANCE)
        INSTANCE(k).elt_ishost = zeros(size(INSTANCE(k).global_elements,1),1);
    end
end




%% Write Flagshyp Input
fid = fopen(NewINP,'w');

fprintf(fid, "3-D EmbeddedElt \nExplicit Analysis	1\nEmbedded Elements	%u\nVolumeCorrection 	0\n",embed_flag);
if n_elt_type == 2
    fprintf(fid,"%u\nhexa8\ntruss2\n%u\n",n_elt_type,n_nodes);
elseif n_elt_type == 1
    if contains(PART(1).eltType, "C3D8")
        fprintf(fid,"%u\nhexa8\n%u\n",n_elt_type,n_nodes);
    elseif contains(PART(1).eltType, "T3D2")
        fprintf(fid,"%u\ntruss2\n%u\n",n_elt_type,n_nodes);
    else
        fprintf("what is this element?\n");
    end
end

for i=1:length(INSTANCE)
    for j=1:size(INSTANCE(i).global_nodes,1)
        node = INSTANCE(i).global_nodes(j,1);
        fprintf(fid, "%-4u %-4u %-4.5f %-4.5f %-4.5f\n",node,BCmatrix(node),INSTANCE(i).global_nodes(j,2),INSTANCE(i).global_nodes(j,3),INSTANCE(i).global_nodes(j,4));
    end
end


info1 = zeros(n_elements(1,1), 10); count1 = 0; 
info2 = zeros(n_elements(2,1), 4); count2 = 0;
for i=1:1:length(INSTANCE) 
    if contains(INSTANCE(i).eltType, "C3D8")
        nelt1 = size(INSTANCE(i).global_elements,1);
        info1([count1+1:nelt1+count1],1) = ones(size(INSTANCE(i).global_elements,1),1); %this is material type. we just assume it's 1
        info1([count1+1:nelt1+count1],2) = INSTANCE(i).elt_ishost;
        info1([count1+1:nelt1+count1],3:10) = INSTANCE(i).global_elements;
        count1 = count1 + size(INSTANCE(i).global_elements,1);
    elseif contains(INSTANCE(i).eltType,"T3D2")
        nelt2 = size(INSTANCE(i).global_elements,1);
        info2([count2+1:nelt2+count2],1) = ones(size(INSTANCE(i).global_elements,1),1); %this is material type. we just assume it's 1
        info2([count2+1:nelt2+count2],2) = INSTANCE(i).elt_ishost;
        info2([count2+1:nelt2+count2],3:4) = INSTANCE(i).global_elements;
        count2 = count2 + size(INSTANCE(i).global_elements,1);
    end
end

fprintf(fid,"%u\n",n_elements(1));
format = [repmat("%3u ",1,11) "\n"] ;
format = "%-3u %-3u %-3u %-3u %-3u %-3u %-3u %-3u %-3u %-3u %-3u \n";
for i=1:n_elements(1)
    fprintf(fid,format,[i info1(i,:)]');
end
if n_elements(2)~=0
    fprintf(fid,"%u\n",n_elements(2));
    format = "%-3u %-3u %-3u %-3u %-3u \n";
    for i=1:n_elements(2)
        fprintf(fid,format,[i info2(i,:)]');
    end
end

%Print material definitions but blank (ReadAbaqusInputFile doesn't read
%material definitions yet)

for i=1:n_elt_type
    fprintf(fid,"1\n");
    fprintf(fid,"1 1 7800.0 76.92e9 115.4e9 NeoHooke\n");
    fprintf(fid,"1 2 7800.0 2e11 0.3 %f 2E6 1 NeoHookeTruss\n",PART(i).cross_section);
end

%Print step time
fprintf(fid,"%.4f\n",STEP.time);

fprintf(fid, "0 %u 0 0.0 0.0 0.0\n",size(perscribed_dof,1));
for i=1:size(perscribed_dof,1)
    fprintf(fid, "%-3u %-3u %-6.5f \n",perscribed_dof(i,:));
end

fprintf(fid,"1 20.0 1.0 25 1.e-6 0.0 0.0 1 1 3 2");
fclose("all");