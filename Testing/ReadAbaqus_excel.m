%Read Abaqus Output files in excel format
function [AbqEHost, AbqETruss, AbqE] = ReadAbaqus_excel(name)

%Open the data file, either a xltx or xlsx
if isfile(strcat(name, '.xltx'))
    C=readcell(strcat(name, '.xltx'));
elseif isfile(strcat(name, '.xlsx'))
    C=readcell(strcat(name, '.xlsx'));
else
    fprintf("bleh where's the file\n");
    fprintf(name); fprintf("\n")
    C=readcell(strcat(name, '.xlsx'));
end

%% Read Energy Data
    %Save Column Titles
    Titles = C(1,:);
    C(1:3,:) = [];
    
    %Count the number of energy output steps
    n_energyout = 0;
    t = (C{1,1});
    while ~ismissing(t)
        n_energyout = n_energyout + 1;
        t = (C{1+n_energyout,1});    
    end
    
    %Get energy data. Easy
    Energy = cell2mat(C(1:n_energyout,1:12));
        AbqE.time = Energy(:,1);
        AbqE.AE = Energy(:,2);
        % AbqE.CD = Energy(:,3);
        % AbqE.MD = Energy(:,4);
        % AbqE.FD = Energy(:,5);
        AbqE.IE = Energy(:,6);
        AbqE.KE = Energy(:,7);
        % AbqE.PD = Energy(:,8);
        AbqE.SE = Energy(:,9);
        AbqE.VD = Energy(:,10);
        AbqE.WK = Energy(:,11);
        AbqE.ETOTAL = Energy(:,12);
    C(1:n_energyout+4,:) = [];

%% Read Stress and Strain Field Output Data
    %Count the number of instances and elements for the field output
    n_inst = 0; n_elt = 0; inst_prev = "-";
    
    %Save Column Titles
    Titles2 = C(1,:); C(1,:) = [];
    
    time = sscanf(C{1,3},"Increment         %*u: Step Time = %f");
    while time == 0
        n_elt = n_elt + 1;
        inst = sscanf(C{n_elt,4}, "%s");
            
        if ~strcmp(inst, inst_prev)
            n_inst = n_inst + 1;
        end
        inst_prev = inst;
        
        %Read cell column 3 to detrmine the step time
        time = sscanf(C{n_elt+1,3},"Increment         %*u: Step Time = %f");
    end
    INSTANCE = struct('name', cell(1, n_inst), 'n_elt', cell(1, n_inst), 'n_node', cell(1, n_inst));
    
    %Count n_elements per instance
    inst_prev = "-"; instn = 0;
        for i=1:n_elt
            inst = sscanf(C{i,4}, "%s");
            element = cell2mat(C(i,5));
            
            %if this is a new instance
            if ~strcmp(inst, inst_prev)
                instn = instn + 1;
                INSTANCE(instn).name = inst;
                INSTANCE(instn).n_elt = 1;
                INSTANCE(instn).n_node = 0;
                
                if strcmp(C{i,14}, "NoValue")
                    INSTANCE(instn).istruss = 1;
                else
                    INSTANCE(instn).istruss = 0;
                end
    
            %Same as previous instance
            else
                INSTANCE(instn).n_elt = INSTANCE(instn).n_elt + 1;
            end
            inst_prev = inst;
        end
    
    %Count the number of output steps based on the number of elements
    n_fieldoutput=0;
    tline=sscanf(C{1,3},"%s");
    
    while ~contains(tline,"Frame")
        n_fieldoutput=n_fieldoutput+1;
        tline=sscanf(C{1+n_fieldoutput*n_elt,3},"%s");
    end
    
    %Get stress and strain
    %Loop through each time increment, then through the number of instances,
    %then number of elements per instance
    % LE(time, ij, element) S(time, ij, element)
    line = 1; time = zeros(n_fieldoutput,1);
    for k=1:n_fieldoutput
        if strcmp(sscanf(C{line,3},"%s"),"Frame")
            break;
        end
        time(k) = sscanf(C{line,3},"Increment         %*u: Step Time = %f");
        for j=1:n_inst
            INSTANCE(j).elements=zeros(INSTANCE(j).n_elt,1);
            for i=1:INSTANCE(j).n_elt
                INSTANCE(j).elements(i) = cell2mat(C(line,5));
                
    %             fprintf("line %u    element %u   k%u j%u i%u\n",line, element, k, j, i);
                    %Save data for a truss
                    if INSTANCE(j).istruss == 1
                        INSTANCE(j).LE(k,1,i) = cell2mat(C(line,12));
                        INSTANCE(j).S(k,1,i) = cell2mat(C(line,18));
                    %Save data for a hexahedron
                    else
                        INSTANCE(j).istruss = 0;
                        INSTANCE(j).LE(k,1:6,i) = cell2mat(C(line,12:17));
                        INSTANCE(j).S(k,1:6,i) = cell2mat(C(line,18:23));
                    end
            line = line + 1;
            end
        end
    end
    
    %Delete data we just read
    C(1:line-1,:) = [];

%% Read Acc,Vel,Disp Field Output Data
    %Save Column Titles
    Titles3 = C(1,:);
    C(1,:) = [];
    
    %Count n_nodes total and n_nodes per instance
    t = sscanf(C{1,3},"Increment         %*u: Step Time = %f");
    inst_prev = "-"; instn = 0; n_node = 0; 
    line = 1;
        while t == 0
            instn = instn + 1;
            inst = sscanf(C{line,4}, "%s");
            node = cell2mat(C(line,5));
            
            %if this is a new instance
            if ~strcmp(inst, inst_prev)
                INSTANCE(instn).n_node = INSTANCE(instn).n_node + 1;
            %Same as previous instance
            else
                instn = instn - 1;
                INSTANCE(instn).n_node = INSTANCE(instn).n_node + 1;
            end
    %         fprintf("n_node %u  node# %u instn%u\n", n_node, node,instn);
            inst_prev = inst;
            n_node = n_node + 1;
            line = line + 1;
            t = sscanf(C{line,3},"Increment         %*u: Step Time = %f");
        end
    
    
    %Get nodal data
    % A RF V U (time, i, node)
    line = 1;
    for k=1:n_fieldoutput
    %     C{line,3}
        t = sscanf(C{line,3},"Increment         %*u: Step Time = %f");
        for j=1:n_inst
            INSTANCE(j).nodes=zeros(INSTANCE(j).n_node,1);
            for i=1:INSTANCE(j).n_node
                INSTANCE(j).nodes(i) = cell2mat(C(line,5));
                
    %             fprintf("line %u    node %u   k%u j%u i%u\n",line, node, k, j, i);
                %Save data
                    INSTANCE(j).A(k,1:3,i) = cell2mat(C(line,12:14));
                    INSTANCE(j).RF(k,1:3,i) = cell2mat(C(line,15:17));
                    INSTANCE(j).U(k,1:3,i) = cell2mat(C(line,18:20));
                    INSTANCE(j).V(k,1:3,i) = cell2mat(C(line,21:23));
    
            line = line + 1;
            end
        end
    end

%% Organize data into AbqEHost, AbqETruss
hostnode = 0; hostelt = 0;
trussnode = 0; trusselt = 0;
AbqETruss.time = time;
AbqEHost.time = time;

for i=1:n_inst
% fprintf("%s\n", INSTANCE(i).name);
    if INSTANCE(i).istruss == 1
       trusselt_new = trusselt + INSTANCE(i).n_elt;
       trussnode_new = trussnode + INSTANCE(i).n_node;
        AbqETruss.LE(:,:,trusselt+1:trusselt_new) = INSTANCE(i).LE;
        AbqETruss.S(:,:,trusselt+1:trusselt_new) = INSTANCE(i).S;

        AbqETruss.A(:,:,trussnode+1:trussnode_new) = INSTANCE(i).A;
        AbqETruss.RF(:,:,trussnode+1:trussnode_new) = INSTANCE(i).RF;
        AbqETruss.U(:,:,trussnode+1:trussnode_new) = INSTANCE(i).U;
        AbqETruss.V(:,:,trussnode+1:trussnode_new) = INSTANCE(i).V;

        AbqETruss.elements(1,trusselt+1:trusselt_new)=INSTANCE(i).elements;
        AbqETruss.nodes(1,trussnode+1:trussnode_new)=INSTANCE(i).nodes;
        trusselt = trusselt_new;
        trussnode = trussnode_new;
        
    else
       hostelt_new = hostelt + INSTANCE(i).n_elt;
       hostnode_new = hostnode + INSTANCE(i).n_node;
        AbqEHost.LE(:,:,hostelt+1:hostelt_new) = [INSTANCE(i).LE(:,1,:) INSTANCE(i).LE(:,4,:) INSTANCE(i).LE(:,5,:) INSTANCE(i).LE(:,2,:) INSTANCE(i).LE(:,6,:) INSTANCE(i).LE(:,3,:)];
        AbqEHost.S(:,:,hostelt+1:hostelt_new) = [INSTANCE(i).S(:,1,:) INSTANCE(i).S(:,4,:) INSTANCE(i).S(:,5,:) INSTANCE(i).S(:,2,:) INSTANCE(i).S(:,6,:) INSTANCE(i).S(:,3,:)];

        AbqEHost.A(:,:,hostnode+1:hostnode_new) = INSTANCE(i).A;
        AbqEHost.RF(:,:,hostnode+1:hostnode_new) = INSTANCE(i).RF;
        AbqEHost.U(:,:,hostnode+1:hostnode_new) = INSTANCE(i).U;
        AbqEHost.V(:,:,hostnode+1:hostnode_new) = INSTANCE(i).V;

        AbqEHost.elements(1,hostelt+1:hostelt_new)=INSTANCE(i).elements;
        AbqEHost.nodes(1,hostnode+1:hostnode_new)=INSTANCE(i).nodes;
    end
    
end
