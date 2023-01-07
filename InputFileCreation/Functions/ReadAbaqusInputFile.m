%Reads information from an Abaqus input file and saves it as structures
%11/29/21

function [BC,CONSTRAINT,ELSET,INSTANCE,NSET,PART,STEP] = ReadAbaqusInputFile(OriginalINP)
% OriginalINP = "FlagshypCube_8h_4t_Combined.inp";
% OriginalINP = "FlagshypCube_8h_0t_ShearAllConstrain.inp";

%Define PART structure array
%I know I don't have to but this is a better way to organize
PART(1).name="name";
PART(1).nodes=zeros(3,3);
PART(1).elements=zeros(3,3);
PART(1).ishost=true;
PART(1).eltType=[];
PART(2)=PART(1);

%Open text file and read data  
   
%Check if file exists and opens sucessfully
   if 1%isfile(OriginalINP)
    fid = fopen(OriginalINP,'r');
        if fid < 0
            %disp('Error')
            fprintf("%s Not Found\n", OriginalINP);
            quit;
        else
            
        %Look for part mesh
        tline = fgetl(fid);
            while ~contains(tline, "*Part, name=")
                tline = fgetl(fid);
            end
            PART(1).name = sscanf(tline, "*Part, name=%s");

            
            tline = fgetl(fid);

        %Read in node coordinates
            Nodes = fscanf(fid, '%*f %*c %f %*c %f %*c %f ',[3, inf]);
            PART(1).nodes = Nodes';
            
        %Read in Element Type
            tline = fgetl(fid);
            PART(1).eltType = tline(1,16:end);

        %Read in element definition 
            if PART(1).eltType == "C3D8R" || PART(1).eltType == "C3D8"
                Elements = fscanf(fid, '%*u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u',[8, inf]);
                PART(1).elements = Elements';
            elseif PART(1).eltType == "T3D2"
                Elements = fscanf(fid, '%*u %*c %u %*c %u',[2, inf]);
                PART(1).elements = Elements';   
            else
                fprintf("Element type not programed \n");
            end

         %Go to *End Part
         while ~contains(tline,"*End Part")
            tline = fgetl(fid);
            if contains(tline,"** Section: ")
                PART(1).material_type = sscanf(tline, "** Section: %s");
                tline = fgetl(fid);
                PART(1).section_type = strip(sscanf(tline,"*%s %s elset=%*s, material=%*s"),',');
                PART(1).section_set = strip(sscanf(tline,"*%s %*s elset=%s material=%*s"),',');
                PART(1).material = sscanf(tline,"*%*s %*s elset=%*s material=%s");
                if strcmp(PART(1).material_type, "TrussSection")
                    PART(1).cross_section = sscanf(fgetl(fid),"%f");
                else
                    PART(1).cross_section = 0;
                end
            end
         end
         tline = fgetl(fid);
         tline = fgetl(fid);
%Get other parts
            nParts = 1;
            while ~contains(tline, "**")
                nParts = nParts + 1;
                i = nParts;
            
                %Look for part mesh
%                 tline = fgetl(fid);
                    while ~contains(tline, "*Part, name=")
                        tline = fgetl(fid);
                    end
                    PART(i).name = sscanf(tline, "*Part, name=%s");
                                        
                    tline = fgetl(fid);
        
                %Read in node coordinates
                    Nodes = fscanf(fid, '%*f %*c %f %*c %f %*c %f ',[3, inf]);
                    PART(i).nodes = Nodes';
                    
                %Read in Element Type
                    tline = fgetl(fid);
                    PART(i).eltType = tline(1,16:end);
        
                %Read in element definition 
                if PART(i).eltType == "C3D8R" || PART(i).eltType == "C3D8"
                    Elements = fscanf(fid, '%*u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u %*c %u',[8, inf]);
                    PART(i).elements = Elements';
                elseif PART(i).eltType == "T3D2"
                    Elements = fscanf(fid, '%*u %*c %u %*c %u',[2, inf]);
                    PART(i).elements = Elements';
                else
                    fprintf("Element type not programed \n");
                end

                while ~contains(tline, "*End Part")
                    tline = fgetl(fid);
                    if contains(tline,"**Section: ") || contains(tline,"** Section: ")
                        PART(i).material_type = sscanf(tline, "** Section: %s");
                        tline = fgetl(fid);
                        PART(i).section_type = strip(sscanf(tline,"*%s %s elset=%*s, material=%*s"),',');
                        PART(i).section_set = strip(sscanf(tline,"*%s %*s elset=%s material=%*s"),',');
                        PART(i).material = sscanf(tline,"*%*s %*s elset=%*s material=%s");
                        if strcmp(PART(i).material_type, "TrussSection")
                            PART(i).cross_section = sscanf(fgetl(fid),"%f");
                        else
                            PART(i).cross_section = 0;
                        end
                    end
                end
                tline = fgetl(fid);
                tline = fgetl(fid);
            end


            %Find the start of the Assembly 
            while ~contains(tline, "** ASSEMBLY")
                tline = fgetl(fid);
            end
            for i=1:3
                tline = fgetl(fid);
            end

            %Read Instances
            nInstance=0;
            tline = fgetl(fid);
            while ~contains(tline, "Nset")
                nInstance = nInstance + 1;
                i = nInstance;
                INSTANCE(i).name = strip(sscanf(tline, "*Instance, name=%s"),",");
                INSTANCE(i).part = sscanf(tline, "*Instance, name=%*s part=%s");

                tline = fgetl(fid);
                if contains(tline, "*End Instance")
                    INSTANCE(i).translate = [0 0 0];
                    INSTANCE(i).rotvec = [0 0 0; 0 0 0];
                    INSTANCE(i).rotangle = 0;
                else
                    INSTANCE(i).translate = sscanf(tline, '%f %*c %f %*c %f ',[3, 1])';
                    tline = fgetl(fid);
                    if contains(tline, "*End Instance")
                        INSTANCE(i).rotvec = [0 0 0; 0 0 0];
                        INSTANCE(i).rotangle = 0;
                    else
                        INSTANCE(i).rotvec = sscanf(tline, '%f %*c %f %*c %f %*c %f %*c %f %*c %f %*c %*f',[3, 2])';
                        INSTANCE(i).rotangle = sscanf(tline,'%*f %*c %*f %*c %*f %*c %*f %*c %*f %*c %*f %*c %f');
                        tline = fgetl(fid);
                    end
                end
            tline = fgetl(fid);
                if contains(tline,"**")
                    tline = fgetl(fid);
            end
            end


            %Read node and element sets in the assembly
            n_Nset=0; n_Elset=0;
            while ~contains(tline, "** Constraint") && ~contains(tline, "*End Assembly")
                if contains(tline,"Nset")
                    n_Nset = n_Nset + 1;
                    NSET(n_Nset).name = strip(sscanf(tline, "*Nset, nset=%s instance=%*s %*s"),',');
                    NSET(n_Nset).instance = strip(sscanf(tline, "*Nset, nset=%*s instance=%s %*s"),',');
                    NSET(n_Nset).generate = sscanf(tline, "*Nset, nset=%*s instance=%*s %s");

                    if NSET(n_Nset).generate =="generate"
                        numbers = fscanf(fid,"%u %*c %u %*c %u",[3,1]);
                        NSET(n_Nset).nodes = [numbers(1):numbers(3):numbers(2)]';
                        tline = fgetl(fid);
                    else
%                         tline = fgetl(fid);
                        NSET(n_Nset).nodes = fscanf(fid, "%u%*c ",[1, inf]);
                    end

                elseif contains(tline,"Elset")
                    n_Elset = n_Elset + 1;
                    ELSET(n_Elset).name = strip(sscanf(tline, "*Elset, elset=%s instance=%*s %*s"),',');
                    ELSET(n_Elset).instance = strip(sscanf(tline, "*Elset, elset=%*s instance=%s %*s"),',');
                    ELSET(n_Elset).generate = sscanf(tline, "*Elset, elset=%*s instance=%*s %s");

                    if ELSET(n_Elset).generate =="generate"
                        numbers = fscanf(fid,"%u %*c %u %*c %u",[3,1]);
                        ELSET(n_Elset).elements = [numbers(1):numbers(3):numbers(2)]';
                        tline = fgetl(fid);
                    else
%                         tline = fgetl(fid);
                        ELSET(n_Elset).elements = fscanf(fid, "%u%*c ",[1, inf]);
                    end
                else
                    fprintf('wtf ReadAbaqusInputFile.m, line 203, not Nset or Elset\n');
                end
                tline = fgetl(fid);

            end

            %Read Constraints to find Embedded Constraint
            tline = fgetl(fid);
            if contains(tline,"*Embed")
                CONSTRAINT.Embed.HostSet = sscanf(tline, "*Embedded Element, host elset=%s");
                tline = fgetl(fid);
                CONSTRAINT.Embed.EmbedSet = sscanf(tline,"%s");
            else
                CONSTRAINT = [];
                tline = fgetl(fid);
            end

            %Read Boundary Conditions
            while ~contains(tline, "** BOUNDARY CONDITIONS")
                tline = fgetl(fid);
            end
            for i=1:2
                tline = fgetl(fid);
            end
            
            nBC = 0;
            while ~contains(tline, "** ------------")
                nBC = nBC + 1;
                if ~contains(tline, "** Name:") && nBC>=2 %if this is a BC with multiple dof fixed
                    BC(nBC).name = BC(nBC-1).name;
                    BC(nBC).type = BC(nBC-1).type;
                    BC(nBC).nset  = strip(sscanf(tline,"%s, "),',');
                    BC(nBC).dof  = sscanf(tline,"%*s %u %*c %u");
                    tline = fgetl(fid);
                    BC(nBC).disp = 0;
                elseif ~contains(tline, "** Name:") && nBC<2
                    fprintf("no fixed BC found\n")
                else
                    BC(nBC).name = sscanf(tline,"** Name: %s Type: %*s");
                    BC(nBC).type = sscanf(tline,"** Name: %*s Type: %s");
                    tline = fgetl(fid);
                    tline = fgetl(fid);
                    BC(nBC).nset  = strip(sscanf(tline,"%s, "),',');
                    BC(nBC).dof  = sscanf(tline,"%*s %u %*c %u");
                    tline = fgetl(fid);
                    BC(nBC).disp = 0;
                end
            end

            %Read Step Time
            while ~contains(tline, "** STEP")
                tline = fgetl(fid);
            end            
            for i=1:4
                tline = fgetl(fid);
            end            
            STEP.time = sscanf(tline, "%*c %f");
            tline = fgetl(fid);
            tline = fgetl(fid);
            STEP.lineardamp = sscanf(tline, "%f %*c %*f");
            STEP.quaddamp = sscanf(tline, "%*f %*c %f");

            %Read Step Boundary Conditions
            while ~contains(tline, "** BOUNDARY CONDITIONS")
                tline = fgetl(fid);
            end
            for i=1:2
                tline = fgetl(fid);
            end
            while ~strcmp(tline, "** ")
                nBC = nBC + 1;
                if ~contains(tline, "** Name:") %if this is a BC with multiple dof specified
                    BC(nBC).name = BC(nBC-1).name;
                    BC(nBC).type = BC(nBC-1).type;
                    tline = fgetl(fid);
                    tline = fgetl(fid);
                    BC(nBC).nset  = strip(sscanf(tline,"%s, "),',');
                    BC(nBC).dof  = sscanf(tline,"%*s %u %*c %u %*c %*f");
                    BC(nBC).disp = sscanf(tline,"%*s %*u %*c %*u %*c %f");
                    tline = fgetl(fid);
                else
                    BC(nBC).name = sscanf(tline,"** Name: %s Type: %*s");
                    BC(nBC).type = sscanf(tline,"** Name: %*s Type: %s");
                    tline = fgetl(fid);
                    tline = fgetl(fid);
                    BC(nBC).nset  = strip(sscanf(tline,"%s, "),',');
                    BC(nBC).dof  = sscanf(tline,"%*s %u %*c %u %*c %*f");
                    BC(nBC).disp = sscanf(tline,"%*s %*u %*c %*u %*c %f");
                    tline = fgetl(fid);
                end
            end



            %Close the inp file
            fclose(fid);

        end
   else
       fprintf("Not a file for some reason\n")
   end 
   
end