function p = find_intersection(X,FEM,e_connectivity,host1,host2)
%Given 2 adjacent host elements, this function finds the shared face and
%can calculate the intersection of that shared face and a truss element
%that spans both hosts


%Get connectivity of host1 and host 2
h1_connect = FEM(1).mesh.connectivity(:,host1);
h2_connect = FEM(1).mesh.connectivity(:,host2);

%Construct the 6 faces of each host element
faces1 = zeros(6,4);
faces1(1,:) = [h1_connect(4) h1_connect(3) h1_connect(2) h1_connect(1)];
faces1(2,:) = [h1_connect(1) h1_connect(2) h1_connect(6) h1_connect(5)];
faces1(3,:) = [h1_connect(2) h1_connect(3) h1_connect(7) h1_connect(6)];
faces1(4,:) = [h1_connect(3) h1_connect(4) h1_connect(8) h1_connect(7)];
faces1(5,:) = [h1_connect(1) h1_connect(5) h1_connect(8) h1_connect(4)];
faces1(6,:) = [h1_connect(5) h1_connect(6) h1_connect(7) h1_connect(8)];

faces2 = zeros(6,4);
%Defined backwards so that the normal of the shared face will be in the
%same direction as on host1
faces2(1,:) = [h2_connect(1) h2_connect(2) h2_connect(3) h2_connect(4)];
faces2(2,:) = [h2_connect(5) h2_connect(6) h2_connect(2) h2_connect(1)];
faces2(3,:) = [h2_connect(6) h2_connect(7) h2_connect(3) h2_connect(2)];
faces2(4,:) = [h2_connect(7) h2_connect(8) h2_connect(4) h2_connect(3)];
faces2(5,:) = [h2_connect(4) h2_connect(8) h2_connect(5) h2_connect(1)];
faces2(6,:) = [h2_connect(8) h2_connect(7) h2_connect(6) h2_connect(5)];

%Compare Host1 faces with Host2 faces to find a shared face
%The node numbers may not be in the same order for each host, but always 
%the same sequence (ie 1 2 3 4 vs 3 4 1 2)

SharedFace = zeros(1,4);
for i1 = 1:6
    for i2 = 1:6
        for j2 = 1:4
%             next_num = mod(j2,4)+1;
%             fprintf("next_num %u\n",next_num);
            if faces1(i1,1) == faces2(i2,j2)
                %If matching nodes are found, check the next 3 nodes. If
                %the index reaches 4 it has to loop back to 1 until all
                %four face nodes have been compared. 
                next_num = mod(j2,4)+1;
                if faces1(i1,2) == faces2(i2,next_num) && faces1(i1,3) == faces2(i2,mod(j2,4)+2)
                    SharedFace = faces1(i1,:);
                end
            end
        end
        
        
    end
end
    
    
    
if SharedFace(1,1) == 0
    fprintf("Shared face of host elements not found. Check mesh\n");
end

%Find intersection of SharedFace and the line between the two truss nodes
plane_line1 = X(:,SharedFace(1,2)) - X(:,SharedFace(1,1));
plane_line2 = X(:,SharedFace(1,3)) - X(:,SharedFace(1,2));
plane_normal = cross(plane_line1,plane_line2);
plane_normal = plane_normal./norm(plane_normal);

truss_line = X(:,e_connectivity(2)) - X(:,e_connectivity(1));
truss_line = truss_line./norm(truss_line);

%Plane eq: (p-p0).n = 0
%Line eq: p = (Lo + ld)
%Solve for d = ((p0-lo).n)/(l.n)
%Intersection point is p = (Lo + ld)
p0 = X(:,SharedFace(1,4)); Lo = X(:,e_connectivity(1));
d = ((p0-Lo)'*plane_normal);
d = d/(truss_line'*plane_normal);

p = Lo + truss_line*d;


end