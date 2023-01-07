function inside = point_in_hexahedron(p,x)
%Function to find if a point is inside of an 8 node hexahedron 
%Based on https://stackoverflow.com/questions/8877872/determining-if-a-point-is-inside-a-polyhedron

%Point of interest p=[p1 p2 p3];
%Coordinates of hex nodes x=[x1 x2 x3 x4 x5 x6 x7 x8; y1 y2 ...] (3x8)
%Flagshyp node connectivity: 1234 5678
%Face definitions, outward normal, CCW
%   4321
%   1265
%   2376
%   3487
%   1584
%   5678
% Define normal by (2-1)x(3-2)

%So what we need to do is check if point p is on the positive or negitive
%sides of each face. If it's on the negitive (oposite the plane normal) of
%all the faces, then it's inside of the hex. 

%Steps
% Loop though faces
%   Calculate face normal: face nodes (n2-n1)x(n3-n2)
%   Choose face node 1 as reference point and define the equation of the
%       plane as normal*(r-n1)=0 where r is any point on the plane
%   Calculate normal*(p-n1). If positive, p is on the + side of the face
%   and is not in the hex. End loop. If negitive (or zero), p is on the - side of the
%   face and could be in the hex. 

%Define a 3x3x12 matrix defining the 6 faces
faces = zeros(3,3,12);
faces(:,:,1) = [x(:,4) x(:,3) x(:,2)];
faces(:,:,2) = [x(:,2) x(:,1) x(:,4)];
faces(:,:,3) = [x(:,1) x(:,2) x(:,6)];
faces(:,:,4) = [x(:,6) x(:,5) x(:,1)];
faces(:,:,5) = [x(:,2) x(:,3) x(:,7)];
faces(:,:,6) = [x(:,7) x(:,6) x(:,2)];
faces(:,:,7) = [x(:,3) x(:,4) x(:,8)];
faces(:,:,8) = [x(:,8) x(:,7) x(:,3)];
faces(:,:,9) = [x(:,1) x(:,5) x(:,8)];
faces(:,:,10) = [x(:,8) x(:,4) x(:,1)];
faces(:,:,11) = [x(:,5) x(:,6) x(:,7)];
faces(:,:,12) = [x(:,7) x(:,8) x(:,5)];

inside=true; 
normals=zeros(3,12);
for i=1:2:12
    A=(faces(:,2,i)-faces(:,1,i));
    B=(faces(:,3,i)-faces(:,2,i));
    n1=cross(A,B);
    normals(:,i)=n1;

    xcc=mean(faces(1,:,i)); ycc=mean(faces(2,:,i)); zcc=mean(faces(3,:,i));
    fc(1,i)=xcc; fc(2,i)=ycc; fc(3,i)=zcc;
    xcc=mean(faces(1,:,i+1)); ycc=mean(faces(2,:,i+1)); zcc=mean(faces(3,:,i+1));
    fc(1,i+1)=xcc; fc(2,i+1)=ycc; fc(3,i+1)=zcc;
    
    A=(faces(:,2,i+1)-faces(:,1,i+1));
    B=(faces(:,3,i+1)-faces(:,2,i+1));
    n2=cross(A,B);
    normals(:,i+1)=n2;

    %Find n*(p-n1)
    d1=(p'-faces(:,1,i));
    d2=(p'-faces(:,1,i+1));
    above1 = n1'*d1;
    above2 = n2'*d2;
%     above1 = above1/(norm(d1)*norm(n1)); %this is supposed to help with numeric stability
%     above2 = above2/(norm(d2)*norm(n2));
    if (above1>1E-10 && above2>1E-10)
        inside=false; 
%         break;
    end
end

%     X=[faces(1,:,1);faces(1,:,2);faces(1,:,3);faces(1,:,4);faces(1,:,5);faces(1,:,6)]'; 
%     Y=[faces(2,:,1);faces(2,:,2);faces(2,:,3);faces(2,:,4);faces(2,:,5);faces(2,:,6)]';
%     Z=[faces(3,:,1);faces(3,:,2);faces(3,:,3);faces(3,:,4);faces(3,:,5);faces(3,:,6)]';
%     
%     XX=[faces(1,:,1);faces(1,:,2);faces(1,:,3);faces(1,:,4);faces(1,:,5);faces(1,:,6)]'; 
%     YY=[faces(2,:,1);faces(2,:,2);faces(2,:,3);faces(2,:,4);faces(2,:,5);faces(2,:,6)]';
%     ZZ=[faces(3,:,1);faces(3,:,2);faces(3,:,3);faces(3,:,4);faces(3,:,5);faces(3,:,6)]';
%     figure(); hold on; view(3); grid on;
%     patch(X,Y,Z,'blue','FaceAlpha', 0.2);
%     patch(XX,YY,ZZ,'red','FaceAlpha', 0.2);
%     plot3(p(1),p(2),p(3),'r.');
% 
%     xc=mean(x(1,:)); yc=mean(x(2,:)); zc=mean(x(3,:));
%     plot3(xc,yc,zc, 'ko');
%     for i=1:12
%         plot3([(normals(1,i)+fc(1,i)) fc(1,i)],[(normals(2,i)+fc(2,i)) fc(2,i)],[(normals(3,i)+fc(3,i)) fc(3,i)], 'g');
%     end
%     for i=1:12
%     plot3([faces(1,:,i) faces(1,1,i)],[faces(2,:,i) faces(2,1,i)],[faces(3,:,i) faces(3,1,i)],'k-','LineWidth',1);
%     end
%     ylim([-0.1 2.01]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     hold off;
end
    

    %% Testing Code
%Function to find if a point is inside of an 8 node hexahedron 
% clear; clc; close all;
% Point of interest p=[p1 p2 p3];
% Given face normals of triangular faces and point on face, determine if the point is inside

% p = [1 0.1 0.5]; 
% x = [1,1,0,0,1,1,0,0;0,1,1,0,0,1,1,0;0,0,0,0,1,1,1,1];
% 
% %For testing
%         F=[1 0 0; 0 1 0; 0 0 1];
%         F=[2 0 0.6; 0 4 0; 0.6 0 1.5];
%         F=[2 -.1 -0.6; -.0 .5 0; -0.6 0 1.7];
% %         x=x-2;
% 
%         xn = zeros(3,8);
%         for i=1:8
%            xn(:,i) = F*x(:,i); 
%         end
%         x=xn;
%So what we need to do is check if point p is on the positive or negitive
%sides of each face. If it's on the negitive (oposite the plane normal) of
%all the faces, then it's inside of the hex. 


% %Define a 3x3x12 matrix defining the 6 faces
% faces = zeros(3,3,12);
% faces(:,:,1) = [x(:,4) x(:,3) x(:,2)];
% faces(:,:,2) = [x(:,2) x(:,1) x(:,4)];
% faces(:,:,3) = [x(:,1) x(:,2) x(:,6)];
% faces(:,:,4) = [x(:,6) x(:,5) x(:,1)];
% faces(:,:,5) = [x(:,2) x(:,3) x(:,7)];
% faces(:,:,6) = [x(:,7) x(:,6) x(:,2)];
% faces(:,:,7) = [x(:,3) x(:,4) x(:,8)];
% faces(:,:,8) = [x(:,8) x(:,7) x(:,3)];
% faces(:,:,9) = [x(:,1) x(:,5) x(:,8)];
% faces(:,:,10) = [x(:,8) x(:,4) x(:,1)];
% faces(:,:,11) = [x(:,5) x(:,6) x(:,7)];
% faces(:,:,12) = [x(:,7) x(:,8) x(:,5)];
% 
% inside=true; 
% normals=zeros(3,12);
% xc=mean(x(1,:)); yc=mean(x(2,:)); zc=mean(x(3,:));
% % p=[xc yc zc];
% fc=zeros(3,12);
% 
% for i=1:2:12
%     A=(faces(:,2,i)-faces(:,1,i));
%     B=(faces(:,3,i)-faces(:,2,i));
%     n1=cross(A,B);
%     normals(:,i)=n1;
%     
%     A=(faces(:,2,i+1)-faces(:,1,i+1));
%     B=(faces(:,3,i+1)-faces(:,2,i+1));
%     n2=cross(A,B);
%     normals(:,i+1)=n2;
%     
%     xcc=mean(faces(1,:,i)); ycc=mean(faces(2,:,i)); zcc=mean(faces(3,:,i));
%     fc(1,i)=xcc; fc(2,i)=ycc; fc(3,i)=zcc;
%     xcc=mean(faces(1,:,i+1)); ycc=mean(faces(2,:,i+1)); zcc=mean(faces(3,:,i+1));
%     fc(1,i+1)=xcc; fc(2,i+1)=ycc; fc(3,i+1)=zcc;
% 
%     %Find n*(p-n1)
%     d1=(p'-faces(:,1,i));
%     d2=(p'-faces(:,1,i+1));
%     above1 = n1'*d1;
%     above2 = n2'*d2;
%     above1 = above1/(norm(d1)*norm(n1)) %this is supposed to help with numeric stability
%     above2 = above2/(norm(d2)*norm(n2))
%     if (above1>1E-10 && above2>1E-10)
%         inside=false; 
% %         disp(above)
% %         break;
%     end
% end
% 
% inside
% 
% %Also for testing
% %Also for testing
% 
%     X=[faces(1,:,1);faces(1,:,2);faces(1,:,3);faces(1,:,4);faces(1,:,5);faces(1,:,6)]'; 
%     Y=[faces(2,:,1);faces(2,:,2);faces(2,:,3);faces(2,:,4);faces(2,:,5);faces(2,:,6)]';
%     Z=[faces(3,:,1);faces(3,:,2);faces(3,:,3);faces(3,:,4);faces(3,:,5);faces(3,:,6)]';
%     
%     XX=[faces(1,:,1);faces(1,:,2);faces(1,:,3);faces(1,:,4);faces(1,:,5);faces(1,:,6)]'; 
%     YY=[faces(2,:,1);faces(2,:,2);faces(2,:,3);faces(2,:,4);faces(2,:,5);faces(2,:,6)]';
%     ZZ=[faces(3,:,1);faces(3,:,2);faces(3,:,3);faces(3,:,4);faces(3,:,5);faces(3,:,6)]';
%   
% 
%     figure(); hold on; view(3); grid on;
%     patch(X,Y,Z,'blue','FaceAlpha', 0.2);
%     patch(XX,YY,ZZ,'red','FaceAlpha', 0.2);
%     plot3(p(1),p(2),p(3),'r.');
% 
%     plot3(xc,yc,zc, 'ko');
%     for i=1:12
%         plot3([(normals(1,i)+fc(1,i)) fc(1,i)],[(normals(2,i)+fc(2,i)) fc(2,i)],[(normals(3,i)+fc(3,i)) fc(3,i)], 'g');
%     end
%     for i=1:12
%     plot3([faces(1,:,i) faces(1,1,i)],[faces(2,:,i) faces(2,1,i)],[faces(3,:,i) faces(3,1,i)],'k-','LineWidth',1);
%     end
% %     ylim([-0.1 -0.06]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     hold off;
% 
% end
%     
   