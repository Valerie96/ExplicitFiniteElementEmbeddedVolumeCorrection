function [max_l] =calc_max_element_size(xlocal,element_type)

% global_nodes    = FEM.mesh.connectivity(:,ielement);
% xlocal          = GEOM.x(:,global_nodes) ;


    switch element_type 
        case'hexa8'
            % examine each each edge of hex 
            % there are 12 edges in a hex
            lengths=zeros(1,12);
    
            % Nodes 1-2
            lengths(1) = sqrt((xlocal(1,2)-xlocal(1,1))^2 + (xlocal(2,2)-xlocal(2,1))^2 + (xlocal(3,2)-xlocal(3,1))^2);
            % Nodes 2-3
            lengths(2) = sqrt((xlocal(1,3)-xlocal(1,2))^2 + (xlocal(2,3)-xlocal(2,2))^2 + (xlocal(3,3)-xlocal(3,2))^2);
            % Nodes 3-4
            lengths(3) = sqrt((xlocal(1,4)-xlocal(1,3))^2 + (xlocal(2,4)-xlocal(2,3))^2 + (xlocal(3,4)-xlocal(3,3))^2);
            % Nodes 4-1
            lengths(4) = sqrt((xlocal(1,1)-xlocal(1,4))^2 + (xlocal(2,1)-xlocal(2,4))^2 + (xlocal(3,1)-xlocal(3,4))^2);
    
            % Nodes 5-6
            lengths(5) = sqrt((xlocal(1,6)-xlocal(1,5))^2 + (xlocal(2,6)-xlocal(2,5))^2 + (xlocal(3,6)-xlocal(3,5))^2);
            % Nodes 6-7
            lengths(6) = sqrt((xlocal(1,7)-xlocal(1,6))^2 + (xlocal(2,7)-xlocal(2,6))^2 + (xlocal(3,7)-xlocal(3,6))^2);
            % Nodes 7-8
            lengths(7) = sqrt((xlocal(1,8)-xlocal(1,7))^2 + (xlocal(2,8)-xlocal(2,7))^2 + (xlocal(3,8)-xlocal(3,7))^2);
            % Nodes 8-5
            lengths(8) = sqrt((xlocal(1,5)-xlocal(1,8))^2 + (xlocal(2,5)-xlocal(2,8))^2 + (xlocal(3,5)-xlocal(3,8))^2);
    
            % Nodes 1-5
            lengths(9) = sqrt((xlocal(1,5)-xlocal(1,1))^2 + (xlocal(2,5)-xlocal(2,1))^2 + (xlocal(3,5)-xlocal(3,1))^2);
            % Nodes 2-6
            lengths(10) = sqrt((xlocal(1,6)-xlocal(1,2))^2 + (xlocal(2,6)-xlocal(2,2))^2 + (xlocal(3,6)-xlocal(3,2))^2);
            % Nodes 3-7
            lengths(11) = sqrt((xlocal(1,7)-xlocal(1,3))^2 + (xlocal(2,7)-xlocal(2,3))^2 + (xlocal(3,7)-xlocal(3,3))^2);
            % Nodes 4-8
            lengths(12) = sqrt((xlocal(1,8)-xlocal(1,4))^2 + (xlocal(2,8)-xlocal(2,4))^2 + (xlocal(3,8)-xlocal(3,4))^2);
    
%             le = min(lengths);
            max_l = max(numbers);
            
        case 'truss2'
%             le = norm(xlocal(:,2) - xlocal(:,1));
            max_l = le;
        
        otherwise
       % error not yet programmed ....
       %disp ('From SolutionEquations/CalculateElementSize.m: ')
       %disp ('element type size calculation is not yet programmed... ')
       
    end


end