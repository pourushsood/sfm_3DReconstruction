function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% AUTHOR: POURUSH SOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 3) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 3) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points

    X = zeros(length(x1),3);
    % Iterating over all points
    for i1 = 1:length(x1)
        a = x1(i1,1);
        b = x1(i1,2);
        c = 1;
        % vec2skew for point in image 1
        x1mat = [0,-c,b;c,0,-a;-b,a,0];
        P1 = K*R1*[eye(3) -C1];
        a = x2(i1,1);
        b = x2(i1,2);
        c = 1;   
        % vec2skew for point in image 2
        x2mat = [0,-c,b;c,0,-a;-b,a,0];
        P2 = K*R2*[eye(3) -C2];
        % Forming the least squares matrix [X]xP=0
        A = vertcat(x1mat*P1,x2mat*P2);
        [~,~,V] = svd(A);
        X(i1,:) = V(1:3,4)/V(4,4);
    end
    


