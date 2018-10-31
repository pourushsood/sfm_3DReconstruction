function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X)
%% Nonlinear_Triangulation
% AUTHOR: POURUSH SOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

%% Perfoming the triangulation for a fixed number of iterations (10 here) for every point in the 3D space
    for j = 1:10
        for i = 1:length(X)
            X(i,:) = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X(i,:));
        end
    end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
% This function computes the change that needs to be made to reach the
% minima through non linear least squares for a single point in the 3D
% space.
    
  % Computing the location of the 3D point in the image
    x1c = K*R1*(X0'-C1);
    x2c = K*R2*(X0'-C2);
    x3c = K*R3*(X0'-C3);
  % Performing the non linear least squares
    f = [x1c(1)/x1c(3) x1c(2)/x1c(3) x2c(1)/x2c(3) x2c(2)/x2c(3) x3c(1)/x3c(3) x3c(2)/x3c(3)]';
    % Forming the Jacobian Matrix
    J = [(Jacobian_Triangulation(C1,R1,K,x1c));(Jacobian_Triangulation(C2,R2,K,x2c));(Jacobian_Triangulation(C3,R3,K,x3c))];
    b = [x1 x2 x3]';
    dx = ((J'*J)\J')*(b-f);
    X =X0 + dx';

function J = Jacobian_Triangulation(C, R, K, X)
% This function computes the Jacobian for a single point in one of the
% images
  % Derivative terms
    dudx = [K(1,1)*R(1,1) + K(1,3)*R(3,1)    K(1,1)*R(1,2) + K(1,3)*R(3,2)     K(1,1)*R(1,3) + K(1,3)*R(3,3)];
    dvdx = [K(2,2)*R(2,1) + K(2,3)*R(3,1)    K(2,2)*R(2,2) + K(2,3)*R(3,2)     K(2,2)*R(2,3) + K(2,3)*R(3,3)];
    dwdx = [R(3,1) R(3,2) R(3,3)];
  % Jacobian
    w = X(3);
    u = X(1);
    v = X(2);
    dfdx = [(w*dudx - u*dwdx)/(w*w)  ;  (w*dvdx - v*dwdx)/(w*w)];
    J = dfdx;

