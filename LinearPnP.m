function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% AUTHOR: POURUSH SOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly


    x = [x ones(length(X),1)];
    % Using inv(K)*x instead of x for numerical stability
    xc = (K\x')';
    j = 1;
    a = zeros(2*length(x),12);
    
    % Forming the matrix corresponding to the least squares problem Ax = 0 
    for i = 1:length(x)
        a(j,:) = [-X(i,1) -X(i,2) -X(i,3) 0 0 0 xc(i,1)*X(i,1) xc(i,1)*X(i,2) X(i,3)*xc(i,1) -1 0 xc(i,1)];
        a(j+1,:) = [0 0 0 -X(i,1) -X(i,2) -X(i,3) xc(i,2)*X(i,1) xc(i,2)*X(i,2) X(i,3)*xc(i,2) 0 -1 xc(i,2)];
        j = j+2;
    end

    % Least squares
    [~,~,V] = svd(a);
    out = V(:,end);
    R = (reshape(out(1:9),[3,3]))';
    t = out(10:12);
    
    % Sign of the determinant
    [U,D,V] = svd(R);
    if det(U*V') > 0
        Rc = U*V';
        tc = t/D(1,1);
    elseif det(U*V') < 0
        Rc = -U*V';
        tc = -t/D(1,1);
    end

    % Computing camera center and rotation
    C = -Rc'*tc;
    R = Rc;




