function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% AUTHOR: POURUSH SOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

% Running the RANSAC iteration 500 times
for i = 1:500
    % Sampling 8 points randomly
    sampind = randsample(length(x1),8);
    a = zeros([8,9]);
    delta = 0.0001; % Error bound for RANSAC
    inliers = 0;
    maxinliers = 0;
    maxF = 0; % Best Fundamental Matrix
    % Forming the A matrix
    for i1 = 1:8
        u1 = x1(sampind(i1),1);
        v1 = x1(sampind(i1),2);
        u2 = x2(sampind(i1),1);
        v2 = x2(sampind(i1),2);
        a(i1,:) = [u1*u2 v1*u2 u2 u1*v2 v1*v2 v2 u1 v1 1];
    end
    % Computing x
    [~,~, V] = svd(a);
    x = V(:,9);
    % Computing F
    F = (reshape(x,[3,3]))';
    % RANSAC Inliers
    for i2 = 1:length(x2)
        if [x2(i2,:) 1]*F*([x1(i2,:) 1])' < delta
            inliers = inliers + 1;
        end
    end
    if inliers > maxinliers
        maxinliers = inliers;
        maxF = F;
    end
end
F = maxF;
[U,D,V] = svd(F);
D(3,3) = 0;
F = U*D*V';
F = F/norm(F); % Normalizing F
