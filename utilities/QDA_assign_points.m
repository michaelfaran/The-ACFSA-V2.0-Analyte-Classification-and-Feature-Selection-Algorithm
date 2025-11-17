function labels = QDA_assign_points(X2, classStats, varargin)
% QDA_assign_points  Classify 2D points using QDA built from classStats
%
%   labels = QDA_assign_points(X2, classStats, 'Res',300,'Trunc',5)
%
% Inputs:
%   X2         : N x 2 array of query points
%   classStats : Kx1 struct with fields:
%                  .meanXY  (1x2)
%                  .covXY   (2x2)  -- symmetric positive (semi)definite
% Name-Value (optional; ignored here, just for signature compatibility):
%   'Res'      : grid resolution (unused)
%   'Trunc'    : truncation in sigmas (unused)
%
% Output:
%   labels     : N x 1 vector of class indices in {1,...,K}
%
% Equal priors are assumed.

% ---- parse (ignored) options just to accept your signature
p = inputParser;
p.addParameter('Res', 300, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Trunc', 5, @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:}); %#ok<NASGU>

% ---- extract class parameters
K = numel(classStats); d = 2;
mu = zeros(K,d);  Sigma = zeros(d,d,K);
for k = 1:K
    mu(k,:) = classStats(k).meanXY(:).';
    S = (classStats(k).covXY + classStats(k).covXY.')/2;  % symmetrize
    % tiny jitter to ensure PD for chol:
    t = trace(S)/max(d,1);
    Sigma(:,:,k) = S + (1e-12 + 1e-9*max(t,1e-12))*eye(d);
end
pi_k = ones(K,1)/K;  % equal priors

% ---- precompute Cholesky + log|Sigma|
R = cell(K,1); logDet = zeros(1,K);
for k = 1:K
    try
        R{k} = chol(Sigma(:,:,k),'upper');  % Sigma = R' * R
    catch
        % if still not PD, nudge a bit more
        bump = 1e-6*trace(Sigma(:,:,k))/d;
        R{k} = chol(Sigma(:,:,k)+bump*eye(d),'upper');
        Sigma(:,:,k) = Sigma(:,:,k)+bump*eye(d);
    end
    logDet(k) = 2*sum(log(diag(R{k})));
end

% ---- compute QDA log-posteriors and pick argmax
N = size(X2,1);
logPost = zeros(N,K);
for k = 1:K
    XC = X2 - mu(k,:);
    % Mahalanobis via triangular solve: (X2-mu)/R, then squared norm
    Y = XC / R{k};
    Q = sum(Y.^2,2);
    logPost(:,k) = -0.5*(Q + logDet(k)) + log(pi_k(k));
end
[~, labels] = max(logPost, [], 2);
end
