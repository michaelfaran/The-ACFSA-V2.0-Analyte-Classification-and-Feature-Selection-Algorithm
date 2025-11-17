function Gaussian_error = Gaussian_Witch_Practice_QDA(classStats, varargin)
% Gaussian_Witch_Practice_QDA  Approx. QDA misclassification prob per class
%
%   Gaussian_error = Gaussian_Witch_Practice_QDA(classStats)
%   Gaussian_error = Gaussian_Witch_Practice_QDA(classStats, 'Res', 300, 'Trunc', 5)
%
% INPUT:
%   classStats : Kx1 struct with fields
%       .meanXY : 1x2 mean vector
%       .covXY  : 2x2 covariance matrix (symmetric, PD)
% OPTIONS (name/value):
%   'Res'   : grid resolution per axis for the integrals (default 300)
%   'Trunc' : truncation in std units for the integration box (default 5)
%
% OUTPUT:
%   Gaussian_error : Kx1 vector, misclassification probability per class
%
% METHOD:
%   - Equal priors.
%   - QDA discriminants from the provided (mu_k, Sigma_k).
%   - For each class k, integrate its N(mu_k,Sigma_k) over the region
%     assigned to class k by QDA; error_k = 1 - mass(correctly classified).
%   - Integration is approximated on a grid inside a ±Trunc·σ box around mu_k.
%
% NOTES:
%   - Tail mass beyond ±5σ is ~1e-6 (ignored).
%   - No toolbox dependencies; uses Cholesky for stable Mahalanobis and log|Σ|.
%
% Michael’s friendly witchcraft 🧙‍♂️—but statistically legit ;)

% -------- options
p = inputParser;
p.addParameter('Res',   300, @(x)isnumeric(x)&&isscalar(x)&&x>=50);
p.addParameter('Trunc',   5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Jitter', 1e-9, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.parse(varargin{:});
res   = p.Results.Res;
trunc = p.Results.Trunc;
jitter= p.Results.Jitter;

% -------- pack class parameters
K = numel(classStats);
d = 2;
mu = zeros(K,d);
Sigma = zeros(d,d,K);
pi_k = ones(K,1)/K;  % equal priors

for k = 1:K
    mu(k,:)     = classStats(k).meanXY(:).';
    S           = (classStats(k).covXY + classStats(k).covXY.')/2; % symmetrize
    % Tiny jitter if needed to ensure PD in Cholesky:
    t = trace(S)/max(d,1);
    Sigma(:,:,k) = S + jitter*max(t,1e-12)*eye(d);
end

% -------- precompute Cholesky, log|Σ|
R = cell(K,1);
logDet = zeros(K,1);
for k = 1:K
    try
        R{k} = chol(Sigma(:,:,k));
    catch
        % if still not PD, nudge a bit more
        bump = 1e-6*trace(Sigma(:,:,k))/d;
        R{k} = chol(Sigma(:,:,k) + bump*eye(d));
        Sigma(:,:,k) = Sigma(:,:,k) + bump*eye(d);
    end
    logDet(k) = 2*sum(log(diag(R{k})));
end

% -------- helper: log posterior scores for many points
    function Slog = logPosterior(Z)
        % returns [N x K] log-scores up to an additive constant
        N = size(Z,1);
        Slog = zeros(N,K);
        for kk = 1:K
            XC = Z - mu(kk,:);
            Y  = XC / R{kk};            % (Z-mu)/R
            Q  = sum(Y.^2,2);           % Mahalanobis^2
            Slog(:,kk) = -0.5*(Q + logDet(kk)) + log(pi_k(kk));
        end
    end

% -------- helper: 2D Gaussian PDF (no toolbox)
twopi = 2*pi;
normConst = zeros(1,K);
for k = 1:K
    normConst(k) = 1/sqrt((twopi)^d * exp(logDet(k))); % = |Σ|^{-1/2} (2π)^{-d/2}
end
    function f = gaussPdf(Z, kidx)
        XC = Z - mu(kidx,:);
        Y  = XC / R{kidx};
        Q  = sum(Y.^2,2);
        f  = normConst(kidx) * exp(-0.5*Q);
    end

% -------- compute per-class error via overlap integral on a grid
Gaussian_error = zeros(K,1);

for k = 1:K
    % Conservative box using largest eigenvalue (rotationally safe)
    ev  = eig(Sigma(:,:,k));
    r   = trunc*sqrt(max(ev));                 % radius per axis
    x1r = linspace(mu(k,1)-r, mu(k,1)+r, res);
    x2r = linspace(mu(k,2)-r, mu(k,2)+r, res);
    [Xg, Yg] = meshgrid(x1r, x2r);
    Z = [Xg(:), Yg(:)];

    % Class pdf (k) over the box
    fk = gaussPdf(Z, k);

    % QDA prediction over the box
    scores = logPosterior(Z);
    [~, predIdx] = max(scores, [], 2);   % 1..K predictions

    % Mass correctly classified (approx), normalized by total mass in box
    % (tail beyond box is negligible for trunc>=5)
    massIn  = sum(fk(predIdx == k));
    massAll = sum(fk);
    pCorrect = massIn / massAll;
    Gaussian_error(k) = 1 - pCorrect;
end
end
