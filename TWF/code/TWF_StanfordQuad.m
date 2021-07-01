%% Implementation of the Truncated Wirtinger Flow (TWF) algorithm when
%% applied to real images

% The standalone code below builds on top of the implementation of the Wirtinger Flow algorithm designed and implemented by E. Candes, X. Li, and M. Soltanolkotabi

%% Set Parameters
if exist('Params')                == 0,  Params.L           = 12;   end
if isfield(Params, 'L')           == 0,  Params.L           = 12;   end             % number of measurements

if isfield(Params, 'alpha_lb')    == 0,  Params.alpha_lb    = 0.3;  end
if isfield(Params, 'alpha_ub')    == 0,  Params.alpha_ub    = 5;    end
if isfield(Params, 'alpha_h')     == 0,  Params.alpha_h     = 5;    end
if isfield(Params, 'alpha_y')     == 0,  Params.alpha_y     = 3;    end 
if isfield(Params, 'mu')          == 0,  Params.mu          = 0.2;  end		% step size / learning parameter
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 50;   end		% number of power iterations


%%  Read Image 

% Below X is n1 x n2 x 3; i.e. we have three n1 x n2 images, one for each of the 3 color channels  
namestr = 'stanford' ;
stanstr = 'jpg'      ;
X       = mat2gray(imread([namestr,'.',stanstr])) ;
n1      = size(X,1)                               ;
n2      = size(X,2)                               ;

%% Make masks and linear sampling operators  

% Each mask has iid entries following the octanary pattern in which the entries are 
% distributed as uniform over {1, -1, i, -i} (phase) 

randn('state',2014);
rand ('state',2014);

L = Params.L;                  % Number of masks  
Masks = zeros(n1,n2,L);  % Storage for L masks, each of dim n1 x n2
m = n1* n2* L;

% Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
for ll = 1:L, Masks(:,:,ll) = randsrc(n1,n2,[1i -1i 1 -1]); end

% Make linear operators; 
A = @(I)  fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is n1 x n2 image, output is n1 x n2 x L array
At = @(Y) sum(Masks .* ifft2(Y), 3) * size(Y,1) * size(Y,2);                       % Input is n1 x n2 x L array, output is n1 x n2 image

%% Prepare structure to save intermediate results 

ttimes   = [10:10:50];        % Iterations at which we will save info 
ntimes   = length(ttimes)+1;   % +1 because we will save info after the initialization 
Xhats    = cell(1,ntimes);
for mm = 1:ntimes, Xhats{mm} = zeros(size(X));  end
Times    = zeros(3,ntimes);

%% Truncated Wirtinger flow  

T = max(ttimes);                    % Max number of iterations
tau0 = 330;                         % Time constant for step size 
mu = @(t) Params.mu;                      % Schedule for step size 

for rgb = 1:3, 
    fprintf('Color band %d\n', rgb)
    x = squeeze(X(:,:,rgb)); % Image x is n1 x n2 
    Y = abs(A(x)).^2;        % Measured data 
    
    % Initialization
    normest = sqrt(sum(Y(:))/numel(Y));         % Estimate norm to scale eigenvector  
    z0 = randn(n1,n2); z0 = z0/norm(z0,'fro');  % Initial guess 
    tic
    % Truncated power iterations
    Ey = abs(Y) <= (Params.alpha_y)^2 * normest^2;
    Ytr = Y .* Ey;
    for tt = 1: Params.npower_iter, 
        z0 = At(Ytr.*A(z0)); z0 = z0/norm(z0,'fro');
    end
    Times(rgb,1) = toc;

    z = normest * z0;                   % Apply scaling 
    Xhats{1}(:,:,rgb) = exp(-1i*angle(trace(x'*z))) * z; % Initial guess after global phase adjustment 
 
    % Loop    
    fprintf('Done with initialization, starting loop\n')
    tic
    for t = 1:T,
        Bz = A(z);
        absBz = abs(Bz);

        normz = norm(z,'fro');
        hz_norm = 1/m* norm(absBz(:).^2 - Y(:), 1); 
        diff_Bz_Y = absBz.^2 - Y;   

        E    =  (absBz  <= Params.alpha_ub * normz) .* ...
                (absBz  >= Params.alpha_lb * normz) .* ...
                (abs(diff_Bz_Y) <= Params.alpha_h * hz_norm / normz * absBz);      
        C    = 2* (diff_Bz_Y) ./ conj(Bz)  .*  E;
        grad = At(C) / numel(C);                    % Wirtinger gradient
        z    = z - mu(t) * grad;  % Gradient update 
        
        ind =  find(t == ttimes);                % Store results 
        if ~isempty(ind), 
             Xhats{ind+1}(:,:,rgb) = exp(-1i*angle(trace(x'*z))) * z; 
             Times(rgb,ind+1) = toc;
        end       
    end    
end
fprintf('All done!\n')

%% Show some results 

iter = [0 ttimes];
Relerrs = zeros(1,ntimes);
for mm = 1:ntimes; 
    fprintf('Mean running times after %d iterations: %.1f\n', iter(mm), mean(Times(:,mm)))
    Relerrs(mm) = norm(Xhats{mm}(:)-X(:))/norm(X(:)); 
    fprintf('Relative error after %d iterations: %f\n', iter(mm), Relerrs(mm))  
    fprintf('\n')
end

for tt = 1:ntimes, 
    figure; imshow(mat2gray(abs(Xhats{tt})),[]);
end
