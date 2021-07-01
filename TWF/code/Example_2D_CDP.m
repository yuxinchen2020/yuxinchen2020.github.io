%% Example of the truncated Wirtinger Flow (TWF) algorithm under 2D Code Diffraction Patterns (CDP)
% The TWF algorithm is presented in the paper
% ``Solving Random Quadratic Systems of Equations Is Nearly as Easy as Solving Linear Systems'' by Y. Chen and E. J. Candès.
% The code below is adapted from implementation of the Wirtinger Flow algorithm designed and implemented by E. Candes, X. Li, and M. Soltanolkotabi

%% Set Parameters
if exist('Params')                == 0,  Params.n1          = 100;  end
if isfield(Params, 'n2')          == 0,  Params.n2          = 100;  end             % signal dimension
if isfield(Params, 'L')           == 0,  Params.L           = 12;   end             % number of measurements
if isfield(Params, 'grad_type')   == 0,  Params.grad_type   = 'TWF_Poiss';  end     % 'TWF_Poiss': Poisson likelihood

if isfield(Params, 'alpha_lb')    == 0,  Params.alpha_lb    = 0.3;  end
if isfield(Params, 'alpha_ub')    == 0,  Params.alpha_ub    = 5;    end
if isfield(Params, 'alpha_h')     == 0,  Params.alpha_h     = 5;    end
if isfield(Params, 'alpha_y')     == 0,  Params.alpha_y     = 3;    end 
if isfield(Params, 'T')           == 0,  Params.T           = 500;  end    	% number of iterations
if isfield(Params, 'mu')          == 0,  Params.mu          = 0.2;  end		% step size / learning parameter
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 50;   end		% number of power iterations

n1          = Params.n1;    
n2          = Params.n2; 
L           = Params.L;         
display(Params)
        
%% Make signal and data (noiseless)
x       = randn(n1,n2)  + 1i * randn(n1,n2); 
Masks   = zeros(n1,n2,L);  % Storage for L masks, each of dim n1 x n2
m       = n1* n2* L;  Params.m = m;

% Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
for ll = 1:L, Masks(:,:,ll) = randsrc(n1,n2,[1i -1i 1 -1]); end

% Make linear operators; 
A = @(I)  fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is n1 x n2 image, output is n1 x n2 x L array
At = @(Y) sum(Masks .* ifft2(Y), 3) * size(Y,1) * size(Y,2);                       % Input is n1 x n2 x L array, output is n1 x n2 image

Y = abs(A(x)).^2; 
%Y  = poissrnd(Y);

%% Check results and Report Success/Failure
[Relerrs] = TWF(Y, x, Params, A, At);
T = Params.T;
fprintf('Relative error after initialization: %f\n', Relerrs(1))
fprintf('Relative error after %d iterations: %f\n', T, Relerrs(T+1))
 
figure, semilogy(0:Params.T,Relerrs) 
xlabel('Iteration'), ylabel('Relative error (log10)'), ...
     title('Relative error vs. iteration count')
