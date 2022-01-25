%% Example of the truncated Wirtinger Flow (TWF) algorithm under 1D Gaussian designs
% The TWF algorithm is presented in the paper
% ``Solving Random Quadratic Systems of Equations Is Nearly as Easy as Solving Linear Systems'' by Y. Chen and E. J. Candès.
% The code below is adapted from implementation of the Wirtinger Flow algorithm designed and implemented by E. Candes, X. Li, and M. Soltanolkotabi

%% Set Parameters
if exist('Params')                == 0,  Params.n2          = 1;    end
if isfield(Params, 'n1')          == 0,  Params.n1          = 1000; end             % signal dimension
if isfield(Params, 'm')           == 0,  Params.m           = 8* Params.n1;  end     % number of measurements
if isfield(Params, 'cplx_flag')   == 0,  Params.cplx_flag   = 0;    end             % real: cplx_flag = 0;  complex: cplx_flag = 1;
if isfield(Params, 'grad_type')   == 0,  Params.grad_type   = 'TWF_Poiss';  end     % 'TWF_Poiss': Poisson likelihood

if isfield(Params, 'alpha_lb')    == 0,  Params.alpha_lb    = 0.3;  end
if isfield(Params, 'alpha_ub')    == 0,  Params.alpha_ub    = 5;    end
if isfield(Params, 'alpha_h')     == 0,  Params.alpha_h     = 5;    end
if isfield(Params, 'alpha_y')     == 0,  Params.alpha_y     = 3;    end 
if isfield(Params, 'T')           == 0,  Params.T           = 500;  end    	% number of iterations
if isfield(Params, 'mu')          == 0,  Params.mu          = 0.2;  end		% step size / learning parameter
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 50;   end		% number of power iterations

n           = Params.n1;    
m           = Params.m;         
cplx_flag	= Params.cplx_flag;  % real-valued: cplx_flag = 0;  complex-valued: cplx_flag = 1;    
display(Params)
        
%% Make signal and data (noiseless)
x = randn(n,1)  + cplx_flag * 1i * randn(n,1); 
Amatrix = (randn(m,n) + cplx_flag * 1i * randn(m,n)) / (sqrt(2)^cplx_flag);
A  = @(I) Amatrix  * I;
At = @(Y) Amatrix' * Y;
y  = abs(A(x)).^2;
%y  = poissrnd(y);

%% Check results and Report Success/Failure
[Relerrs] = TWF(y, x, Params, A, At);
T = Params.T;
fprintf('Relative error after initialization: %f\n', Relerrs(1))
fprintf('Relative error after %d iterations: %f\n', T, Relerrs(T+1))
 
figure, semilogy(0:Params.T,Relerrs) 
xlabel('Iteration'), ylabel('Relative error (log10)'), ...
     title('Relative error vs. iteration count')
