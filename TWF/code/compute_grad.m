%%  Compute the truncated gradient based on the Poisson log-likelihood function

function grad = compute_grad(z, y, Params, A, At)
    m = Params.m;
    yz = A(z);
    Kt = 1/m* norm(abs(yz(:)).^2 - y(:), 1); 
        
    if strcmp(Params.grad_type,'TWF_Poiss') == 1   % truncated gradient / Wirtinger flow
        % truncation rules
        Eub =  abs(yz) / norm(z)   <= Params.alpha_ub;
        Elb =  abs(yz) / norm(z)   >= Params.alpha_lb;
        Eh  =  abs(y - abs(yz).^2) <= Params.alpha_h * Kt / norm(z) * abs(yz);
        
        grad  = 1/m* At( 2* ( abs(yz).^2-y ) ./ (abs(yz).^2) .*yz ...
                          .* Eub .* Elb .* Eh );    % truncated Poisson gradient
                      
    elseif strcmp(Params.grad_type,'WF_Poiss') == 1    % untruncated gradient / Wirtinger flow
        grad  = 1/m* At( 2* ( abs(yz).^2-y ) ./ (abs(yz).^2) .*yz ); % Poisson gradient
    end
      
