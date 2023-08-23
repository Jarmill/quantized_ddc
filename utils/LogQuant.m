%from the codegen routine in https://www.mathworks.com/help/control/ug/absolute-stability-for-quantized-system.html
function y = LogQuant(u,rho)

% Parameters for the quantizer
u0 = 1;

% Logarithmic quantizer
if (u == 0)
    y = 0;
elseif (u > 0)
        y = LogQuantPositive(u,u0,rho);
    else
        y = - LogQuantPositive(-u,u0,rho);
    end  
end

