function fitparams = fitirfgamma(xdata, data, initparams)
%
% fitparams = fitexpmat(xdata, data, initparams)
%
% Fits 'data' to exponential curve y=A*exp(-x/tau)+y0 using initial
% parameters in initparams. 
%
% xdata is a column vector of x values for the data.
%
% 'data' is a colunm vector or a matrix of column vectors of y data points
% to be fit. 
%
% 'initparams' is 3xn matrix of initial parameters, where n is the number
% of columns in 'data'.  Each column of initparams = [tau; A; y0]. 
%
% 'fitparams' has the same format as 'initparams' and are the best fit
% results of each column in 'data'.  
% 
%  


    % nested function call defining chi-square calculation for a single
    % exponential.
    function [chisqu] = chisqu_exp(params)
        % pull out fitting parameters
        x0 = params(1);
        p = params(2)^2;
        c = params(3)^2;
        A = params(4);
%        y0 = params(5);
        y0=0;
        
        xrightid1=find(xdata > x0,1);

        % compute fitted exponential, residual and chi-square.
        fit = A/gamma(p)*[zeros(xrightid1-1,1); exp(-(xdata(xrightid1:end)'-x0)*c).*((xdata(xrightid1:end)'-x0)*c).^(p-1)] + y0;
        residual = onedata - fit;
        chisqu = sum(residual.^2);
    end % chisqu_exp

[numx numfits] = size(data);


% error checking
if any(size(initparams) ~= [4 numfits])
    error('initparams must be a 4-by-n matrix');
end

fitparams = [];
model = @chisqu_exp; 
for i=1:numfits
    
    onedata = data(:,i);
    
    % chi-square minimization call
    options=optimset('MaxFunEval', 100000, 'MaxIter', 100000);    
    estimates = fminsearch(model, initparams(:,i), options);



    fitparams = [fitparams estimates];
end % for loop

end % fitexpmat


