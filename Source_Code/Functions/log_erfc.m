function[d] = log_erfc(x)

d     = nan(size(x));
indx_high  = x>27.22;
indx_low   = x<=27.22;
d(indx_low)     = log(erfc(x(indx_low)));
x_high = x(indx_high);
% d(indx_high) = -x(indx_high).^2 - log(x(indx_high)) + log(sqrt(pi) + 1./(2*x(indx_high)));  % o(x^-2)
% d(indx_high) = -x_high.^2 - log(x_high) - 1./(2*x_high.^2) + log(sqrt(pi) + 1./(2*x_high)); % o(x^-4)
% d(indx_high) = -x_high.^2 - log(x_high) - 1./(2*x_high.^2) - 3./(4*x_high.^4) + log(sqrt(pi) + 1./(2*x_high)); % o(x^-6)
d(indx_high) = -x_high.^2 - log(x_high) - 1./(2*x_high.^2) - 3./(4*x_high.^4) - 15./(8*x_high.^6) - 105./(16*x_high.^8)+ log(sqrt(pi) + 1./(2*x_high)); % o(x^-11)


