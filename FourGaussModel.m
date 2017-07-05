function F = FourGaussModel(x,freq)

% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% Gaussian 2 has same width as 1, 0.75 height and shift+0.19ppm
% x(5-7) = Lactate
% x(8) = offset
% x(9) = slope
% x(10) = quadratic

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    0.66*x(1)*exp(x(2)*(freq-(x(3)+0.19)).*(freq-(x(3)+0.19))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) +...
    x(4)*exp(x(5)*(freq-(x(6)-0.055)).*(freq-(x(6)-0.055))) +...
    x(7) + x(8)*(freq - x(6)) + x(9)*(freq-x(6)).^2 ;
