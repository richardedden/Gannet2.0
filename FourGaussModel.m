function F = FourGaussModel(x,freq)

% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% Peak 2 has same width as 1, 0.66*height and shift+0.19ppm
% x(4-6) = lactate
% x(7) = offset
% x(8) = slope
% x(9) = quadratic

F = x(1) * exp(x(2) * (freq-x(3)).^2) + ...
    0.66*x(1) * exp(x(2) * (freq-(x(3)+0.19)).^2) + ...
    x(4) * exp(x(5) * (freq-x(6)).^2) + ...
    x(4) * exp(x(5) * (freq-(x(6)-0.055)).^2) + ...
    x(7) + x(8)*(freq-x(6)) + x(9)*(freq-x(6)).^2;