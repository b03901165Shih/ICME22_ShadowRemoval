function [angles, magnit] = gradient_profile(gray_img)

% Create the operators for computing image derivative at every pixel.
hx = 0.5*[-1,0,1];
hy = hx';

% Compute the derivative in the x and y direction for every pixel.
dx = filter2(hx, padarray(double(gray_img), [1,1]));
dy = filter2(hy, padarray(double(gray_img), [1,1]));

% Remove the 1 pixel border.
dx = dx(2 : (size(dx, 1) - 1), 2 : (size(dx, 2) - 1));
dy = dy(2 : (size(dy, 1) - 1), 2 : (size(dy, 2) - 1)); 

% Convert the gradient vectors to polar coordinates (angle and magnitude).
angles = atan2(dy, dx);
magnit = ((dy.^2) + (dx.^2)).^.5;
end

