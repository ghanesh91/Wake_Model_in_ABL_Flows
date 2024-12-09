function [x, y] = computeCircle(centerX, centerY, radius, numPoints)
    % Compute the angles for plotting
    theta = linspace(0, 2*pi, numPoints);
    
    % Compute the x and y coordinates of the circle
    x = radius * cos(theta) + centerX;
    y = radius * sin(theta) + centerY;
end