%% Rotation Matrix

%% Description
% The function returns the rotiation matrix that gives the new set of
% coordinates as a function of the old ones
%% Calling Syntax
% function [matrix] = rotation_matrix(axis,angle)
%
%% I/O Variables
% |IN Integer|  *axis*: _axis of rotation_
%
% |IN Double|  *angle*: _angle of rotation (rad)_ $\theta$

%% Example
% axis = 1;
% angle = 30;
% 
% [matrix] = rotation_matrix(axis,angle);   

%%
function [matrix] = rotation_matrix(axis,angle)
    if axis == 1 || axis == 'X'
        matrix = [1,0,0;0,cos(angle),sin(angle);0,-sin(angle),cos(angle)];
    elseif axis == 2 || axis == 'Y'
        matrix = [cos(angle),0,-sin(angle);0,1,0;sin(angle),0,cos(angle)];
    elseif axis == 3 || axis == 'Z'
        matrix = [cos(angle),sin(angle),0;-sin(angle),cos(angle),0;0,0,1];
    end
end