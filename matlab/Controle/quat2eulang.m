function [Z1,Y1,Z2] = quat2eulang(q0,q1,q2,q3,seq)
    %% Quaternion to Euler Angles (DCM)
    L = length(q2);
    q = [ones(L,1)*q0,q1,q2,q3];
    [Z1,Y1,Z2] = quat2angle(q,seq);
end