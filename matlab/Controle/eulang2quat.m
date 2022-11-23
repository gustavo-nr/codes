function q = eulang2quat(ang1,ang2,ang3,rot_seq)
    %% Euler Angles to Direct Cosine Matrix (DCM)
    DCM1 = rotation_matrix(rot_seq(1),ang1);
    DCM2 = rotation_matrix(rot_seq(2),ang2);
    DCM3 = rotation_matrix(rot_seq(3),ang3);
    DCM = DCM3*DCM2*DCM1;
    
    %% DCM to Quaternion - Sheppard's Method
    % For this step, it's used the Sheppard's Method
    %%
    % 1. Find the largest value between $q_0^2,\,q_1^2,\,q_2^2, \,q_3^2$
    q0_2 = 1/4*(1 + trace(DCM));
    q1_2 = 1/4*(1 + 2*DCM(1,1) - trace(DCM));
    q2_2 = 1/4*(1 + 2*DCM(2,2) - trace(DCM));
    q3_2 = 1/4*(1 + 2*DCM(3,3) - trace(DCM));
    q_2 = [q0_2 q1_2 q2_2 q3_2];
    %%
    % 2. Compute the remaining EPs using
    [M,I] = max(q_2);
    if I == 1
        q0 = sqrt(q0_2);
        q1 = (DCM(2,3)-DCM(3,2))/(4*q0);
        q2 = (DCM(3,1)-DCM(1,3))/(4*q0);
        q3 = (DCM(1,2)-DCM(2,1))/(4*q0);
    end
    if I == 2
        q1 = sqrt(q1_2);
        q0 = (DCM(2,3)-DCM(3,2))/(4*q1);
        q2 = (DCM(1,2)+DCM(2,1))/(4*q1);
        q3 = (DCM(3,1)+DCM(1,3))/(4*q1);
    end
    if I == 3
        q2 = sqrt(q2_2);
        q0 = (DCM(3,1)-DCM(1,3))/(4*q2);
        q1 = (DCM(1,2)+DCM(2,1))/(4*q2);
        q3 = (DCM(2,3)+DCM(3,2))/(4*q2);
    end
    if I == 4
        q3 = sqrt(q3_2);
        q0 = (DCM(1,2)-DCM(2,1))/(4*q3);
        q1 = (DCM(3,1)+DCM(1,3))/(4*q3);
        q2 = (DCM(2,3)+DCM(3,2))/(4*q3);
    end
    q = [q0 q1 q2 q3];
end
