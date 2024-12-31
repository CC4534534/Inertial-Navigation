function eulerAngles = C2Euler123(Cm)
    % Convert rotation matrix Cm to Euler angles (123 convention)
    
    % Calculate the Euler angles
    roll = atan2(Cm(2,3), Cm(3,3));        % ArcTan(Cm(2,3), Cm(3,3))
    pitch = asin(-Cm(1,3));                 % ArcSin(-Cm(1,3))
    yaw = atan2(Cm(1,2), Cm(1,1));          % ArcTan(Cm(1,2), Cm(1,1))
    
    % Return the Euler angles as a vector
    eulerAngles = [roll, pitch, yaw];
end
