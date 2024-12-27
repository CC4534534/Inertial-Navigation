function Euler3212C = Euler3212C(q)
    % Calculate the rotation matrix using a 3-2-1 Euler angle sequence
    % Input: q - A vector of 3 Euler angles [q1, q2, q3]
    % Output: C - The 3x3 rotation matrix

    % Precompute trigonometric values
    st1 = sin(q(1));
    ct1 = cos(q(1));
    st2 = sin(q(2));
    ct2 = cos(q(2));
    st3 = sin(q(3));
    ct3 = cos(q(3));

    % Construct the rotation matrix
    Euler3212C = [
        ct2 * ct1,               ct2 * st1,               -st2;
        st3 * st2 * ct1 - ct3 * st1, st3 * st2 * st1 + ct3 * ct1, st3 * ct2;
        ct3 * st2 * ct1 + st3 * st1, ct3 * st2 * st1 - st3 * ct1, ct3 * ct2
    ];
end