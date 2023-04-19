function Rx(θ)
    return [1   0       0;
            0   cosd(θ)  sind(θ);
            0   -sind(θ)  cosd(θ)]
end

function Ry(θ)
    return [cosd(θ)  0  sind(θ);
            0       1  0;
            -sind(θ) 0  cosd(θ)]
end

function Rz(θ)
    return [cosd(θ)  sind(θ) 0;
            -sind(θ) cosd(θ) 0;
            0       0      1]
end