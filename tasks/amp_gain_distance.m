function [wavefield,factor1,factor2] = amp_gain_distance(wavefield0,src_pos,x,y,z,mul1,mute_angle,mul2)
%amp_gain_distance Amplitude gain according to distance from source
%   Assume a Spherical propagation
    disp('amp_gain_distance')
    [nz,ny,nx]=size(wavefield0);
    srcx = src_pos(1);
    srcy = src_pos(2);
    srcz = src_pos(3);
    if nargin < 3
        x = 1:nx;
        y = 1:ny;
        z = 1:nz;
    end
    if nargin < 6
        mul1 = 2;
    end
    if nargin < 7
        mute_angle = pi/6;
    end
    if nargin < 8
        mul2 = 0;
    end
    factor1 = zeros(size(wavefield0));
    factor2 = ones(size(wavefield0));
    dist0 = (z(end)-srcz)*1/3;
    for i=1:nx
        ix = x(i);
        for j=1:ny
            iy = y(j);
            for k=1:nz
                iz = z(k);
                zdist = abs(iz-src_pos(3));
%                 zdist = norm([ix,iy,iz]-src_pos);
                factor1(k,j,i) = (zdist/dist0)^mul1;
                pos_angle = atan(norm([ix,iy]-[srcx,srcy])/(iz-srcz));
                if pos_angle > mute_angle && pos_angle >= 0
                    factor2(k,j,i) = ((pi/2-pos_angle)/(pi/2-mute_angle))^mul2;
                elseif pos_angle < 0
                    factor2(k,j,i) = 0;
                end
            end
        end
    end
    wavefield = wavefield0 .* factor1 .*factor2;
end

