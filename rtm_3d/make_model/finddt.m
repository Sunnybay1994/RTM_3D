function dtmax = finddt(epr_min, miur_min, dx, dy, dz)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmin = epr_min * ep0;
    mumin = miur_min * mu0;
    dtmax = 6.0 / 7.0 * sqrt(epmin * mumin / (dx^-2 + dy^-2 + dz^-2));
end

