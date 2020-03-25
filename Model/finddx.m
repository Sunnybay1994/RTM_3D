function dxmax = finddx(epr_max, miur_max, fmax)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmax = epr_max * ep0;
    mumax = miur_max * mu0;
    wlmin = 1 / (fmax * sqrt(epmax * mumax));
    dxmax = wlmin;
end

