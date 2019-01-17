function F_n = calc_spring_forces(R_n,V_n,F_n,KS,KD,springs)
    for i=1:length(springs)
        r = R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:);
        d = norm(r);
        v = V_n(springs(i).fromPI,:) - V_n(springs(i).toPI,:);
        f = ( (KS*(d-springs(i).L) + KD*dot(v,r)/d)/d)*r; % f_ba not f_ab
        F_n(springs(i).toPI,:) = F_n(springs(i).toPI,:) + f;
        F_n(springs(i).fromPI,:) = F_n(springs(i).fromPI,:) - f;
    end
end

