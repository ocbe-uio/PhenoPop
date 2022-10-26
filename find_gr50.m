function [gr50_conc] = find_gr50(params,Conc) %params = alpha b E n
    min_concentration = Conc(1);
    max_concentration = Conc(length(Conc));
    % Find concentration c s.t. GRmax - GR(c) = 1/2 * (GRmax - GRmin)
    %fun = @(c) ratefunc(params,min_concentration) + ratefunc(params,max_concentration) - 2 * ratefunc(params,c);

    initial_guess = Conc(fix(length(Conc)/2));
    GR_minconc = ratefunc(params,min_concentration);
    GR_maxconc = ratefunc(params,max_concentration);

    syms x;
    eqn = GR_maxconc + GR_minconc - 2 * ratefunc(params,x) == 0;
    %eqn = simplify(eqn)
    foo = simplify(eqn);
    vpa_range = [0.0000001 max_concentration];
    gr50_conc = vpasolve(eqn, x, vpa_range);

    if numel(gr50_conc) ~= numel([1])
        "No GR50 value was found; Setting GR50 value to max_concentration"
        params;
        Conc;
        GR_minconc; 
        GR_maxconc;
        eqn;
        vpa_range;
        gr50_conc = max_concentration;
        % Plot left hand side
        %figure()
        %lhs = children(eqn);
        %lhs = lhs(1);
        %fplot(lhs, [min_concentration  max_concentration])
        %set(gca, 'XScale', 'log')
        %set(gca, 'YScale', 'log')
        %saveas(gcf, [pwd, '/article/find_gr50_lhs.png'])
    end

end 
