function [r_f, r_delta, C_S_surf, solutionType] =...
    solveSteadyStateDimensionless3(c1, c2, c3)


%build the r_delta
R_delta_fun = @(R_f)(sqrt(R_f.^2 + c1./R_f.^2));


%build function of first solution of C at surf
C_S_surf_1 =...
    @(R_f, R_delta)...
    (1/c2*(R_delta.^2.*(1-log(R_delta./R_f)) - R_f.^2));


%build function of second solution of C at surf
C_S_surf_2 =...
    @(R_f, R_delta)...
    1 - c3 / c2 *...
    (R_delta.^2 - R_f.^2);


% build equation to solve
f_to_solve = @(R_f, R_delta)...
    (C_S_surf_2(R_f, R_delta) - C_S_surf_1(R_f, R_delta));
f_to_solve_R_f = @(R_f)(f_to_solve(R_f, R_delta_fun(R_f)));

% compute the solutions
try
    r_f = fzero(f_to_solve_R_f, 1);
catch
    warning ('invalid solution for r_f');
    r_f = NaN;
end
r_delta = R_delta_fun(r_f);
C_S_surf = C_S_surf_2(r_f, r_delta);

% determine if this is valid solution
if isnan(r_f)
    solutionType = -4; %not solved
    return;
elseif (C_S_surf < 0)
    warning ('invalid solution C_S_surf_2 < 0');
    r_f = NaN;
    r_delta = NaN;
    C_S_surf = NaN;
    solutionType = -3; %biofilm thickness less than R
    return;
elseif (r_f > 1)
    warning ('invalid solution r_f > 1');
    r_f = NaN;
    r_delta = NaN;
    C_S_surf = NaN;
    solutionType = -2; %biofilm thickness less than R
    return;
elseif (r_f > r_delta)
    warning ('invalid solution r_f > r_delta');
    r_f = NaN;
    r_delta = NaN;
    C_S_surf = NaN;
    solutionType = -1; %biofilm thickness less than delta
    return;
elseif (r_delta > 1)
    warning ('fully penetrated biofilm');
    r_delta = 1;
    C_S_surf = C_S_surf_2(r_f, r_delta);
    solutionType = 1; %fully penetrated biofilm
    return;
end
solutionType = 2; %valid


