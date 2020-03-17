function [nbus, nl, ng, A, X, Pmax, Pmin, Fmax, Fmin, Pd, gens, buses, lines]  = ...
    loadData(mpc)

define_constants

gens = find(~isload(mpc.gen));
nbus  = size(mpc.bus, 1);                    % Number of Buses
buses = mpc.bus(:, 1);
ng    = size(mpc.gen, 1);                    % Number of Generators
Jn = zeros(ng, 1);
for j = 1:ng
    Jn(j) = mpc.gen(j, 1);
end

nl   = size(mpc.branch, 1);                 % Number of Lines
X    = mpc.branch(:,BR_X);                  % Reactance


Pmax = mpc.gen(:, PMAX);                   % Real power maximum
Pmin = zeros(ng, 1);                       % Real power minimum
Pd   = mpc.bus(:, PD); 

lines = [mpc.branch(:, F_BUS) mpc.branch(:, T_BUS)];

% Flow limits
Fmax = 100.*ones(nl, 1);
Fmin = -Fmax;

A = zeros(nbus, nl);
for i = 1:nbus
    for ll = 1:nl
        k = lines(ll, 1);
        l = lines(ll, 2);
        if i == k
            A(i,ll) = 1;
        elseif i == l
            A(i,ll) = -1;
        end
    end
end

end