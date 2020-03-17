define_constants

%mpc = loadcase('case5arroyoB');
mpc = loadcase('case24_ieee_rts');
mpc = loadcase('case30');
%mpc = loadcase('case57');
%mpc = loadcase('case118');

[nbus, nl, ng, A, X, Pmax, Pmin, Fmax, Fmin, Pd, gens, buses, lines] =...
    loadData(mpc);
    
%% Parameters
delL   = -pi/2;
delH   = pi/2;
muL    = -20000;
muH    = 20000;
%% Primal Variables
v      = binvar(nl,   1);        % Interdiction variable
s      = sdpvar(nbus, 1);        % Load Shed
Pg     = sdpvar(ng, 1);        % Generation at generator g
Pf     = sdpvar(nl,   1);        % Flow on line l
delta  = sdpvar(nbus, 1);        % Nodal phase angle

t      = sdpvar(nl,   1);        % Variable equal to u*v
h      = sdpvar(nl,   1);        % Slack Variable to linearize u*v
sF     = sdpvar(nl,   1);        % Slack Variable to linearize v*deltaFrom
sT     = sdpvar(nl,   1);        % Slack Variable to linearize v*deltaTo
zF     = sdpvar(nl,   1);        % Variable equal to v*deltaFrom
zT     = sdpvar(nl,   1);        % Variable equal to v*deltaTo
%% Dual Multipliers
mu     = sdpvar(nl,   1);        % Relates power flow and phase angles

lam    = sdpvar(nbus, 1);        % Power Balance

omH    = sdpvar(nl,   1);        % Line capacities
omL    = sdpvar(nl,   1);        % Line Capacities 

thetaH = sdpvar(ng,   1);        % Generator Capcities
thetaL = sdpvar(ng,   1);        % Generator Capcities

alphaH = sdpvar(nbus, 1);        % Max Curtailment
alphaL = sdpvar(nbus, 1);        % Min Curtailment

%% Compslackness vars
wOmH    = binvar(nl,   1);        
wOmL    = binvar(nl,   1);   

wThH    = binvar(ng,   1);        
wThL    = binvar(ng,   1);         

wAlH    = binvar(nbus, 1);         
wAlL    = binvar(nbus, 1);        

M = 2000;
%% Create the set of incident generators for each node
genLoc = cell(1, ng);
inGens = cell(1, nbus);

% Get generator locations
for i = 1:ng
    ns = mpc.gen(i,1);
    genLoc{1, i} = unique(ns);
end

% Get generators located at each bus
for i = 1:nbus
    ns = [];
    for j = 1:ng
        if mpc.gen(j, GEN_BUS) == i
            ns = [ns, j];
        end
    end
    inGens{1, i} = ns;
end 
%% Constraints
Sbar = .25*sum(Pd);
Curtailment = [sum(s) >= Sbar];
 
OPF           = [];
CompSlackness = [];

for l = 1:nl
    fb = lines(l, 1);   % Fbus for line l
    tb = lines(l, 2);   % Tbus for line l
    OPF = [OPF, Pf(l) == 1./X(l).*(zF(l) - zT(l)),...                % (15)
                zF(l) == delta(fb) - sF(l),...                      % (16)
                zT(l) == delta(tb) - sT(l),...                      % (17)
                delL*v(l) <= zF(l) <= delH*v(l),...                 % (18)
                delL*v(l) <= zT(l) <= delH*v(l),...                 % (19)
                delL*(1-v(l)) <= sF(l) <= delH*(1-v(l)),...         % (20)
                delL*(1-v(l)) <= sT(l) <= delH*(1-v(l)),...         % (21)
                -Fmax(l) <= Pf(l) <= Fmax(l),...                    % (23)
                t(l) == mu(l) - h(l),...                            % (27)
                muL*v(l) <= t(l) <= muH*v(l),...                    % (28)
                muL*(1-v(l)) <= h(l) <= muH*(1-v(l)),...            % (29)
                A(:, l)'*lam - mu(l) - omL(l) + omH(l) == 0];       % (31)      
end

for i = 1:nbus
    OPF = [OPF, 0 <= s(i) <= Pd(i),...                              % (25)
                1 - lam(i) - alphaL(i) + alphaH(i) == 0,...         % (32)
                (1./X)'.*A(i, :)*t == 0];                           % (26)
    incGens = inGens{1, i};
    if isempty(incGens) == 0 
        OPF = [OPF, sum(Pg(incGens))  - A(i,:)*Pf + s(i) == Pd(i)]; % (22)
    end
end

for i = 1:ng
    OPF = [OPF, Pmin(i) <= Pg(i) <= Pmax(i)] ;                        % (24)
end

for j = 1:ng
    Nn = genLoc{1, j};
    for nn = 1:length(Nn)
        jn = Nn(nn);
        OPF = [OPF,  -lam(jn) - thetaL(j) +  thetaH(j) == 0];
    end
end

Bounds = [omL >= 0, omH >= 0, thetaL >= 0, thetaH >= 0, alphaL >= 0,...
          alphaH >= 0];
      
for l = 1:nl
    CompSlackness = [CompSlackness, omL(l) <= M*wOmL(l),...
                                    Pf(l) + Fmax(l) <= M*(1-wOmL(l)),...
                                    omH(l) <= M*wOmH(l),...
                                    Fmax(l) - Pf(l) <= M*(1-wOmH(l)),...
                                    wOmL(l) + wOmH(l) <= 1];
end

for j = 1:ng
    CompSlackness = [CompSlackness, thetaL(j) <= M*wThL(j),...
                                    Pg(j) - Pmin(j) <= M*(1-wThL(j)),...
                                    thetaH(j) <= M*wThH(j),...
                                    Pmax(j) - Pg(j) <= M*(1-wThH(j)),...
                                    wThL(j) + wThH(j) <= 1];
end

for n = 1:nbus
    CompSlackness = [CompSlackness, alphaL(n) <= M*wAlL(n),...
                                    s(n) <= M*(1-wAlL(n)),...
                                    alphaH(n) <= M*wAlH(n),...
                                    Pd(n) - s(n) <= M*(1-wAlH(n)),...
                                    wAlL(n) + wAlH(n) <= 1];
end
%% Optimize
Constraints = [Curtailment, OPF, CompSlackness, Bounds];

Objective = -sum(v);

options = sdpsettings('verbose',1,'solver','gurobi');
%options.gurobi.TimeLimit = 60*30;                             % Time limit = 30mins
sol = optimize(Constraints,Objective,options);
vopt = value(v)
sopt = value(s)
Stot = sum(sopt)
time = sol.solvertime

 