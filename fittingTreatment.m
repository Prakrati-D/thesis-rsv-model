%initializations
rng default

%generates 1000 samples of initial guesses for model parameters
initial_guess =  lhsdesign(7,1000);

%day 3, 10 and 16 are from first infection, day 5 and 10 from second
%infection
td = [3; 10; 16; 5; 10];

Output_param = [];

% for parallel computing
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

% for parallel computing
parfor i=1:1000
    %testing one parameter set at a time from 1000 samples
    p0 = initial_guess(:,i);

    %calls the globalsearch algorithm
    param = runGS(p0);

    %stores optimised parameter set for each iteration and corresponding
    %cost function value
    if(size(param)~=0)
        [output] = param;
        cost_fcn = myobjective(param);
        [output] = [output;cost_fcn]
        Output_param = [Output_param output];
    end
end

% collects smallest cost function value and corresponding optimised
% parameter set
writematrix(Output_param,'TreatmentModelOutput.csv');
smallest_optimised_param_costfunction = min(Output_param(8,:));
disp(smallest_optimised_param_costfunction)
col_num_of_smallest_op_param = find(Output_param(8,:)==smallest_optimised_param_costfunction);
optimised_param = Output_param(:, col_num_of_smallest_op_param).'; %this matrix would include the cost function at the end for reference
disp(optimised_param)


%runs globalsearch for each initial parameter guess
function r = runGS(p0)
    
    %lower and upper bounds for model parameter search
    lb = [0 0 0 0 0 0 0];
    ub = [100000 100000 100000 1000000 100000 100000 1];
    problem = createOptimProblem('fmincon','objective',@myobjective,'x0',p0,'lb',lb,'ub',ub);
    gs = GlobalSearch('PlotFcn',[]);
    r = run(gs,problem);
end

%defines the objective function
function SSE = myobjective(p)

RSV = readtable("VDNE treatment.csv");
D8 = RSV.D(RSV.Day==8); D8=D8(~isnan(D8));
D15 = RSV.D(RSV.Day==15); D15=D15(~isnan(D15));
D21 = RSV.D(RSV.Day==21); D21=D21(~isnan(D21));
D20 = RSV.D(RSV.Day==20); D20=D20(~isnan(D20));
D25 = RSV.D(RSV.Day==25); D25=D25(~isnan(D25));
E8 = RSV.E(RSV.Day==8); E8=E8(~isnan(E8));
E15 = RSV.E(RSV.Day==15); E15=E15(~isnan(E15));
E21 = RSV.E(RSV.Day==21); E21=E21(~isnan(E21));
E20 = RSV.E(RSV.Day==20); E20=E20(~isnan(E20));
E25 = RSV.E(RSV.Day==25); E25=E25(~isnan(E25));
V8 = RSV.V(RSV.Day==8); V8=V8(~isnan(V8));
V15 = RSV.V(RSV.Day==15); V15=V15(~isnan(V15));
V21 = RSV.V(RSV.Day==21); V21=V21(~isnan(V21));
V20 = RSV.V(RSV.Day==20); V20=V20(~isnan(V20));
V25 = RSV.V(RSV.Day==25); V25=V25(~isnan(V25));

yd = [median(V8) median(D8) median(E8);
      median(V15) median(D15) median(E15);
      median(V21) median(D21) median(E21);
      median(V20) median(D20) median(E20);
      median(V25) median(D25) median(E25)]; %ydata [3; 10; 16; 5; 10]

%considering weight for each error to be inverse of the variance in
%experimental values
weights = [1/var(yd(:,1)) 1/var(yd(:,2)) 1/var(yd(:,3))];
td = [0; 3; 10; 16; 0; 5; 10];

%ys has ODE solutions for both infections
ys = solveODE(p,td,[[p(7) 0 0];...
    [p(7) median(D15) median(E15)]]);

SSE = sum(sum([weights(1)*((yd(:,1)-ys(:,1)).^2), weights(2)*((yd(:,2)-ys(:,2)).^2), weights(3)*((yd(:,3)-ys(:,3)).^2)]));
end

function C = solveODE(p,td,y0)
    deltaBE = 0;
    kappaB = 0;
    kappaE = 0;
    deltaE = 0;
    Dmax = 0;
    flag=0;
    [~,Cv1]=ode23s(@DifEq,td(1:4),y0(1,:)); %first infection
    [~,Cv2]=ode23s(@DifEq,td(5:7),y0(2,:)); %second infection
    function dY = DifEq(td,y)
        dydt = zeros(3,1);
        
        dydt(1) = p(1)*y(1)*(1+p(2)*y(2)) - p(3)*y(1); %V(t)
        
        dydt(2) =  p(4)*y(1)*(Dmax-y(2)) + deltaBE*y(3)*(Dmax-y(2)) - kappaB*y(2); %D(t)
        
        if (y(2) < p(5)) && flag==0
            dydt(3) = kappaE - deltaE*y(3);
        elseif (y(2) >= p(5)) || flag==1
            %setting flag=1 is changing the switch
            flag=1;
            dydt(3) = kappaE + p(6)*y(2) - deltaE*y(3); %E(t)
        end
        
        dY = dydt;
    end
    C = [Cv1(2:4,:); Cv2(2:3,:)]; %merging two solutions in one matrix for optimisation
end