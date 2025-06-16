%initializations
rng default

%generates 1000 samples of initial guesses for model parameters
initial_guess =  lhsdesign(5,500);

Output_param = [];

% for parallel computing
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

% for parallel computing
parfor i=1:500
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
writematrix(Output_param,'ControlModelOutput.csv');
smallest_optimised_param_costfunction = min(Output_param(6,:));
disp(smallest_optimised_param_costfunction)
col_num_of_smallest_op_param = find(Output_param(6,:)==smallest_optimised_param_costfunction);
optimised_param = Output_param(:, col_num_of_smallest_op_param).'; %this matrix would include the cost function at the end for reference
disp(optimised_param)


%runs globalsearch for each initial parameter guess
function r = runGS(p0)
    
    %lower and upper bounds for model parameter search
    lb = [0 7 0 0 0];
    ub = [100000 14 100000 100000 1000000];
    problem = createOptimProblem('fmincon','objective',@myobjective,'x0',p0,'lb',lb,'ub',ub);
    gs = GlobalSearch('PlotFcn',[]);
    r = run(gs,problem);
end

%defines the objective function
function SSE = myobjective(p)
PBS = readtable("VDNE control.csv");
D8 = PBS.D(PBS.Day==8); D8=D8(~isnan(D8));
D15 = PBS.D(PBS.Day==15); D15=D15(~isnan(D15));
D21 = PBS.D(PBS.Day==21); D21=D21(~isnan(D21));
D20 = PBS.D(PBS.Day==20); D20=D20(~isnan(D20));
D25 = PBS.D(PBS.Day==25); D25=D25(~isnan(D25));
E8 = PBS.E(PBS.Day==8); E8=E8(~isnan(E8));
E15 = PBS.E(PBS.Day==15); E15=E15(~isnan(E15));
E21 = PBS.E(PBS.Day==21); E21=E21(~isnan(E21));
E20 = PBS.E(PBS.Day==20); E20=E20(~isnan(E20));
E25 = PBS.E(PBS.Day==25); E25=E25(~isnan(E25));

yd = [median(D8) median(E8);
      median(D15) median(E15);
      median(D21) median(E21);
      median(D20) median(E20);
      median(D25) median(E25)]; %ydata [3; 10; 16; 5; 10]

%considering weight for each error to be inverse of the variance in
%experimental values
weights = [1/var(yd(:,1)) 1/var(yd(:,2))];
td = [0; 3; 10; 16; 0; 5; 10];

%ys has ODE solutions for both infections
ys = solveODE(p,td,[[0 0];...
    [yd(2,1) yd(2,2)]]);

SSE = sum(sum([weights(1)*((yd(:,1)-ys(:,1)).^2), weights(2)*((yd(:,2)-ys(:,2)).^2)]));
end

function C = solveODE(p,td,y0)
[~,Cv1]=ode23s(@DifEq,td(1:4),y0(1,:)); %first infection
[~,Cv2]=ode23s(@DifEq,td(5:7),y0(2,:)); %second infection
    function dY = DifEq(td,y)
        dydt = zeros(2,1);
        
        dydt(1) =  p(1)*y(2)*(p(2)-y(1)) - p(3)*y(1); %D(t)
        
        dydt(2) = p(4) - p(5)*y(2); %E(t)

        dY = dydt;
    end
C = [Cv1(2:4,:); Cv2(2:3,:)]; %merging two solutions in one matrix for optimisation
end