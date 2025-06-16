%initializations
rng default

lb = [0.015 4E-5 0.087 2.04E4 8 1E-7 0.025 6 2.2 0.09 2.2 4.5E-7];
ub = [1.5 4E-3 8.7 2.04E6 16 1E7 2.5 10 220 9 220 4.5E-5];

%generates 1000 samples of initial guesses for model parameters
[x,initial_guess] =  lhsdesign_modifie(1,lb,ub);

Output_param = [];

% for parallel computing
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

% for parallel computing
parfor i=1:1
    %testing one parameter set at a time from 500 samples
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
writematrix(Output_param,'RawModelOutput.csv');
smallest_optimised_param_costfunction = min(Output_param(13,:));
disp(smallest_optimised_param_costfunction)
col_num_of_smallest_op_param = find(Output_param(13,:)==smallest_optimised_param_costfunction);
optimised_param = Output_param(:, col_num_of_smallest_op_param).'; %this matrix would include the cost function at the end for reference
disp(optimised_param)


%runs globalsearch for each initial parameter guess
function r = runGS(p0)
    
    %lower and upper bounds for model parameter search
    
    lb = [0.015 4E-5 0.087 2.04E4 8 1E-7 0.025 6 2.2 0.09 2.2 4.5E-7];
    ub = [1.5 4E-3 8.7 2.04E6 16 1E7 2.5 10 220 9 220 4.5E-5];
    problem = createOptimProblem('fmincon','objective',@myobjective,'x0',p0,'lb',lb,'ub',ub);
    gs = GlobalSearch('PlotFcn',[]);
    r = run(gs,problem);
end

%defines the objective function
function SSE = myobjective(p)
PBS = readtable("VDNE control.csv");
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
V8 = RSV.V(RSV.Day==8); V8=V8(~isnan(V8)); V8 = V8(V8~=0);
V15 = RSV.V(RSV.Day==15); V15=V15(~isnan(V15)); V15 = V15(V15~=0);
V21 = RSV.V(RSV.Day==21); V21=V21(~isnan(V21)); 
V20 = RSV.V(RSV.Day==20); V20=V20(~isnan(V20)); V20 = V20(V20~=0);
V25 = RSV.V(RSV.Day==25); V25=V25(~isnan(V25)); V25 = V25(V25~=0);
a = 30;
b = 100;
yRSV = [median(log(V8)) median(D8) median(E8);
        median(log(V15)) median(D15) median(E15);
        log(1E-14) median(D21) median(E21);
        median(log(V20)) median(D20) median(E20);
        median(log(V25)) median(D25) median(E25)]; %ydata [3; 10; 16; 5; 10]
yRSV(:,2) = yRSV(:,2)/a;
yRSV(:,3) = yRSV(:,3)/b;
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

yPBS = [median(D8) median(E8);
      median(D15) median(E15);
      median(D21) median(E21);
      median(D20) median(E20);
      median(D25) median(E25)]; %ydata [3; 10; 16; 5; 10]

yPBS(:,1) = yPBS(:,1)/a;
yPBS(:,2) = yPBS(:,2)/b;
%considering weight for each error to be inverse of the variance in
%experimental values
wRSV = [1/(var(yRSV(:,1))) 1/(var(yRSV(:,2))) 1/(var(yRSV(:,3)))];
wPBS = [1/(var(yPBS(:,1))) 1/(var(yPBS(:,2)))];

td = [3; 10; 16; 0; 5; 10];

%ys has ODE solutions for both infections
ysRSV = solveODE(p,td,[[log(p(12)) yRSV(1,2) yRSV(1,3)]; ...
    [log(p(12)*2.3 + exp(yRSV(2,1))) yRSV(2,2) yRSV(2,3)]]);
ysPBS = solveODE(p,td,[[0 yPBS(1,1) yPBS(1,2)]; ...
    [0 yPBS(2,1) yPBS(2,2)]]);

SSE = sum(sum([wRSV(1)*((yRSV(:,1)-ysRSV(:,1)).^2), wRSV(2)*((yRSV(:,2)-ysRSV(:,2)).^2), ...
    wRSV(3)*((yRSV(:,3)-ysRSV(:,3)).^2), wPBS(1)*((yPBS(:,1)-ysPBS(:,2)).^2), wPBS(2)*((yPBS(:,2)-ysPBS(:,3)).^2)]));
end

function C = solveODE(p,td,y0)
    flag=0;
    [~,Cv1]=ode23s(@DifEq,td(1:3),y0(1,:)); %first infection
    [~,Cv2]=ode23s(@DifEq,td(4:6),y0(2,:)); %second infection
    function dY = DifEq(td,y)
        dydt = zeros(3,1);
        
        dydt(1) = p(1)*(1+p(2)*y(2)) - p(3); %V(t)
        
        dydt(2) =  p(4)*exp(y(1))*(p(5)-y(2)) + p(6)*y(3)*(p(5)-y(2)) - p(7)*y(2); %D(t)
        
        if (y(2) < p(8)) && flag==0
            dydt(3) = p(9) - p(10)*y(3);
        elseif (y(2) >= p(8)) || flag==1
            %setting flag=1 is changing the switch
            flag=1;
            dydt(3) = p(9) + p(11)*y(2) - p(10)*y(3); %E(t)
        end
        
        dY = dydt;
    end
    C = [Cv1(1:3,:); Cv2(2:3,:)]; %merging two solutions in one matrix for optimisation
end

function [X_normalized, X_scaled] = lhsdesign_modifie(n,min_ranges_p,max_ranges_p)
p=length(min_ranges_p);
[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(p,n);
OFFSET=ones(p,n);

for i=1:p
    SLOPE(i,:)=ones(1,n).*slope(i);
    OFFSET(i,:)=ones(1,n).*offset(i);
end
X_normalized = lhsdesign(p,n);

X_scaled=SLOPE.*X_normalized+OFFSET;
end