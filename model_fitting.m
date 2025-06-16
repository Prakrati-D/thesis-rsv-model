%% Script with all functions for running optimisation problem

% initializations
rng default
% 
% % number of parameter samples to test
n=100;
% 
%lower and upper bounds for model parameter search
load('bounds','lb','ub');

% generates n samples of initial guesses for model parameters
[x,initial_guess] =  lhsdesign_modifie(n,lb,ub);
% 
% % initialise store size of the optimisted parameter set
output_p = [];
% 
% for parallel computing
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end
 
% for parallel computing and storing the output of each iteration
parfor i = 1:n
    % testing one parameter set at a time from n samples
    p0 = initial_guess(:,i);

    % calls the globalsearch algorithm
    param = runGS(p0);

    % stores optimised parameter set for each iteration and corresponding
    % cost function value
    if(size(param)~=0)
        [output] = param;
        cost = objective(param);
        [output] = [output; cost];
        output_p = [output_p output];
    end
end

% collects smallest cost function value and corresponding optimised
% parameter set
writematrix(output_p,'ModelOutput.csv');
cost_min = min(output_p(end,:));
cost_min_col = find(output_p(end,:)==cost_min);
optimised_p = output_p(1:end -1, cost_min_col).'; 

%% Global search function
% runs globalsearch for each initial parameter guess

% % t = timeit(@runGS);
% w = warning('query','last');
% id = w.identifier;
% warning('off',id)

function r = runGS(p0)
    %lower and upper bounds for model parameter search
    load('bounds','lb','ub');
    tic
    % set up optimisation problem
%     problem = createOptimProblem('fmincon','objective',@objective,'x0',p0,'lb',lb,'ub',ub);
%     gs = GlobalSearch('PlotFcn',[]);
%     [r, fval, exitflag, output] = run(gs,problem);
%     fval
%     exitflag
%     output
%     [r, fval, exitflag, output] = ga(@objective,15,[],[],[],[],lb,ub);
    [r, fval, exitflag, output] = patternsearch(@objective,p0,[],[],[],[],lb,ub);
%     opts = optimoptions("fmincon",'MaxFunctionEvaluations',10000);
%     [r, fval, exitflag, output] = surrogateopt(@objective,lb,ub,[],[],[],[],[]);
    fval
    exitflag
    output
    toc
end

%% Define object function
% Sum of square errors
 
function SSE = objective(p)

    % load data used for objective function
    load('data_RSV.mat','yRSV','wRSV');
    load('data_PBS.mat','yPBS','wPBS');

    td = [0; 3; 10; 16; 0; 5; 10];
    
    % ys has ODE solutions for both infections
    
    ysRSV = solveODERSV(p,td,[[log(p(12)) p(13) p(14)]; ...
    [log(p(15) + exp(yRSV(2,1))) yRSV(2,2) yRSV(2,3)]]);
    ysPBS = solveODEPBS(p,td,[[0 p(13) p(14)]; ...
    [0 yPBS(2,1) yPBS(2,2)]]);
    maxD1 = solveODE(p,linspace(0,16),[[log(p(12)) p(13) p(14)]; ...
    [log(p(15) + exp(yRSV(2,1))) yRSV(2,2) yRSV(2,3)]], 1/5);
    maxD2 = solveODE(p,linspace(0,16),[[log(p(12)) p(13) p(14)]; ...
    [log(p(15) + exp(yRSV(2,1))) yRSV(2,2) yRSV(2,3)]],5);
    SSE1 = sum(sum([wRSV(1)*((yRSV(:,1)-ysRSV(:,1)).^2), wRSV(2)*((yRSV(:,2)-ysRSV(:,2)).^2), ...
    wRSV(3)*((yRSV(:,3)-ysRSV(:,3)).^2), wPBS(1)*((yPBS(:,1)-ysPBS(:,2)).^2), wPBS(2)*((yPBS(:,2)-ysPBS(:,3)).^2)]));SSE1
    SSE2 = ((maxD1-p(8))^2)*wRSV(2);SSE2
    SSE3 = ((p(8)-maxD2)^2)*wRSV(2);SSE3
    SSE = sum([SSE1,SSE2,SSE3]);
%     SSE = sum(sum([wRSV(1)*((yRSV(:,1)-ysRSV(:,1)).^2), wRSV(2)*((yRSV(:,2)-ysRSV(:,2)).^2), ...
%     wRSV(3)*((yRSV(:,3)-ysRSV(:,3)).^2), wPBS(1)*((yPBS(:,1)-ysPBS(:,2)).^2), wPBS(2)*((yPBS(:,2)-ysPBS(:,3)).^2)]));
end

%% ODE solver function
% function that solves ODEs using built in solvers ode15s/ode23s

function sols = solveODERSV(p,td,y0)
    
    % switch initially off
    flag = 0;
    
    % outputs from ode solver
    [~,Cv1] = ode23s(@DifEqs,td(1:4),y0(1,:)); %first infection
    [~,Cv2] = ode23s(@DifEqs,td(5:7),y0(2,:)); %second infection
    
    function dY = DifEqs(t,y)
        
        % initialise the size of the system
        dydt = zeros(3,1);
        
        % variables
        V = y(1); % V(t) virus
        D = y(2); % D(t) damage
        E = y(3); % E(t) eosinophils

        % ODE equations
        dydt(1) = p(1)*(1+p(2)*D) - p(3); %V(t)
        
        dydt(2) =  p(4)*exp(V)*(p(5)-D) + p(6)*E*(p(5)-D) - (p(7))*D; %D(t)
        
        if (y(2) < p(8)) && flag==0
            dydt(3) = p(9) - p(10)*E;
        elseif (y(2) >= p(8)) || flag==1
            %setting flag=1 is changing the switch
            flag=1;
            dydt(3) = p(9) + p(11)*D - p(10)*E; %E(t)
        end
        
        dY = dydt;
    end

    % merge solutions in one matrix for optimisation
    sols = [Cv1(2:4,:); Cv2(2:3,:)];

end

function sols = solveODEPBS(p,td,y0)
    
    % switch initially off
    flag = 0;
    
    % outputs from ode solver
    [~,Cv1] = ode23s(@DifEqs,td(1:4),y0(1,:)); %first infection
    [~,Cv2] = ode23s(@DifEqs,td(5:7),y0(2,:)); %second infection
    
    function dY = DifEqs(t,y)
        
        % initialise the size of the system
        dydt = zeros(3,1);
        
        % variables
        V = y(1); % V(t) virus
        D = y(2); % D(t) damage
        E = y(3); % E(t) eosinophils

        % ODE equations
        dydt(1) = p(1)*(1+p(2)*D)*V - p(3)*V; %V(t)
        
        dydt(2) =  p(4)*V*(p(5)-D) + p(6)*E*(p(5)-D) - (p(7))*D; %D(t)
        
        if (y(2) < p(8)) && flag==0
            dydt(3) = p(9) - p(10)*E;
        elseif (y(2) >= p(8)) || flag==1
            %setting flag=1 is changing the switch
            flag=1;
            dydt(3) = p(9) + p(11)*D - p(10)*E; %E(t)
        end
        
        dY = dydt;
    end

    % merge solutions in one matrix for optimisation
    sols = [Cv1(2:4,:); Cv2(2:3,:)];

end

function sols = solveODE(p,td,y0,k)
    
    % switch initially off
    flag = 0;
    
    % outputs from ode solver
    [~,Cv1] = ode23s(@DifEqs,td,y0(1,:)); %first infection
    [~,Cv2] = ode23s(@DifEqs,td,y0(2,:)); %second infection
    
    function dY = DifEqs(t,y)
        
        % initialise the size of the system
        dydt = zeros(3,1);
        
        % variables
        V = y(1); % V(t) virus
        D = y(2); % D(t) damage
        E = y(3); % E(t) eosinophils

        % ODE equations
        dydt(1) = p(1)*(1+p(2)*D) - p(3); %V(t)
        
        dydt(2) =  p(4)*exp(V)*(p(5)-D) + p(6)*E*(p(5)-D) - (p(7)*k)*D; %D(t)
        
        if (y(2) < p(8)) && flag==0
            dydt(3) = p(9) - p(10)*E;
        elseif (y(2) >= p(8)) || flag==1
            %setting flag=1 is changing the switch
            flag=1;
            dydt(3) = p(9) + p(11)*D - p(10)*E; %E(t)
        end
        
        dY = dydt;
    end

    % merge solutions in one matrix for optimisation
    if(k<1)
      maxD = max(Cv1(:,2));
    elseif(k>1)
      maxD = max(Cv2(:,2));
    end
    sols = maxD;
end

%% Parameter sampling function
% function generates n samples of parameters from the specified range

function [X_normalized, X_scaled] = lhsdesign_modifie(n,min_ranges_p,max_ranges_p)
p = length(min_ranges_p);

[M,N] = size(min_ranges_p);
    if M < N
        min_ranges_p = min_ranges_p';
    end
    
[M,N] = size(max_ranges_p);
    if M < N
        max_ranges_p = max_ranges_p';
    end

slope = max_ranges_p-min_ranges_p;
offset = min_ranges_p;

SLOPE = ones(p,n);
OFFSET = ones(p,n);

    for i = 1:p
        SLOPE(i,:) = ones(1,n).*slope(i);
        OFFSET(i,:) = ones(1,n).*offset(i);
    end

X_normalized = lhsdesign(p,n);
X_scaled = SLOPE.*X_normalized+OFFSET;

end