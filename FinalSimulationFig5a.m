%simulations with optimised parameters
%p = [p(1) p(2) p(3) p(4) p(5) p(6) p(7) p(8) p(9)];

yBarrierDay3 = [1/0.839648935 1/0.353323191 1/0.115518237 1/0.188301372 1/0.176568128 1/0.186769508 1/0.337745948 1/0.03208497 1/0.080471812 1/0.054285354];
yBarrierDay10 = [1/1.276306677 1/0.610758183 1/0.615412799 1/0.374977825 1/1.267745426 1/0.660307652 1/0.17837704 1/0.841418693 1/0.636951896 1/0.483848203];
yBarrierDay16 = [1/1.090913165 1/1.081785198 1/0.999439071 1/1.044157971 1/1.219417999 1/1.165326546 1/0.8858519 1/1.160685532 1/1.809803744 1/1.13153952];
yBarrierDay20 = [1/0.115797611 1/0.141492786 1/0.202708497 1/0.17244626 1/0.08480893 1/0.112396627 1/0.223579564];
yBarrierDay25 = [1/0.841057034 1/1.209916175 1/0.311331235 1/0.815515205 1/0.98539008 1/1.048210749 1/1.292862041];
 
yEosDay3 = [4.792613351 1.096345106 0.965338606 0.500882474 0.786028585 0.206569863 1.204833601 13.89785519 13.96226479 0.87457368];
yEosDay10 = [4.946314523 1.682851648 8.744956467 1.596681907 1.15400261 3.105117682 2.196844295 6.276676902 11.08503063 0.544576469];
yEosDay16 = [0 0 0 0 0.190774433 0.214160749 0 0 0.40491282 0.256533682];
yEosDay20 = [113.9537344 193.5868035 182.514377 141.0068497 156.5307501 32.15323472 21.96016202 38.27545754];
yEosDay25 = [48.16071071 280.6401386 72.49751595 101.3586246 17.80923603 28.0879391 23.59473214];

yVirusDay3 = log([5.91063E-07 3.74451E-07 2.372E-07 1.5356E-07 6.83663E-07 8.18975E-07 7.22344E-07 9.28624E-07 2.74216E-07]);
yVirusDay10 = log([1.24E-09 0.000000011 1.06E-08 3.21E-08 9.7E-09 1.45E-08 1E-11 1E-11 1E-11 1E-11]);
yVirusDay16 = log([1E-11 1E-11 1E-11 1E-11 1E-11 1E-11 1E-11 1E-11 1E-11 1E-11]);
yVirusDay20 = log([1.64471E-08 2.69464E-07 5.33486E-08 1.38817E-07 2.45703E-07 3.04392E-09 1.0785E-06 1.92277E-07]);
yVirusDay25 = log([5.1661E-08 0.000000001 7.228E-09 4.922E-09 8.4843E-08 5.573E-09 0.000000001 1.234E-08]);


yd = [log(5.91063E-07) 5.508864668 1.030841856;
      log(5.47E-09) 1.597451426 2.650980988;
      log(1E-11) 0.900207483 0;
      log(1.65547E-07) 7.067498103 127.480292;
      log(6.4005E-09) 1.014826535 48.16071071]; %ydata [3; 10; 16; 5; 10]

param = [0.15 0.0004 0.87 2.04E5 1.19e-06 0.25 21.9 0.89 4.56e-06 6.56];
% param = [2.61E-5 0.49 1.49 403.58 8.6E-7 0.2 3.47 0.28 0.02 9.15];

td = [8 15 21 20 25];
tspan = linspace(0,16);
[yFitI, a] = solveFinalSol(param, tspan, [0 3 10 16], [log(param(9)) 0 1]);

[yFitII, b] = solveFinalSol(param, tspan, [0 5 10], [log(param(9)+exp(a(3,1))) a(3,2) a(3,3)]);
% yFitII(:,1) = exp(yFitII(:,1));

ymax = [min([min(yFitII(:,1)) min(yFitI(:,1))]) 10 max([max(yFitII(:,3)) 150])]; 
Names = {'log(V(t))', 'D(t)', 'E(t)'};

%loop to generate simulations for first infection
for i=1:3
    plot(linspace(5,21),yFitI(:,i),'b-','MarkerSize',20,'LineWidth',4)
    hold on;
    plot(linspace(15,30),yFitII(:,i),'g-','MarkerSize',20,'LineWidth',4)
    hold on
    plot(td(1:3),yd(1:3,i),'bo','MarkerSize',20,'LineWidth',4)
    hold on
    plot(td(4:5),yd(4:5,i),'go','MarkerSize',20,'LineWidth',4)
    if(i==1)
        ylim([ymax(i) 0]);
    elseif(i==2)
        hold on
        plot(5,0,'bo','MarkerSize',20,'LineWidth',4)
        hold on
        yline(param(10))
        ylim([0 ymax(i)]);
    elseif(i==3)
        hold on
        plot(5,1,'bo','MarkerSize',20,'LineWidth',4)
        ylim([0 ymax(i)]);
    end
    xlabel('Age of mice (days)');
    ylabel(Names{i});
    ax = gca;
    ax.FontSize = 15;
    xlim([0 30]);
    xticks(0:5:30);
    xtickangle(0);
%     ylim([0 ymax(i)]);
    hold off;
    legend('$1^{st} Infection$', '$2^{nd} Infection$', 'Interpreter', 'latex');
    saveas(gcf, strcat('~/Library/Mobile Documents/com~apple~CloudDocs/Research/PhD/Asthma/RSV/Model/',Names{i}, '.png'));
end

function [y,ySpec] = solveFinalSol(p, tspan, td, y0)
    flag=0;
    [~,y] = ode23s(@DifEq,tspan,y0);
    [~,ySpec] = ode23s(@DifEq,td,y0);
    function dY = DifEq(tspan,y)
        dydt = zeros(3,1);

        dydt(1) = p(1)*(1+p(2)*y(2)) - p(3);

        dydt(2) =  p(4)*exp(y(1))*(13.48-y(2)) + p(5)*y(3)*(13.48-y(2)) - p(6)*y(2); %D(t)
        
        if (y(2) < p(10)) && flag==0
            dydt(3) = 0;
        elseif (y(2) >= p(10)) || flag==1
            flag=1;
            dydt(3) = p(7)*y(2) - p(8)*y(3); %E(t)
        end
        dY = dydt;
    end
end