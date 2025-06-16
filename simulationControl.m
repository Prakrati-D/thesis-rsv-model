read_data

load('data_PBS.mat','yPBS','wPBS','yPBSmd');
    
param = [0.8829108210727442 0.0022197538338062594 1.5855337824641154 96479.30905559342 1204.4009093540556 3.066471440462994e-09 0.5932602880626062 100.00000000000004 150.77389797772875 1.6815434782131855 46.91994852372584 8.197388672722167e-07 72.70233764974128 6.386317175899307e-09 3.5762895623721156e-06];

td = [8 15 21 20 25];
tspan = linspace(0,16);
[yFitI, a] = solveFinalSol(param, tspan, [0 3 10 16], [0 0 param(14)]);
yFitI(:,2)

[yFitII, b] = solveFinalSol(param, tspan, [0 5 10], [0 a(3,2) a(3,3)]);
yFitII(:,2)

ymax = [1 400 max([max(yFitII(:,3)) 9000])]; 
Names = {'Virus V(t)', 'Barrier Damage D(t)', 'Eosinophils E(t)'};

%loop to generate simulations for first infection
for i=1:3
    if(i==1)
        plot(1,1,'b-','MarkerSize',20,'LineWidth',4)
        hold on
        ax = gca;
        ax.YLim = [-11 -5];
        ax.YTick = [-11 -10 -9 -8 -7 -6 -5];
        ax.YTickLabel = {'10^{-11}','10^{-10}', '10^{-9}', '10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}'};
        set(gca, 'YMinorTick','on');
        set(gca, 'XTickLabel',[]);
    else
    plot(linspace(5,21),yFitI(:,i),'b--','MarkerSize',20,'LineWidth',4)
    hold on;
    plot(linspace(15,30),yFitII(:,i),'g--','MarkerSize',20,'LineWidth',4)
    hold on
    errorbar(td(1:3),yPBS(1:3,i-1),yPBSmd(1:3,i-1),'LineWidth',3,'Color','b','LineStyle','none')
    hold on
    plot(td(1:3),yPBS(1:3,i-1),'bo','MarkerSize',12,'LineWidth',3)
    hold on
    errorbar(td(4:5),yPBS(4:5,i-1),yPBSmd(4:5,i-1),'LineWidth',3, 'Color','g','LineStyle','none')
    hold on
    plot(td(4:5),yPBS(4:5,i-1),'go','MarkerSize',12,'LineWidth',3)
    hold on
    if(i==2)
        yline(param(8),'Color','r','LineWidth',2)
        ylim([0 400]);
        set(gca, 'XTickLabel',[]);
    elseif(i==3)
        ylim([0 10000])
        yticks([0 2000 4000 6000 8000 10000]);
        yticklabels({'0', '2x10^{3}', '4x10^{3}', '6x10^{3}', '8x10^{3}', '1x10^{4}'})
%         xlabel('Age of mice (days)');
        xlim([0 30]);
    
    xtickangle(0);
    end
    end
    hold on
    xline(5,'LineStyle','--','LineWidth',2)
    hold on
    xline(15,'LineStyle','--','LineWidth',2)
    xlim([0 30]);
    xticks(0:5:30);
    ax = gca;
    ax.FontSize = 25;
    ax.LineWidth = 1.5;
    hold off;
    grid on
    saveas(gcf, strcat('~/Library/Mobile Documents/com~apple~CloudDocs/Research/PhD/Asthma/RSV/Model/',Names{i}, '.png'));
end

function [y,ySpec] = solveFinalSol(p, tspan, td, y0)
    flag=0;
    [~,y] = ode15s(@DifEq,tspan,y0);
    [~,ySpec] = ode15s(@DifEq,td,y0);
    function dY = DifEq(t,y)

        dydt = zeros(3,1);
        
        % variables
        V = 0; % V(t) virus
        D = y(2); % V(t) damage
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
end