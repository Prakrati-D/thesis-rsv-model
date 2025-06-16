%% Script to read in data
% read data from csv files and organise into samples

opts = detectImportOptions('./data_rsv_iwoffii.csv');

opts = setvartype(opts, {'V', 'D', 'N'}, {'double', 'double', 'double'});

A = readtable("./data_rsv_iwoffii.csv", opts);

%%
%retrieving all data

V.RSV_day20 = A.V(A.Group=="RSV" & A.Day==20); V.RSV_day20 = V.RSV_day20(~isnan(V.RSV_day20));

D.PBS_day20 = A.D(A.Group=="PBS" & A.Day==20); D.PBS_day20 = D.PBS_day20(~isnan(D.PBS_day20));
D.RSV_day20 = A.D(A.Group=="RSV" & A.Day==20); D.RSV_day20 = D.RSV_day20(~isnan(D.RSV_day20)); 

E.PBS_day20 = A.E(A.Group=="PBS" & A.Day==20); E.PBS_day20 = E.PBS_day20(~isnan(E.PBS_day20));
E.RSV_day20 = A.E(A.Group=="RSV" & A.Day==20); E.RSV_day20 = E.RSV_day20(~isnan(E.RSV_day20));

N.PBS_day6 = A.N(A.Group=="PBS" & A.Day==6); N.PBS_day6 = N.PBS_day6(~isnan(N.PBS_day6));
N.RSV_day6 = A.N(A.Group=="RSV" & A.Day==6); N.RSV_day6 = N.RSV_day6(~isnan(N.RSV_day6));

N.PBS_day16 = A.N(A.Group=="PBS" & A.Day==16); N.PBS_day16 = N.PBS_day16(~isnan(N.PBS_day16));
N.RSV_day16 = A.N(A.Group=="RSV" & A.Day==16); N.RSV_day16 = N.RSV_day16(~isnan(N.RSV_day16));

N.PBS_day20 = A.N(A.Group=="PBS" & A.Day==20); N.PBS_day20 = N.PBS_day20(~isnan(N.PBS_day20));
N.RSV_day20 = A.N(A.Group=="RSV" & A.Day==20); N.RSV_day20 = N.RSV_day20(~isnan(N.RSV_day20));

%% finding scales for each variable

V.scale_day20 = max([max(V.RSV_day20)]);
D.scale_day20 = max([max(D.PBS_day20), max(D.RSV_day20)]);
E.scale_day20 = max([max(E.PBS_day20), max(E.RSV_day20)]);
N.scale_day6 = max([max(N.PBS_day6), max(N.RSV_day6)]);
N.scale_day16 = max([max(N.PBS_day16), max(N.RSV_day16)]);
N.scale_day20 = max([max(N.PBS_day20), max(N.RSV_day20)]);

%% scaled medians

V.RSV_day20_median = median((V.RSV_day20));

D.PBS_day20_median = median(D.PBS_day20);
D.RSV_day20_median = median(D.RSV_day20); 

E.PBS_day20_median = median(E.PBS_day20);
E.RSV_day20_median = median(E.RSV_day20);

% N.PBS_day6_median = median(N.PBS_day6./N.scale_day6);
% N.RSV_day6_median = median(N.RSV_day6./N.scale_day6);
% 
% N.PBS_day16_median = median(N.PBS_day16./N.scale_day16);
% N.RSV_day16_median = median(N.RSV_day16./N.scale_day16);
% 
% N.PBS_day20_median = median(N.PBS_day20./N.scale_day20);
% N.RSV_day20_median = median(N.RSV_day20./N.scale_day20);

%% scaled SEM of the medians

V.RSV_SEM = mad(V.RSV_day20,1);

D.PBS_SEM = mad(D.PBS_day20,1);
D.RSV_SEM = mad(D.RSV_day20,1);

E.PBS_SEM = mad(E.PBS_day20,1);
E.RSV_SEM = mad(E.RSV_day20,1);

% N.PBS_day6_SEM = mad(N.PBS_day6./N.scale_day6,1);
% N.RSV_day6_SEM = mad(N.RSV_day6./N.scale_day6,1);
% 
% N.PBS_day16_SEM = mad(N.PBS_day16./N.scale_day16,1);
% N.RSV_day16_SEM = mad(N.RSV_day16./N.scale_day16,1);
% 
% N.PBS_day20_SEM = mad(N.PBS_day20./N.scale_day20,1);
% N.RSV_day20_SEM = mad(N.RSV_day20./N.scale_day20,1);

%% fold-change medians

% V.fold_day20 = V.AI_RSV_day20_median / V.RSV_day20_median;

D.fold_day20 = D.RSV_day20_median / D.PBS_day20_median;

E.fold_day20 = E.RSV_day20_median / E.PBS_day20_median;

% N.fold_day6 = N.RSV_day6_median / N.PBS_day6_median;
% 
% N.fold_day16 = N.RSV_day16_median / N.PBS_day16_median;
% 
% N.fold_day20 = N.RSV_day20_median / N.PBS_day20_median;

%% fold-change SEM

% V.fold_SEM = error_propagation(V.AI_RSV_day20_median, V.AI_RSV_SEM, V.RSV_day20_median, V.RSV_SEM);

D.fold_SEM = error_propagation(D.RSV_day20_median, D.RSV_SEM, D.PBS_day20_median, D.PBS_SEM);

E.fold_SEM = error_propagation(E.RSV_day20_median, E.RSV_SEM, E.PBS_day20_median, E.PBS_SEM);

% N.fold_SEM_day6 = error_propagation(N.RSV_day6_median, N.RSV_day6_SEM, N.PBS_day6_median, N.PBS_day6_SEM);
% 
% N.fold_SEM_day16 = error_propagation(N.RSV_day16_median, N.RSV_day16_SEM, N.PBS_day16_median, N.PBS_day16_SEM);
% 
% N.fold_SEM_day20 = error_propagation(N.RSV_day20_median, N.RSV_day20_SEM, N.PBS_day20_median, N.PBS_day20_SEM);
% 
%%

V.time = 20;
D.time = 20;
E.time = 20;
N.time = [6, 16, 20];

function SE = error_propagation(mean_1, error_1, mean_2, error_2)
    fc = mean_1 / mean_2;
    relative_error_1 = error_1 / mean_1 ;
    relative_error_2 = error_2 / mean_2 ;
    SE = fc * sqrt(relative_error_1^2 + relative_error_2^2);
end