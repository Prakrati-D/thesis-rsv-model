%% Script to read in data
% read data from csv files and organise into samples

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

yRSV = [median(log10(V8)) median(D8) median(E8);
        median(log10(V15)) median(D15) median(E15);
        log10(3E-11) median(D21) median(E21);
        median(log10(V20)) median(D20) median(E20);
        median(log10(V25)) median(D25) median(E25)]; %ydata [3; 10; 16; 5; 10]

yRSVmd = [mad(log10(V8),1) mad(D8,1) mad(E8,1);
          mad(log10(V15),1) mad(D15,1) mad(E15,1);
          mad(log10(V21),1) mad(D21,1) mad(E21,1);
          mad(log10(V20),1) mad(D20,1) mad(E20,1);
          mad(log10(V25),1) mad(D25,1) mad(E25,1)];

%%

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

yPBSmd = [mad(D8,1) mad(E8,1);
          mad(D15,1) mad(E15,1);
          mad(D21,1) mad(E21,1);
          mad(D20,1) mad(E20,1);
          mad(D25,1) mad(E25,1)];


%considering weight for each error to be inverse of the variance in
%experimental values
wRSV = [1/(var(yRSV(:,1))) 1/(var(yRSV(:,2))) 1/(var(yRSV(:,3)))];
wPBS = [1/(var(yPBS(:,1))) 1/(var(yPBS(:,2)))];


% save the variable yRSV as yRSV.mat
save('data_RSV', 'yRSV', 'wRSV', 'yRSVmd')
save('data_PBS', 'yPBS', 'wPBS', 'yPBSmd')