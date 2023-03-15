clc;
clear all; 
close all; 

%% Abbreviation
CTL = 1; VER = 2; HOR = 3; CTL2 = 4; time = 5;
SOL = 1; TA = 2; ANG = 3; FSR = 4; FOO = 5; 
VEL = 1; ACC = 2; 
ProtoAll = [CTL, VER, HOR]; 

%% Files 
filepath = "C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/data/test/test.mat"; 

%% Load Data and Acquisition Set_up from Mr. Kick 
[data{VER,1:4}] = load_EMG(filepath); 

sweep_length = 10;              % Signal length in second
Fs = 2000;                      % Samples per second
dt = 1/Fs;                      % Seconds per sample
pre_trig = 4;                   % Pre-trigger 
N = Fs*sweep_length;            % Total number of sample

%% Functions
msToSec = @(x) x*10^-3; % ms to sec 
secToMs = @(x) x*10^3; % ms to sec 

%% Filter
[data{VER,FSR}] = rescale(filt_FSR( data{VER,FSR} ));  


%% Rescale Ankle values
a = 3.5480; 
b = 20.6570; 

for i = 1:size(data{VER,ANG},1) %sweeps
    data{VER,ANG}(i,:) = data{VER,ANG}(i,:)*a+b; % [ankle angle]
end
 

%% Align
included = [2,9,12,15,17,24,27,30,32,36,38];
before = -2000;
after  = 2000; 
first_Pk = [12108,12005,11998,11829,12071,11776,11888,11912,11957,11941,11867];
second_Pk = [45, 80, 79, 49, 50, 72, 48, 85, 64, 72, 63]; 
total_pk = first_Pk% + second_Pk; 

for i = 1:numel(included)
    idx_arr = before+total_pk(i):after+total_pk(i);
    align{VER,ANG}(i,:) = data{VER,ANG}(included(i),idx_arr); 
    align{VER,FSR}(i,:) = data{VER,FSR}(included(i),idx_arr); 
end 
align{VER,time} = linspace(secToMs(before*dt), secToMs(after*dt-dt), numel(idx_arr) ); 


%% Find colapse time
for i = 1:numel(included) % loop through sweeps
    FSR_data = align{VER,FSR}(i,:);
    edge_indexes(i) = find(edge(FSR_data)); 
end 


%% Speed and aceleration
span = 10;  % sets the span of the moving average to span.

for i = 1:size(align{VER,ANG},1) %sweeps
    align{VER,ANG}(i,:) = smooth(align{VER,ANG}(i,:), span, 'moving'); % [ankle angle]
    align{VER,VEL}(i,:) = diff(align{VER,ANG}(i,:))/(dt*10^3);     % [deg/sec]
    align{VER,ACC}(i,:) = diff(smooth(align{VER,VEL}(i,:), span, 'moving'))/(dt*10^3);     % [deg/sec]
end
avg_edge = secToMs(mean(edge_indexes+before)*dt); 



%% Plot
sweep = 7;  
patchcolor = "black"; 
FaceAlpha = 0.1; 

figure; sgtitle("")

subplot(311); hold on; ylim([-20 30])
    yyaxis left; ylabel("Position" + newline + "[Deg]")
        x = [-200 200 200 -200];
        y = [-20 -20 30 30];
        patch(x,y,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
    %plot(align{VER,time}, align{VER,ANG}(sweep,:) )
    plot(align{VER,time}, mean(align{VER,ANG},1) , "color", "black")

    yyaxis right; ylabel("Fourth step" + newline + "push down")
    plot(align{VER,time}, align{VER,FSR}(sweep,:), "LineWidth",1.4, "Color","red" )

    ax = gca;
    ax.YAxis(1).Color = 'black';
    ax.YAxis(2).Color = 'red';

subplot(312); hold on;  xlim([-220 220]); ylim([-12 0])
    yyaxis left; ylabel("Position" + newline + "[Deg]"); 
    x = [-200 200 200 -200];
    y = [-20 -20 30 30];
    patch(x,y,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
    plot(align{VER,time}, mean(align{VER,ANG},1) , "color", "black")
    
    yyaxis right; ylabel("Fourth step" + newline + "push down") 
    plot([avg_edge avg_edge],[0 1], "LineWidth",1.4, "Color","red")
    ax = gca;
    ax.YAxis(1).Color = 'black';
    ax.YAxis(2).Color = 'red';


subplot(313); hold on; xlim([-220 220])
    yyaxis left;  ylabel("Velocity" + newline + "[Deg/ms]");  ylim([-0.6 0.3])
    x = [-200 200 200 -200];
    y = [-20 -20 30 30];
    patch(x,y,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
    plot(align{VER,time}(1:end-1), mean(align{VER,VEL},1), "color", "black")

    yyaxis right; ylabel("Fourth step" + newline + "push down") 
    plot([avg_edge avg_edge],[0 1], "LineWidth",1.4, "Color","red")
    ax = gca;
    ax.YAxis(1).Color = 'black';
    ax.YAxis(2).Color = 'red';

    xlabel("Time after estimated perturbation [ms] ")


%offset = zeros(1,size(data{VER,SOL},1)); 
%save('C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/offset/offset_test_ver','offset')