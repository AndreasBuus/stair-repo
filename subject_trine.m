clc; 
clear all; 
close all; 

SubjectName = "Trine"; 

% Folders path 
folderpath_10sem = "C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/data/Trine 21.04.2022/";

% File Path 
filepath_CTL = folderpath_10sem + "control1.mat";
filepath_CTL_2 = folderpath_10sem + "control2.mat"; 
filepath_HOR = folderpath_10sem + "horisontal1.mat"; 
filepath_VER = folderpath_10sem + "vertical1.mat"; 

%% Abbreviation

% protocol abbreviation types
CTL = 1; VER = 2; HOR = 3; CTL2 = 4; 
ProtoAll = [CTL, VER, HOR, CTL2];

% sensor abbreviation type
SOL = 1; TA = 2; ANG = 3; FSR = 4; time = 5; 

% sensor abbreviation subtypes
VEL = 1; ACC = 2; 

data    = cell(3,4);
% example: data{protocol, sensor}(sweep, data number)


%% Load data and Acquisition Set-Up from Mr Kick
% Functions used: [load_EMG()]

data    = cell(numel(ProtoAll),4); % data{protocol, EMG}(sweep, data_num)

% Load control 
[data{CTL,1:4}] = load_EMG(filepath_CTL); clear filepath_CTL
[data{CTL2,1:4}] = load_EMG(filepath_CTL_2); clear filepath_CTL_2
[data{VER,1:4}] = load_EMG(filepath_VER); clear filepath_VER
[data{HOR,1:4}] = load_EMG(filepath_HOR); clear filepath_HOR

% Acquisition Set-Up
sweep_length = 10;              % Signal length in second
Fs = 2000;                      % Samples per second
dt = 1/Fs;                      % Seconds per sample
pre_trig = 4;                   % Pre-trigger 
N = Fs*sweep_length;            % Total number of samples per signal

% Exclude data 
exclude_CTL = [];               % excluded control sweeps
exclude_CTL2= [];               % excluded control sweeps
exclude_VER = [10,20,26];       % excluded horizontal sweeps
exclude_HOR = [20,35];          % excluded horizontal sweeps


%% Perturbation and exclude trails
% Functions used: [sort_sweeps()]

HOR_perturbation = [1,8,11,14,16,20,22,25,29,31,33,38,39,42,43,45,52,54,55,57,60]; % Checked

VER_perturbation = [1,8,12,14,17,19,22,29,31,33,34,37,38,41,43,45,53,56,58,59];    % Checked 
VER_perturbation_delay    = VER_perturbation([3,4,8,9,10,12,13,14,15,17,18,20]);
VER_perturbation_notDelay = VER_perturbation([2,5,6,7,11,16,19]);


[VER_yes, VER_no] = sort_sweeps(size(data{VER,SOL},1), VER_perturbation,  exclude_VER); 
[VER_delay]    = sort_sweeps(size(data{VER,SOL},1), VER_perturbation_delay   ,  exclude_VER); 
[VER_notDelay] = sort_sweeps(size(data{VER,SOL},1), VER_perturbation_notDelay,  exclude_VER); 
[HOR_yes, HOR_no] = sort_sweeps(size(data{HOR,SOL},1), HOR_perturbation,  exclude_HOR); 

for i = [SOL, TA, FSR, ANG]
    data{CTL,i}(exclude_CTL,:) = []; 
    data{CTL2,i}(exclude_CTL2,:) = []; 

    data{VER,i}(exclude_VER,:) = []; 
    data{HOR,i}(exclude_HOR,:) = []; 
end 

%% Load offset 

raw = data;                     % save raw data 
% offset = zeros(1,size(data{HOR,SOL},1)); 

load('C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/offset/offset_trine_ver'); 
offset_ver = offset; 
load('C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/offset/offset_trine_hor'); 
offset_hor = offset; 


%% Filtrering and detrend (similar to MR. kick)

fc = 40;                            % Cutoff frequency for LowPass filter
order = 1;                          % Filter order 
[b,a] = butter(order,fc/(Fs/2));    % Filter coefficient

for proto = ProtoAll
    [data{proto,SOL}, data{proto,TA}] = rectify_filter(data{proto,SOL}, data{proto,TA}, b, a);  % rectify and filter EMG
    [data{proto,FSR}] = func_filt_FSR(data{proto,FSR}, "test", 0 , "limit_pct", 0.95); 
end


%% Begining of Stand- and Swingphase
% find change in FSR signal
step_index = cell(3,1);
error_index = cell(3,1);

for proto = ProtoAll
    [step_index{proto}, error_index{proto}] = func_step_index(data{proto,FSR});
end

sweep = 1; [step_index{CTL}(sweep,:)] = func_step_index_corr('FSR', data{CTL,FSR}, 'direction', 'right', 'edge_num',8, 'move_num', 2, 'sweep',sweep);
sweep = 4; [step_index{CTL}(sweep,:)] = func_step_index_corr('FSR', data{CTL,FSR}, 'direction', 'right', 'edge_num',8, 'move_num', 2, 'sweep',sweep);

sweep = 11; [step_index{HOR}(sweep,:)] = func_step_index_corr('FSR', data{HOR,FSR}, 'direction', 'right', 'edge_num',8, 'move_num', 2, 'sweep',sweep);
sweep = 14; [step_index{HOR}(sweep,:)] = func_step_index_corr('FSR', data{HOR,FSR}, 'direction', 'right', 'edge_num',8, 'move_num', 2, 'sweep',sweep);
sweep = 31; [step_index{HOR}(sweep,:)] = func_step_index_corr('FSR', data{HOR,FSR}, 'direction', 'right', 'edge_num',8, 'move_num', 2, 'sweep',sweep);


if true 
    % >>>> TEST CODE <<<<
    proto = VER; % Checked: CTL[x], HOR[x], VER[x], CTL2[x]
    loop = true; sweep = 1; 
    prompt = "Continue, press >c<" + newline + "Quite, press >q<"+ newline + "Change sweep number, press >t<"+ newline;
    figure(2); 

    while loop == true
        clc
        sgtitle("Sweep: " + sweep + ". Protocol: " + proto)
        hold off;
        plot(data{proto, ANG}(sweep,:))
        hold on;
        plot(data{proto, FSR}(sweep,:))

        plot([8000,8000],[-1 6], "LineWidth",3, "Color", "red")

        [rise, fall] = func_find_edge(0); 
        plot([step_index{proto}(sweep, rise), step_index{proto}(sweep,rise)],[1 4], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,rise), step_index{proto}(sweep,fall)],[2.5 2.5], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,fall), step_index{proto}(sweep,fall)],[1 4], "LineWidth",2, "Color", "red")

        [rise, fall] = func_find_edge(2); 
        plot([step_index{proto}(sweep,rise), step_index{proto}(sweep,rise)],[1 4], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,rise), step_index{proto}(sweep,fall)],[2.5 2.5], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,fall), step_index{proto}(sweep,fall)],[1 4], "LineWidth",2, "Color", "red")

        [rise, fall] = func_find_edge(4); 
        plot([step_index{proto}(sweep, rise), step_index{proto}(sweep,rise)],[1 4], "LineWidth",2, "Color", "green")
        plot([step_index{proto}(sweep,rise), step_index{proto}(sweep,fall)],[2.5 2.5], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,fall), step_index{proto}(sweep,fall)],[1 4], "LineWidth",2, "Color", "blue")
    
        [rise, fall] = func_find_edge(6); 
        plot([step_index{proto}(sweep, rise), step_index{proto}(sweep,rise)],[1 4], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,rise), step_index{proto}(sweep,fall)],[2.5 2.5], "LineWidth",2, "Color", "black")
        plot([step_index{proto}(sweep,fall), step_index{proto}(sweep,fall)],[1 4], "LineWidth",2, "Color", "red")

        [rise, fall] = func_find_edge(7); 
        plot([step_index{proto}(sweep, rise), step_index{proto}(sweep,rise)],[1 4], "LineWidth",2, "Color", "yellow")

        correctInput = false; 
        while correctInput == false
            str = input(prompt, 's');
            if strcmp(str,"q")
                disp("Loop stopped")
                loop = false; correctInput = true; 
            elseif strcmp(str,"t")
                sweep = input("New sweep number: ")-1; 
                correctInput = true; 
            elseif strcmp(str,"c") %, sweep == size(data{proto, SOL}))
                correctInput = true; 
            end 
            if correctInput == false
                warning("Input not accepted")
            end
        end
        sweep = sweep + 1;
        if sweep > size(data{proto,SOL},1)
            loop = false; 
        end 
    end
    close 2
end

%% Save Data 

type{1} = VER_yes;   type{2} = VER_no; 
type{3} = HOR_yes;   type{4} = HOR_no; 
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName + "_data", 'data')
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName+"_type",'type')
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName+"_step",'step_index')

disp("Processed: " + SubjectName)