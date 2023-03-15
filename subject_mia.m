clc; 
clear all; 
close all; 

SubjectName = "Mia"; 

% Folders path 
folderpath_10sem = "C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/data/Mia 20.4.2022/";

% File Path 
filepath_CTL = folderpath_10sem + "mia_control1_(20_4_2022).mat";
filepath_CTL_2 = folderpath_10sem + "mia_control2_(20_4_2022).mat"; 
filepath_HOR = folderpath_10sem + "mia_horisontal1_(20_4_2022).mat"; 
filepath_VER = folderpath_10sem + "mia_vertical1_(20_4_2022).mat"; 

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

% include function folder
addpath("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/FunctionFiles")

%% Load data and Acquisition Set-Up from Mr Kick
% Functions used: [load_EMG()]


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
exclude_CTL = [15,16];               % excluded control sweeps
exclude_CTL2= [];               % excluded control sweeps
exclude_VER = [];               % excluded horizontal sweeps
exclude_HOR = [];     % excluded horizontal sweeps


%HOR_perturbation = [7,8,16,18,30,33,35,38,39,41,49,52,54,56,59,60];
HOR_perturbation = [7,8,11,13,16,18,21,25,27,28,30,33,35,38,39,41,49,52,54,56,59,60];
VER_perturbation = [5,8,10,11,15,17,20,24,28,30,32,35,37,40,42,44,52,54,56,58];


[VER_yes, VER_no] = sort_sweeps(size(data{VER,SOL},1), VER_perturbation,  exclude_VER); 
[HOR_yes, HOR_no] = sort_sweeps(size(data{HOR,SOL},1), HOR_perturbation,  exclude_HOR); 


for i = [SOL, TA, FSR, ANG]
    data{CTL,i}(exclude_CTL,:) = []; 
    data{CTL2,i}(exclude_CTL2,:) = []; 

    data{VER,i}(exclude_VER,:) = []; 
    data{HOR,i}(exclude_HOR,:) = []; 
end 

%% Load offset 

raw = data;                     % save raw data 
load('C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/offset/offset_mia_ver'); 
offset_ver = offset; 
load('C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/offset/offset_mia_hor'); 
offset_hor = offset; 

%% Filtrering and detrend (similar to MR. kick)
% Functions used: [rectify_filter()], [filt_FSR()]

fc = 40;                            % Cutoff frequency for LowPass filter
order = 1;                          % Filter order 
[b,a] = butter(order,fc/(Fs/2));    % Filter coefficient

% rectify and filter EMG. Remove noise in FSR 
for proto = ProtoAll
    [data{proto,SOL}, data{proto,TA}] = rectify_filter(data{proto,SOL}, data{proto,TA}, b, a);  
end

%gap_size2
%gap_size3

[data{CTL,FSR}] = func_filt_FSR(data{proto,FSR}, "test", 0 , "limit_pct", 0.95, "gap_size2", 20, 'gap_size3',20); 
[data{VER,FSR}] = func_filt_FSR(data{proto,FSR}, "test", 0 , "limit_pct", 0.95); 
[data{HOR,FSR}] = func_filt_FSR(data{proto,FSR}, "test", 0 , "limit_pct", 0.95); 
[data{CTL2,FSR}] = func_filt_FSR(data{proto,FSR}, "test", 0 , "limit_pct", 0.95); 



%% 𝔹𝕖𝕘𝕚𝕟𝕚𝕟𝕘 𝕠𝕗 𝕊𝕥𝕒𝕟𝕕- 𝕒𝕟𝕕 𝕊𝕨𝕚𝕟𝕘𝕡𝕙𝕒𝕤𝕖
% Functions used: [step_index_pos()], [step_index_corr()], [step_corr()]
% find change in FSR signal
step_index = cell(3,1);
error_index = cell(3,1);

for proto = ProtoAll
    [step_index{proto}, error_index{proto}] = func_step_index(data{proto,FSR});
end

if true 
    % >>>> TEST CODE <<<<
    proto = CTL; % CTL[], HOR[x], VER[], CTL2
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

% 
% step_index    = cell(3,1);
% for PROTO = ProtoAll
%     [step_index{PROTO}] = step_index_pos(data{PROTO,FSR});     % find change in FSR signal
% end
% 
% sweep = 12; step_corr{HOR}(sweep,[5,4]) = [11863, 14120];
% sweep = 13; step_corr{HOR}(sweep,[5,4]) = [11863, 14120];
% sweep = 15; step_corr{HOR}(sweep,[5,4]) = [12035, 13890];
% sweep = 16; step_corr{HOR}(sweep,[5,4]) = [12035, 13890];
% sweep = 23; step_corr{HOR}(sweep,[5,4]) = [11939, 13594];
% 
% 
% for sweep = [2] %Tilbage
%     temp = step_index{HOR}(sweep,:);
%     step_corr{HOR}(sweep,4)= temp(6); 
%     step_corr{HOR}(sweep,5)= temp(7);         
% end
% 
% for sweep = [16,19,20,21,24,25,27,32,34,35,38,39,40,41,45,51,54,55,56,60] %Frem
%     temp = step_index{HOR}(sweep,:);
%     step_corr{HOR}(sweep,4)= temp(2); 
%     step_corr{HOR}(sweep,5)= temp(3);         
% end
% 
% 
% step_corr{VER}(5,5)=11774; step_corr{VER}(5,4)=13323; 
% step_corr{VER}(6,5)=12166; step_corr{VER}(6,4)=13916; 
% step_corr{VER}(7,5)=11953; step_corr{VER}(7,4)=13403; 
% step_corr{VER}(17,5)=12129; step_corr{VER}(17,4)=13759; 
% step_corr{VER}(57,5)=12393; step_corr{VER}(57,4)=14406; 
% 
% 
% [step_index{HOR}] = step_index_corr(data{HOR,FSR}, step_index{HOR}, step_corr{HOR}); 
% [step_index{VER}] = step_index_corr(data{VER,FSR}, step_index{VER}, step_corr{VER}); 
% 
% sweep = 7; step_index{HOR}(sweep,[7,6,5,4]) = [9870,11348,11990,15291]; 
% sweep = 18; step_index{HOR}(sweep,[7,6,5,4]) = [9987,11453,12152,15052]; 
% sweep = 30; step_index{HOR}(sweep,[7,6,5,4]) = [9870,11460,12167,14731]; 
% sweep = 36; step_index{HOR}(sweep,[7,6,5,4]) = [9996,11348,12184,13864]; 
% sweep = 45; step_index{HOR}(sweep,[4]) = [13798]; 
% sweep = 46; step_index{HOR}(sweep,[7,6,5,4]) = [9870,11357,11749,13844]; 
% sweep = 47; step_index{HOR}(sweep,[7,6,5,4]) = [9913,11417,12166,13723]; 
% sweep = 52; step_index{HOR}(sweep,[7,6,5,4]) = [9790,11379,12142,13923]; 
% sweep = 57; step_index{HOR}(sweep,[7,6,5,4]) = [9659,11263,12119,13741]; 
% 
% 
% step_index{VER}(9,7)=9987; step_index{VER}(9,6)=11312; step_index{VER}(9,5)=12108; step_index{VER}(9,4)=13555; 
% step_index{VER}(10,7)=9997; step_index{VER}(10,6)=11381; step_index{VER}(10,5)=12030; step_index{VER}(10,4)=13758; 
% step_index{VER}(11,7)=9803; step_index{VER}(11,6)=11332; step_index{VER}(11,5)=11655; step_index{VER}(11,4)=13704; 
% step_index{VER}(15,7)=9926; step_index{VER}(15,6)=11318; step_index{VER}(15,5)=12130; step_index{VER}(15,4)=13658; 
% sweep = 16; step_index{HOR}(sweep,[7,6,5,4]) = [9917,11170,11973,13608]; 
% 
% step_index{VER}(18,7)=9881; step_index{VER}(18,6)=11215; step_index{VER}(18,5)=11999; step_index{VER}(18,4)=13585; 
% step_index{VER}(19,7)=9689; step_index{VER}(19,6)=11063; step_index{VER}(19,5)=11685; step_index{VER}(19,4)=13356; 
% step_index{VER}(20,7)=9859; step_index{VER}(20,6)=11159; step_index{VER}(20,5)=11691; step_index{VER}(20,4)=13500; 
% step_index{VER}(22,7)=9728; step_index{VER}(22,6)=11127; step_index{VER}(22,5)=11502; step_index{VER}(22,4)=13337; 
% step_index{VER}(24,7)=9467; step_index{VER}(24,6)=11132; step_index{VER}(24,5)=11870; step_index{VER}(24,4)=13660; 
% step_index{VER}(25,7)=9826; step_index{VER}(25,6)=11144; step_index{VER}(25,5)=11829; step_index{VER}(25,4)=13461; 
% step_index{VER}(28,7)=9858; step_index{VER}(28,6)=11433; step_index{VER}(28,5)=11935; step_index{VER}(28,4)=13574; 
% step_index{VER}(29,7)=9913; step_index{VER}(29,6)=11472; step_index{VER}(29,5)=12297; step_index{VER}(29,4)=13805; 
% step_index{VER}(30,7)=9691; step_index{VER}(30,6)=11284; step_index{VER}(30,5)=11821; step_index{VER}(30,4)=13685; 
% step_index{VER}(31,4)=13580; 
% step_index{VER}(32,7)=10181; step_index{VER}(32,6)=11633; step_index{VER}(32,5)=12696; step_index{VER}(32,4)=14132; 
% step_index{VER}(33,7)=9878; step_index{VER}(33,6)=11342; step_index{VER}(33,5)=11665; step_index{VER}(33,4)=13588; 
% step_index{VER}(35,7)=9824; step_index{VER}(35,6)=11378; step_index{VER}(35,5)=12077; step_index{VER}(35,4)=13807; 
% step_index{VER}(36,7)=10029; step_index{VER}(36,6)=11529; step_index{VER}(36,5)=11935; step_index{VER}(36,4)=13972; 
% step_index{VER}(38,7)=10006; step_index{VER}(38,6)=11733; step_index{VER}(38,5)=12127; step_index{VER}(38,4)=14122; 
% step_index{VER}(39,7)=9934; step_index{VER}(39,6)=11485; step_index{VER}(39,5)=12203; step_index{VER}(39,4)=13865; 
% step_index{VER}(40,7)=9992; step_index{VER}(40,6)=11487; step_index{VER}(40,5)=12234; step_index{VER}(40,4)=13877; 
% step_index{VER}(41,7)=9680; step_index{VER}(41,6)=11690; step_index{VER}(41,5)=12231; step_index{VER}(41,4)=14168; 
% step_index{VER}(43,7)=9519; step_index{VER}(43,6)=11328; step_index{VER}(43,5)=12179; step_index{VER}(43,4)=13725; 
% step_index{VER}(44,7)=10001; step_index{VER}(44,6)=11358; step_index{VER}(44,5)=12254; step_index{VER}(44,4)=13705; 
% step_index{VER}(45,7)=9890; step_index{VER}(45,6)=11354; step_index{VER}(45,5)=11951; step_index{VER}(45,4)=13623; 
% step_index{VER}(47,7)=9803; step_index{VER}(47,6)=11480; step_index{VER}(47,5)=12220; step_index{VER}(47,4)=13865; 
% step_index{VER}(48,7)=9522; step_index{VER}(48,6)=11489; step_index{VER}(48,5)=11850; step_index{VER}(48,4)=13931; 
% step_index{VER}(49,7)=9839; step_index{VER}(49,6)=11378; step_index{VER}(49,5)=12076; step_index{VER}(49,4)=13724; 
% step_index{VER}(50,7)=9481; step_index{VER}(50,6)=11285; step_index{VER}(50,5)=11909; step_index{VER}(50,4)=13717; 
% step_index{VER}(52,7)=9495; step_index{VER}(52,6)=11195; step_index{VER}(52,5)=11558; step_index{VER}(52,4)=13500; 
% step_index{VER}(53,7)=9697; step_index{VER}(53,6)=11179; step_index{VER}(53,5)=12068; step_index{VER}(53,4)=13419; 
% step_index{VER}(54,7)=9565; step_index{VER}(54,6)=11449; step_index{VER}(54,5)=12373; step_index{VER}(54,4)=13852; 
% step_index{VER}(56,7)=9619; step_index{VER}(56,6)=11488; step_index{VER}(56,5)=11902; step_index{VER}(56,4)=14018; 
% 
% 
% %CTL -tilbage
% for sweep = [2,5]
%     temp = step_index{CTL}(sweep,:); 
%     step_corr{CTL}(sweep,4)= temp(6); 
%     step_corr{CTL}(sweep,5)= temp(7);         
% end
% 
% %CTL - frem
% for sweep = [11,18]
%     temp = step_index{CTL}(sweep,:); 
%     step_corr{CTL}(sweep,4)= temp(2); 
%     step_corr{CTL}(sweep,5)= temp(3);         
% end
% 
% %CTL2 -tilbage
% for sweep = []
%     temp = step_index{CTL2}(sweep,:); 
%     step_corr{CTL2}(sweep,4)= temp(6); 
%     step_corr{CTL2}(sweep,5)= temp(7);         
% end
% 
% %CTL2 - frem
% for sweep = [1]
%     temp = step_index{CTL2}(sweep,:); 
%     step_corr{CTL2}(sweep,4)= temp(2); 
%     step_corr{CTL2}(sweep,5)= temp(3);         
% end
% 
% 
% [step_index{CTL}] = step_index_corr(data{CTL,FSR}, step_index{CTL}, step_corr{CTL}); 
% [step_index{CTL2}] = step_index_corr(data{CTL2,FSR}, step_index{CTL2}, step_corr{CTL2}); 
% 
% 
% % -- var: "Stand_duration" 
% % Duration of stand phase of tread 4
% % stand_duration{PROTO}(sweep)
% stand_DUR    = cell(3,1);
% for PROTO = ProtoAll
%     for i = 1:size(data{PROTO,FSR},1)
%         stand_DUR{PROTO}(i) = (step_index{PROTO}(i,4) - step_index{PROTO}(i,5))*dt; 
%     end
% end
% 
% 
% if false 
%     % >>>> TEST CODE <<<<
%     proto = CTL2; % checked VER, HOR, CTL 
%     loop = true; sweep = 1; 
%     prompt = "Continue, press >c<" + newline + "Quite, press >q<"+ newline + "Change sweep number, press >t<"+ newline;
%     
%     while loop == true
%         clc
%         figure; hold on; sgtitle("sweep: " + sweep)
%         plot(data{proto, ANG}(sweep,:))
%         plot(data{proto, FSR}(sweep,:))
%         plot([8000, 8000],[-2 6], "LineWidth",3, "Color", "black")
%         plot([step_index{proto}(sweep,4), step_index{proto}(sweep,4)],[0 5], "LineWidth",2, "Color", "green")
%         plot([step_index{proto}(sweep,5), step_index{proto}(sweep,5)],[0 5], "LineWidth",2, "Color", "blue")
%         plot([step_index{proto}(sweep,6), step_index{proto}(sweep,6)],[0 5], "LineWidth",2, 'color', 'black')
%         plot([step_index{proto}(sweep,7), step_index{proto}(sweep,7)],[0 5], "LineWidth",2, "Color", 'black')
% 
%         
%         correctInput = false; 
%         while correctInput == false
%             str = input(prompt, 's');
%             if strcmp(str,"q")
%                 disp("Loop stopped")
%                 loop = false; correctInput = true; 
%             elseif strcmp(str,"t")
%                 sweep = input("New sweep number: ")-1; 
%                 correctInput = true; 
%             elseif strcmp(str,"c") %, sweep == size(data{proto, SOL}))
%                 correctInput = true; 
%             end 
%     
%             if correctInput == false
%                 warning("Input not accepted")
%             end
%         end
%         close all; 
%         sweep = sweep + 1;
%         if sweep > size(data{proto,SOL},1)
%             loop = false; 
%         end 
%     end
% end

%% Save Data 

type{1} = VER_yes;   type{2} = VER_no; 
type{3} = HOR_yes;   type{4} = HOR_no; 
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName + "_data", 'data')
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName+"_type",'type')
save("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/" + SubjectName+"_step",'step_index')

disp("Processed: " + SubjectName)