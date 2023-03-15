clc 
clear; 
close all; 

%% Folders 
fprintf('script: Folder . . . '); tic

names = ["thomas", "benedikte", "andreas", "andrew", "gritt", "maria", "trine", "trine2"]%, "Christian", "Soeren"]; % MIA
% CTL2: andrew, Gritt, thomas, Mia, Trine, Trine 2

% Define data path
addpath("C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/FunctionFiles")
folderpath_preprocessed_data = "C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/data_preprocessed/"; 

% preallocation
total_data = cell(1,1,numel(names));
total_type = cell(1,1,numel(names)); 
total_step = cell(1,1,numel(names)); 

% load data 
for i = 1:length(names)
    load(folderpath_preprocessed_data + names(i) + "_data.mat");   total_data{:,:,i} = data; % size(data) = [3,4]
    % example: data{protocol, sensor}(sweep, data number)

    load(folderpath_preprocessed_data + names(i) + "_type.mat");   total_type{:,:,i} = type; % size(type) = [1,4]
    % example: type{1:2} = [VER_yes, VER_no]; type{3:4} = [HOR_yes, HOR_no]; 

    load(folderpath_preprocessed_data + names(i) + "_step.mat");   total_step{:,:,i} = step_index; % size(step_index) = [3,1]
    % example: step_index{protocol}(sweep, step)
end 

fprintf('done [ %4.2f sec ] \n', toc);

%% Abbreviation
fprintf('script: Abbreviation . . . '); tic

% protocol abbreviation types
CTL = 1; VER = 2; HOR = 3; CTL2 = 4; 
proto_all = [CTL, VER, HOR];

% Sensor abbreviation type
SOL = 1; TA = 2; ANG = 3; FSR = 4; time = 5; VEL = 6; ACC = 7;

% Plotting labels 
labels = ["Soleus"; "Tibialis"; "Position"; "";  ""; "Velocity"; "Acceleration"];
labels_ms = ["Soleus"+newline+"[\muV]";"Tibialis"+newline+"[\muV]"; "Position"+newline+"[Deg]";  "";  "Time"+newline+"[ms]"; "Velocity"+newline+"[Deg/ms]";"Acceleration"+newline+"[Deg/ms^2]"];
labels_sec = ["Soleus"+newline+"[\muV]";"Tibialis"+newline+"[\muV]"; "Position"+newline+"[Deg]"; "";  "Time"+newline+"[sec]";"Velocity"+newline+"[Deg/s]";"Acceleration"+newline+"[Deg/s^2]"];

% Global arrays
align_with_obtions = ["second_begin", "four_begin", "six_begin"];
steps_tested = [2,4,6];

% Global function
msToSec = @(x) x*10^-3;         % Ms to sec 
secToMs = @(x) x*10^3;          % Sec to ms 

fprintf('done [ %4.2f sec ] \n', toc);

%% Acquisition Set-Up
fprintf('script: Acquisition Set-Up . . . '); tic; 

sweep_length = 10;              % Signal length in second
Fs = 2000;                      % Samples per second
dt = 1/Fs;                      % Seconds per sample
pre_trig = 4;                   % Pre-trigger 
N = Fs*sweep_length;            % Total number of samples per signal

fprintf('done [ %4.2f sec ] \n', toc);

%% Readjust data to Local Peak instead of FSR 
fprintf('script: Readjust data to Local Peak instead of FSR . . .'); tic
readjust = true; 
show_gui = false; 
gui_subject = 4; 
protocol = CTL; 

if readjust 
    for subject = 1:numel(names) % loop through subjects
       
        % Load data for subject 
        data = total_data{1,1,subject}; 
        step_index = total_step{1,1,subject};
        
        for sweep = 1:size(data{protocol,ANG},1) % loop through sweeps        
            for step = 1:3  % loop through steps   
                
                % Find template around foot-strike
                [rise_num, ~] = func_find_edge(steps_tested(step));   
                rise_index = step_index{protocol}(sweep, rise_num); 
                array = rise_index-400:rise_index+400;
                template = data{protocol,ANG}(sweep, array); 
                signal = data{protocol,ANG}(sweep, :); 
            
                % Peak inside template
                [pks, locs] = findpeaks(template, 'MinPeakDistance', 200);
                locs = locs + array(1);
                
                % Find the peak that follow the condition 
                the_pks = 0; the_loc = 0;  % peaks and locations
                for i = 1:numel(locs) % loop through locations
                    if pks(i) > signal(locs(i) - 200) && pks(i) > signal(locs(i)+200)
                        the_pks = pks(i);
                        the_loc = locs(i);
                    end 
                end 

                % Update step_index or throw error
                if the_pks == 0     % non found: error
                    msg = "\n Error: no peak found. Subject: " + subject + ". Sweep: " + sweep + ". Step: " + steps_tested(step) + " "; 
                    fprintf(2,msg); 
                else                % no error 
                    step_index{protocol}(sweep, rise_num) = the_loc;
                end 
            end 
        end
        total_step{1,1,subject} = step_index; % update step_index
    end 
    fprintf('done [ %4.2f sec ] \n', toc);
else 
    fprintf('disable \n');
end 


if show_gui
    fprintf('script: Re-adjust gui - [Waiting for user input]')
    data = total_data{1,1,gui_subject}; 
    step_index = total_step{1,1,gui_subject};
    offset = []; 
    readjustFSR
    pause
%     if ~empty(offset)
%         total_step{1,1,gui_subject} = offset; 
%     end 
   fprintf('gui done \n')
end 

%% remove saturated data 
fprintf('script: Remove saturated data . . .'); tic

remove_saturated = true;    % enable or disable 
threshold = -10;            % remove ANG data if lower than 
span = [0, 20];             % ms 

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
if remove_saturated  

    % From ms to samples 
    span = msToSec(span)*Fs;        % from ms to sample
    exc_ctl = cell(size(names));    % exclude Control
    
    for sub = 1:numel(names) % loop through subjects
        % Load data 
        data = total_data{1,1,sub}; 
        step_index = total_step{1,1,sub};
       
        % Find saturated idxs
        for sweep = 1:size(data{CTL,FSR},1) % loop through sweeps 
            for step = 1:3 % loop through steps 
                [rise] = func_find_edge(steps_tested(step)); 
                edge = step_index{CTL}(sweep,rise);
                y = data{CTL,ANG}(sweep, span(1)+edge:span(2)+edge); 
                if any(find(y<threshold))
                    exc_ctl{sub} = unique([exc_ctl{sub}, sweep]);  
                end
            end 
        end
    
        % Remove saturated idxs
        step_index{CTL}(exc_ctl{sub},:) = []; 
        for i = [SOL, TA, FSR, ANG]
            data{CTL,i}(exc_ctl{sub},:) = []; 
        end 
    
        % Save data
        total_data{1,1,sub} = data; 
        total_step{1,1,sub} = step_index;
    end
    
    fprintf('done [ %4.2f sec ] \n', toc);
else 
    fprintf('disable \n');
end 

%% Normalize EMG (make as a function instead) 
fprintf('script: Normalize EMG  . . . '); tic; 

normalize = true;      % enable or disable
span = 20;             % how big is the smooth span.
normalizing_step = 0;  % which step is the data being normalized to [0,2,4]?

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
if normalize 
    for sub = 1:length(names) % loop through subjects  
        % Load data 
        data = total_data{1,1,sub};
        step_index = total_step{1,1,sub};

        if ~any(strcmp(names(sub), ["Christian", "Soeren"])) % christian and soeren only completed the prebaseline protocol 
            [data] = func_normalize_EMG(step_index, data, 'protocols', [CTL,VER,HOR],  'normalize_to_step', normalizing_step,'span', span);
        else
            [data] = func_normalize_EMG(step_index, data, 'protocols', CTL,  'normalize_to_step', normalizing_step,'span', span);
        end
        total_data{1,1,sub} = data; % update data
    end
    fprintf('done [ %4.2f sec ] \n', toc);
else 
    fprintf('disable \n');
end 

%% Speed and aceleration (make as a function instead) 
fprintf('script: Speed and aceleration . . . '); tic; 

plot_data = false;          % enable or disable plot
span_position = 5;          % inc. sample in guassian filter span 
span_velocity = 5;          % inc. sample in guassian filter span 
span_acceleration = 10;     % inc. sample in gaussian filter span 

for sub = 1:length(names) % subjects
    data = total_data{1,1,sub};  % load data 
    for proto = proto_all % protocols 
    
        % Position 
        data{proto,ANG} = data{proto,ANG}.*4+25; % rescale the signal; 
        pos = data{proto,ANG};
        %data{proto,ANG} = smoothdata(data{proto,ANG}, 2, 'gaussian', span_position);    % gaussian smoothing

        % velocity
        diffs1 = diff(data{proto,ANG}, 1, 2)./(dt*10^3);            % [deg/sample]
        diffs1 = padarray(diffs1, [0 1], 'post');                   % zeropadding           
        data{proto,VEL} = smoothdata(diffs1, 2, 'gaussian', span_velocity);    % gaussian smoothing

        % acceleration
        diffs2 = diff(data{proto,VEL}, 1, 2)./(dt*10^3);            % [deg/sample^2]
        diffs2 = padarray(diffs2, [0 1], 'post');                   % zeropadding
        data{proto,ACC} = smoothdata(diffs2, 2, 'gaussian', span_acceleration);    % gaussian smoothing
       
        % Need to plot before and after 
        if plot_data == 1 && sub == 1
            step_index = total_step{1,1,sub};
            sweep = 1; dur = 1000; before = 500;
            [rise, ~] = func_find_edge(4);
            rise_index = step_index{CTL}(sweep,rise);
            display_array = rise_index-before:rise_index+dur;

            figure; hold on
            subplot(311)
            plot(display_array, pos(1,display_array),'-o', display_array, data{proto,ANG}(sweep,display_array),'-x')
            legend(["raw", "filtered"])

            subplot(312)
            plot(display_array, diffs1(1,display_array),'-o',display_array, data{proto,VEL}(sweep,display_array),'-x')
            legend(["raw", "filtered"])

            subplot(313)
            plot(display_array, diffs2(1,display_array),'-o', display_array, data{proto,ACC}(sweep,display_array),'-x')
            legend(["raw", "filtered"])

        end
    end 
    total_data{1,1,sub} = data; 
end
fprintf('done [ %4.2f sec ] \n', toc);

%% TASK0.1: Simplify the matlab script
% Make the script accessible to be navigated by Andrew and Thomas. 
% - clearify parameters meant to be adjusted. 
% - make passive code as function and remove.

%% TASK0.2 Show average sweep for single subject
fprintf('script: TASK0.2 Show average sweep  . . . ');

show_plot = true;      % Disable or enable plot
subject = 2;            % Obtions: 1:8
proto = CTL;            % Obtions: CTL, VER, HOR 
str_sen = ["Position", "Soleus", "Tibialis"];    % Obtions: "Soleus", "Tibialis","Position", "Velocity", "Acceleration"; 
show_FSR = true; 

align_bool = false;      % Should the data be aligned with step specified in "Align with specific Stair step"
    alignWithStep = 'second_begin'; 
    before = 200; 
    after = 100; 

% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
if show_plot
fprintf(' plot subject: >> ' + names(subject) + ' << \n' );

if align_bool
    % Load data and align with defined step 
    data = total_data{1,1,subject}; 
    step_index = total_step{1,1,subject};
    temp = cell(3,7); 
    [temp{proto,:}] = func_align(step_index{proto}, data{proto,[1:4,6:7]}, 'sec_before', msToSec(before), 'sec_after', msToSec(after), 'alignStep', alignWithStep);
    clear data
    data = temp;
    x_axis = data{proto,time};
    x_axis = secToMs(x_axis);
    str_xlabel = "Time [ms]" + newline + "Data normalized to step four"; 
else 
    % Load data 
    data = total_data{1,1,subject}; 
    x_axis = linspace(-4, 6-dt, N); 
    x_axis = secToMs(x_axis); 
    str_xlabel = "Time [ms]" + newline + "Data normalized to Force-Platform";
end 

switch proto
    case HOR
        type = total_type{:,:,subject};
        yes = type{3}; no = type{4}; str_title = "Horizontal perturbation"; 
    case VER
        type = total_type{:,:,subject};
        yes = type{1}; no = type{2}; str_title = "Vertical perturbation"; 
    case {CTL}
        str_title = "Pre-baseline Control";
end

figure; 
    sgtitle(str_title + " - subject " + names(subject)); 
    for i = 1:size(str_sen,2) % check and plot data included in str_sen
        switch str_sen(i) 
            case "Soleus"
                sensor_modality = SOL; 
                str_ylabel = "Soleus" + newline + "[\muV]"; 
            case "Tibialis"
                sensor_modality = TA; 
                str_ylabel = "Tibialis"+newline+"[\muV]"; 
            case "Position"
                sensor_modality = ANG;
                str_ylabel = "Position" + newline + "[Deg]" + newline + "[Dorsal] <  > [Plantar]"; 
            case "Velocity"
                sensor_modality = VEL; 
                str_ylabel = "Velocity" + newline + "[Deg/s]"; 
            case "Acceleration"
                sensor_modality = ACC; 
                str_ylabel = "Acceleration " + newline + "[Deg/s^2]"; 
             otherwise
                disp("ERROR" + newline + "String: >>" + str_sen(i) + "<< is not registered.")
        end
        
        switch proto
            case {VER, HOR} 
                subplot(size(str_sen,2)*100 + 10 + i); hold on; ylabel(str_ylabel); 
                plot(x_axis, mean(data{proto,sensor_modality}(no,:),1), "LineWidth",3, "color",[0.75, 0.75, 0.75])
                plot(x_axis, mean(data{proto,sensor_modality}(yes,:),1), "LineWidth",1, "color","black")

                if show_FSR
                    y_fsr = rescale(mean(data{proto,FSR},1)); 
                    yyaxis right; ylabel("Phase"); ylim([-0.1 1.1])
                    plot(x_axis, y_fsr, "color",	"yellow");
                end
                
            case {CTL, CTL2}
                subplot(size(str_sen,2)*100 + 10 + i); hold on; ylabel(str_ylabel); 
                % plt std around mean 
                y = mean(data{proto,sensor_modality},1); 
                std_dev = std(data{proto,sensor_modality});
                curve1 = y + std_dev;
                curve2 = y - std_dev;
                x2 = [x_axis, fliplr(x_axis)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, [0.75, 0.75, 0.75], 'LineStyle','none'); 
                plot(x_axis, y, 'LineWidth', 2, "color","black"); 

                if show_FSR % show FSR if enabled 
                    y_fsr = rescale(mean(data{proto,FSR},1)); 
                    yyaxis right; ylabel("Phase"); ylim([-0.1 1.1])
                    plot(x_axis, y_fsr, "color",	"yellow");
                end
        end
        if i == size(str_sen,2)
            xlabel(str_xlabel)
        end 
    end
else 
    fprintf('disable \n');
end 


%% TASK0.3 Show individual sweep for single subject
% undersøg for metodisk fejl.
fprintf('script: TASK0.3 Show individual sweep  . . . ');

show_plot = false;       % Disable or enable plot
subject = 2;            % Obtions: 1:9
proto = CTL;            % Obtions: CTL
sensor_modality = ANG;  % Obtions: SOL, TA, ANG, VEL, ACC
before = 50;            % before foot strike included in ms
after = 0;              % after foot lift-off included in ms

% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
if show_plot
fprintf('plot subject: >> ' + names(subject) + ' << \n' );
data = total_data{1,1,subject}; 
step_index = total_step{1,1,subject};
x_axis_total = linspace(-4, 6-dt, N); % time axis 
sweepNum = size(data{proto,ANG},1);   % total sweep size

figure; % Begin plot
    loop = true; sweep = 1; 
    while loop == true
        clc % clear cmd promt 
        sgtitle("Sweep: " + sweep) % display current sweep in promt

        subplot(211);  
        % plot data 
        yyaxis left;
        plot(x_axis_total, mean(data{proto,sensor_modality},1), '-', "LineWidth",3, "color",[0.75, 0.75, 0.75]) % mean plt 
        hold on 
        plot(x_axis_total, data{proto,sensor_modality}(sweep,:), '-', "LineWidth",1, "color","black") % sweep plt 
        hold off
      
        % plt formalia 
        xlabel("Time"+newline+"[sec]")
        ylabel(labels_ms(sensor_modality))
        title(['Black graph, Sweep data. {\color{gray} Gray graph, Mean data [n=' num2str(sweepNum) '].}'])

        if show_FSR
            y_fsr = rescale(data{proto,FSR}(sweep,:)); 
            yyaxis right; ylabel("Phase"); ylim([-0.1 1.1]); 
            plot(x_axis_total,y_fsr, "color",	"red"); 
        end
    
        % Define the data for each step 
        clear y
        for k = 1:3
            clear data_align
            data_align = cell(3,7); 
            [data_align{proto,:}] = func_align(step_index{proto}, data{proto,[1:4,6:7]}, 'sec_before', msToSec(before), 'sec_after', msToSec(after), 'alignStep', align_with_obtions(k));
            x_axis = data_align{proto,time};
            y{k,2} = secToMs(x_axis);
            y{k,1} = data_align{proto,sensor_modality}(sweep,:); 
            y{k,3} = mean(data_align{proto,sensor_modality},1);
        end 
    
        % plot data for each step 
        for k = 1:3 % loop through steps
            subplot(233+k); hold on; 
            if sensor_modality == or(SOL,TA)
                ylim([0 ceil(max([max(y{1,1}),  max(y{2,1}), max(y{3,1})])/100)*100])
            end
            plot(y{k,2}, y{k,3}, "LineWidth",2, "color", [0.75, 0.75, 0.75]) % Mean
            plot(y{k,2}, y{k,1}, "LineWidth",1, "color","black") % Sweep
            ylim auto
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            plot([0, 0],[YL(1) YL(2)], "--","LineWidth",1, "Color", "red") % plt fsr 
    
            % plt formalia 
            xlabel("Time"+newline+"[ms]")
            title("Step: " + steps_tested(k))
            ylabel(labels_ms(sensor_modality))
        end

        % Wait for user input
        correctInput = false; 
        prompt = "Continue, press >c<" + newline + "Quite, press >q<"+ newline + "Change sweep number, press >t<"+ newline;    
        while correctInput == false     % Wait for correct user input
            str = input(prompt, 's');   % Save user input
            if strcmp(str,"q")          % If pressed - Quite the loop 
                disp("Loop stopped")
                loop = false; correctInput = true; 
            elseif strcmp(str,"t")      % If pressed - Change sweep number
                sweep = input("New sweep number: ")-1; 
                correctInput = true; 
            elseif strcmp(str,"c")      % If pressed - Continue to next sweep
                correctInput = true; 
            end 
    
            if correctInput == false
                warning("Input not accepted")
            end
        end
        sweep = sweep + 1; % Continue to next sweep
        if sweep > size(data{proto,sensor_modality},1) % stop loop if max sweep reached
            loop = false; 
        end 
    end
else 
    fprintf('disable');
end 


%% TASK1.1: FC regression correlation with Soleus activity. (Seperate steps, single subject) 
% Find individuelle outliner og undersøg for refleks response. 

fprintf('\nscript: TASK1.1  . . . '); tic

%for subject = 1:10
show_plot = true;           % plot the following figures for this task
subject = 2;                % subject to analyse
proto = CTL;                % Only works for pre and post baseline trials
dep_sensory_modality = SOL; 
before = 100; 
after = 50;
xlimits = [-100 100];

% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
if show_plot
if ~readjust
    msg = "Be Aware: Data re-adjustment is disabled. Manuel defined search-bars applied . . . "; 
    fprintf(2,msg); 
else 
    fprintf(' Data re-adjustment enabled. General search bar applied . . .')
end

% General search bars: 
dep_off = 39; dep_len = 20; 
step = 2; predict_search(step,:) = [0,20]; depend_search(step,:) = [predict_search(step,1)+dep_off,predict_search(step,1)+dep_off+dep_len];    % ms 
step = 4; predict_search(step,:) = [0,20]; depend_search(step,:) = [predict_search(step,1)+dep_off,predict_search(step,1)+dep_off+dep_len];    % ms 
step = 6; predict_search(step,:) = [0,20]; depend_search(step,:) = [predict_search(step,1)+dep_off,predict_search(step,1)+dep_off+dep_len];    % ms 

% Preallocation
predictor_value = cell(size(names)); 
depended_value = cell(size(names)); 

for sub = 1:length(names)   % loop through subjects
    if ~readjust    % readjust disabled. Manuel search-bars applied
        switch sub 
        case 1
            predict_search(2,:) = [1.5,1.5+20];  depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-0.5,20-0.5]; depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-2,20-2];     depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms 
        case 2
            predict_search(2,:) = [-12,20-12];   depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [3,20-3];      depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];  % ms 
            predict_search(6,:) = [-4.5,20-4.5]; depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_off];   % ms 
        case 3
            predict_search(2,:) = [0.5, 20+0.5]; depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];   % ms 
            predict_search(4,:) = [1.5, 1.5+20]; depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [3, 3+20];     depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];   % ms 
        case 4
            predict_search(2,:) = [-4.5, 20-4.5]; depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-5.5, 20-5.5]; depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-6.5, 20-6.5]; depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms 
        case 5
            predict_search(2,:) = [2, 2+20];      depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];   % ms 
            predict_search(4,:) = [2.5, 2.5+20];  depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [2.5, 2.5+20];  depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];   % ms 
        case 6 
            predict_search(2,:) = [-2, 18];       depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [1, 21];        depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [0, 20];        depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms 
        case 7
            predict_search(2,:) = [-13.5, 20-13.5]; depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-15.5, 20-15.5]; depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-20.5, -0.5]; depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];   % ms 
        case 8
            predict_search(2,:) = [-20,0];        depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-9.5,20-9.5];  depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-9.5,20-9.5];  depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms  
        case 9
            predict_search(2,:) = [-53.5,20-53.5]; depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-7,20-7];       depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-41,20-41];     depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms  
        case 10
            predict_search(2,:) = [-69,20-69];     depend_search(2,:) = [predict_search(2,1)+dep_off,predict_search(2,1)+dep_off+dep_len];    % ms 
            predict_search(4,:) = [-56,20-56];     depend_search(4,:) = [predict_search(4,1)+dep_off,predict_search(4,1)+dep_off+dep_len];    % ms 
            predict_search(6,:) = [-61,20-61];      depend_search(6,:) = [predict_search(6,1)+dep_off,predict_search(6,1)+dep_off+dep_len];    % ms  
        end 
    end

    % remember values for later plt
    if sub == subject 
        predict_search_plt = predict_search; 
        depend_search_plt = depend_search; 
    end

    % load data from defined subject defined by loop
    data = total_data{1,1,sub};  
    step_index = total_step{1,1,sub};
    
    % find predictor and depended 
    for step = [2,4,6]                      % loop through steps
        for i = 1:size(data{proto,1},1)     % loop through sweeps    
            % from ms to sample
            predict_search_index = floor(msToSec(predict_search(step,:))*Fs);
            depend_search_index = floor(msToSec(depend_search(step,:))*Fs);
    
            % find rise index for the given step and define window 
            [rise, ~] = func_find_edge(step);
            rise_index = step_index{proto}(i,rise);
            predict_search_array = predict_search_index(1)+rise_index : predict_search_index(2)+rise_index; 
            depend_search_array = depend_search_index(1)+rise_index : depend_search_index(2)+rise_index; 
              
            % find average values
            %predictor_value(step,i) = mean(data{proto,pre_sensory_modality}((i), predict_search_array),2); 
            predictor_value{sub}(step, i) = (data{proto,ANG}((i),predict_search_array(1)) - data{proto,ANG}((i),predict_search_array(end))) / diff(predict_search(step,:)); 
            depended_value{sub}(SOL, step, i) = mean(data{proto,SOL}((i), depend_search_array),2);
            depended_value{sub}(TA, step, i) = mean(data{proto,TA}((i), depend_search_array),2);
        end 
    end
end 

% load data for plot, defined by >> subject <<
data = total_data{1,1,subject};  
step_index = total_step{1,1,subject};

screensize = get(0,'ScreenSize');
width = screensize(3);
height = screensize(4);

% plot figure
figSize = [50 50 width-200 height-200]; % where to plt and size
figure('Position', figSize); % begin plot 
    sgtitle("TASK 1.1 Data: Prebaseline. Subject: " + subject + ". [n = " + size(data{proto,1},1) + "]."); 

    predict_search_plt = predict_search_plt; 
    depend_search_plt = depend_search_plt; 

    % patch properties 
    y_pat = [-1000 -1000 1000 1000];
    patchcolor = "blue"; 
    FaceAlpha = 0.4; 
    
    for k = 1:3 % loop through steps       
        x_pat_pre = [predict_search_plt(steps_tested(k),1) predict_search_plt(steps_tested(k),2) predict_search_plt(steps_tested(k),2) predict_search_plt(steps_tested(k),1)];
        x_pat_dep = [depend_search_plt(steps_tested(k),1) depend_search_plt(steps_tested(k),2) depend_search_plt(steps_tested(k),2) depend_search_plt(steps_tested(k),1)];
    
        data_plot = cell(3,7); 
        [data_plot{proto,:}] = func_align(step_index{proto}, data{proto,[1:4,6:7]}, 'sec_before', msToSec(before), 'sec_after', msToSec(after), 'alignStep', align_with_obtions(k));

        subplot(5,3,0+k); hold on; % Ankel 
            % formalia setup
            ylabel("Position");
            title("Step " + k*2)
            subtitle("Ankel" +" "+ predict_search_plt(steps_tested(k),1) + " : "+predict_search_plt(steps_tested(k),2)+"ms")
            xlim(xlimits)

            % plt data 
            y = mean(data_plot{proto,ANG},1); 
            std_dev = std(data_plot{proto,ANG});
            curve1 = y + std_dev;
            curve2 = y - std_dev;
            x_axis = data_plot{proto,time};
            x_axis = secToMs(x_axis);
            x2 = [x_axis, fliplr(x_axis)];
            inBetween = [curve1, fliplr(curve2)];   
            fill(x2, inBetween, [0.75, 0.75, 0.75], 'LineStyle','none'); 
            plot(x_axis, y, 'LineWidth', 2, "color","black"); 
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat_pre,y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
    
        subplot(5,3,3+k); hold on; % velocity
            % formalia setup
            ylabel(labels_ms(VEL));
            subtitle(labels(VEL) +" "+ predict_search_plt(steps_tested(k),1) + " : "+predict_search_plt(steps_tested(k),2)+"ms")
            xlim(xlimits)
            
            % plt data 
            y = mean(data_plot{proto,VEL},1); 
            std_dev = std(data_plot{proto,VEL});
            curve1 = y + std_dev;
            curve2 = y - std_dev;
            x_axis = data_plot{proto,time};
            x_axis = secToMs(x_axis);
            x2 = [x_axis, fliplr(x_axis)];
            inBetween = [curve1, fliplr(curve2)];   
            fill(x2, inBetween, [0.75, 0.75, 0.75], 'LineStyle','none'); 
            plot(x_axis, y, 'LineWidth', 2, "color","black"); 
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat_pre,y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
            
        subplot(5,3,6+k); hold on; % soleus
            % formalia setup
            ylabel(labels_ms(SOL)); 
            subtitle(labels(SOL) +" "+ depend_search_plt(steps_tested(k),1) + " : "+depend_search_plt(steps_tested(k),2)+"ms")            
            xlim(xlimits)
    
            % plt data 
            y = mean(data_plot{proto,SOL},1); 
            std_dev = std(data_plot{proto,SOL});
            curve1 = y + std_dev;
            curve2 = y - std_dev;
            x_axis = data_plot{proto,time};
            x_axis = secToMs(x_axis);
            x2 = [x_axis, fliplr(x_axis)];
            inBetween = [curve1, fliplr(curve2)];   
            fill(x2, inBetween, [0.75, 0.75, 0.75], 'LineStyle','none'); 
            plot(x_axis, y, 'LineWidth', 2, "color","black"); 
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat_dep,y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")

        subplot(5,3,9+k); hold on; % tibíalis
             % formalia setup
            ylabel(labels_ms(TA)); 
            subtitle(labels(TA) +" "+ depend_search_plt(steps_tested(k),1) + " : "+depend_search_plt(steps_tested(k),2)+"ms")            
            xlim(xlimits)
    
            % plt data 
            y = mean(data_plot{proto,TA},1); 
            std_dev = std(data_plot{proto,TA});
            curve1 = y + std_dev;
            curve2 = y - std_dev;
            x_axis = data_plot{proto,time};
            x_axis = secToMs(x_axis);
            x2 = [x_axis, fliplr(x_axis)];
            inBetween = [curve1, fliplr(curve2)];   
            fill(x2, inBetween, [0.75, 0.75, 0.75], 'LineStyle','none'); 
            plot(x_axis, y, 'LineWidth', 2, "color","black"); 
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat_dep,y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")

        subplot(5,3,12+k); hold on; % plot regression 
            % formalia 
            xlabel("Pos(start-end)/window size")
            ylabel("Avg. "+labels(dep_sensory_modality))
            
            % plt data 
            depended = []; predictor = [];  
            depended = nonzeros(squeeze(depended_value{subject}(dep_sensory_modality, steps_tested(k),:))); 
            predictor = nonzeros(squeeze(predictor_value{subject}(steps_tested(k),:)));
            mdl = fitlm(predictor, depended);
            b = table2array(mdl.Coefficients(1,1)); 
            a = table2array(mdl.Coefficients(2,1)); 
            p_value = table2array(mdl.Coefficients(2,4)); 
            linearReg = @(x) x*a + b;     
            plot(predictor, depended,"x","color", "blue")          
            plot(predictor, linearReg(predictor), "color", "red")
            subtitle("Sweep average" + newline + "p-value: " + round(p_value,2))
    end 

    % Define the file name and path to save the PNG file
    filename = "subject"+subject+"-"+labels(dep_sensory_modality)+".png";
    filepath = 'C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/png files/task1 - EMG vs Pos(start-end) - Readjusted/';
    fullpath = fullfile(filepath, filename);
    saveas(gcf, fullpath, 'png');
fprintf('done [ %4.2f sec ] \n', toc);
else
fprintf('disable \n');
end

%% TASK1.2: FC regression correlation with Soleus activity (All Step, single subject) 
% Same task as task1, but all steps considered as the same. 

fprintf('script: TASK1.2  . . . '); tic
show_plot = true; 


if show_plot
marker = ["*",".","x"];    
color = [[0 0 1];[0.5 0 0.5];[1 .1 0]]; 

figure; % begin plot
sgtitle("TASK 1.2 Data: Prebaseline. Subject: " + subject); 

for sensory_type = [SOL,TA]
    subplot(1,2,sensory_type); hold on

    depended = []; predictor = [];  
    for k = 1:3
        depended(k,:) = nonzeros(squeeze(depended_value{subject}(sensory_type, steps_tested(k),:))); 
        predictor(k,:) = nonzeros(squeeze(predictor_value{subject}(steps_tested(k),:)));
        plot(predictor(k,:), depended(k,:), marker(k),"color", "blue")
    end
    
    % linear regression 
    depended_all_step = [depended(1,:) depended(2,:) depended(3,:)];
    predictor_all_step = [predictor(1,:) predictor(2,:) predictor(3,:)];
    mdl = fitlm(predictor_all_step, depended_all_step);   % <--- SIG
    b = table2array(mdl.Coefficients(1,1)); 
    a = table2array(mdl.Coefficients(2,1));
    p_value = table2array(mdl.Coefficients(2,4)); 
    linearReg = @(x) x*a + b; 
    plot(predictor_all_step, linearReg(predictor_all_step), "color", "red")
    
    % plt formalia 
    legend(["Data: Step 2", "Data: Step 4", "Data: Step 6","Fit"])
    xlabel("Pos(start-end)/window size")
    ylabel("Avg. " + labels(sensory_type))
    title("p-value: " + round(p_value,2))
end 

filename = "subject"+subject+"-all steps.png";
filepath = 'C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/png files/task1 - EMG vs Pos(start-end) - Readjusted/';
fullpath = fullfile(filepath, filename);
saveas(gcf, fullpath, 'png');

fprintf('done [ %4.2f sec ] \n', toc);
else
fprintf('disable \n');
end


%% TASK1.3: FC regression correlation with EMG (Seperate steps, All subject) 
fprintf('script: TASK1.3  . . .'); tic
show_plot = true; 

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
if show_plot
    marker = ["*",".","x"];    
    color = [[0 0 1];[0.5 0 0.5];[1 .1 0]]; 

    % re-arrange data from cell to struct. 
    data_reg = struct; 
    data_reg.dep_steps = cell(2,3);     % depended sortet each step, [sol,ta]
    data_reg.pre_steps = cell(1,3);     % predictor sortet each step, [vel]
    data_reg.pre = [];                  % depended sortet all step, [sol,ta]
    data_reg.dep = cell(2,1);           % predictor sortet all step, [vel]

    % re-arrange data from cell to array and save in struct
    for sub = 1:numel(names)
        for step = 1:3 
            for EMG = [SOL,TA]
                data_reg.dep_steps{EMG,step} = [data_reg.dep_steps{EMG,step}, nonzeros(squeeze(depended_value{sub}(EMG, steps_tested(step),:)))' ]; 
                data_reg.dep{EMG} = [data_reg.dep{EMG}, nonzeros(squeeze(depended_value{sub}(EMG, steps_tested(step),:)))' ];
            end
            data_reg.pre = [data_reg.pre,  nonzeros(squeeze(predictor_value{sub}(steps_tested(step),:)))'];
            data_reg.pre_steps{step} = [data_reg.pre_steps{step}, nonzeros(squeeze(predictor_value{sub}(steps_tested(step),:)))'];
        end 
    end 

    % defining plt window size
    screensize = get(0,'ScreenSize');
    width = screensize(3);
    height = screensize(4);
    figSize = [100 100 width-200 height-200]; % where to plt and size
    figure('Position', figSize); % begin plot 
    sgtitle("TASK 1.3 All subjects (Subjects displayed in different colors)"); 
    
    
    % plt subject data in different colors 
    for sub = 1:numel(names)
        for sensory_type = [SOL,TA] % loop 
            depended = []; predictor = [];  
            for k = 1:3 % loop steps 
                depended(k,:) = nonzeros(squeeze(depended_value{sub}(sensory_type, steps_tested(k),:))); 
                predictor(k,:) = nonzeros(squeeze(predictor_value{sub}(steps_tested(k),:)));
                
                % plt the individuel steps
                subplot(2, 5, 5*(sensory_type-1)+k); hold on; % 5*(s-1)+k={1,2,3,6,7,8}, k={1,2,3}, s={1,2}            
                plot(predictor(k,:), depended(k,:),'x');
                title("Step " + steps_tested(k))
                xlabel("Pos(s-e)/w")
                ylabel(labels(sensory_type))
            end

            % plt the combined steps 
            subplot(2,5, 4+5*(sensory_type-1):5+5*(sensory_type-1)); hold on % 4+5*(s-1)={4,9}, 5+5*(s-1)={5,10} s={1,2}
            plot(predictor(:)', depended(:)','x')
            title("All steps")
            subtitle("P-value " + "R^2")
            ylabel(labels(sensory_type))
            xlabel("Pos(s-e)/w")
        end
    end 

    % set y-limits and x-limits on all subplots 
    ax = findobj(gcf, 'type', 'axes');
    ylims = get(ax, 'YLim'); 
    xlims = get(ax, 'XLim'); 
    [~, idx_y_ta] = max(cellfun(@(x) diff(x), ylims(1:4)));
    [~, idx_y_sol] = max(cellfun(@(x) diff(x), ylims(5:8)));
    [~, idx_x_ta] = max(cellfun(@(x) diff(x), xlims(1:4)));
    [~, idx_x_sol] = max(cellfun(@(x) diff(x), xlims(5:8)));
    idx_y_sol = idx_y_sol+4; 
    idx_x_sol = idx_x_sol+4; 
    for i = 1:numel(ax)
        if any(i == 1:4)
            set(ax(i), 'YLim', ylims{idx_y_ta})
            set(ax(i), 'XLim', xlims{idx_x_ta})
        elseif any(i == 5:8)  
            set(ax(i), 'YLim', ylims{idx_y_sol})
            set(ax(i), 'XLim', xlims{idx_x_sol})
        end
    end
    
    % plot regression 
    for sensory_type = [SOL,TA] % loop 
        for step = 1:3 % loop steps  
            mdl = fitlm(data_reg.pre_steps{step}, data_reg.dep_steps{sensory_type, step}); 
            b = table2array(mdl.Coefficients(1,1)); 
            a = table2array(mdl.Coefficients(2,1)); 
            p_value = table2array(mdl.Coefficients(2,4)); 
            r2 = mdl.Rsquared.Adjusted; 
            linearReg = @(x) x*a + b; 
            subplot(2, 5, 5*(sensory_type-1)+step); hold on; 
            plot(data_reg.pre_steps{step}, linearReg(data_reg.pre_steps{step}), "color", "red")
            subtitle("P-value " + p_value + ". R^2 " + r2)
        end

        mdl = fitlm(data_reg.pre, data_reg.dep{sensory_type}); 
        b = table2array(mdl.Coefficients(1,1)); 
        a = table2array(mdl.Coefficients(2,1)); 
        p_value = table2array(mdl.Coefficients(2,4)); 
        r2 = mdl.Rsquared.Adjusted; 
        linearReg = @(x) x*a + b; 
        subplot(2,5, 4+5*(sensory_type-1):5+5*(sensory_type-1)); hold on
        plot(data_reg.pre_steps{step}, linearReg(data_reg.pre_steps{step}), "color", "red")
        subtitle("P-value " + p_value + ". R^2 " + r2)
    end


filename = "all subject (1-8)";
filepath = 'C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/png files/task1 - EMG vs Pos(start-end) - Readjusted/';
fullpath = fullfile(filepath, filename);
saveas(gcf, fullpath, 'png');

fprintf('done [ %4.2f sec ] \n', toc);
else
fprintf('disable \n');
end

%% Task1.5 FC regression correlation with Soleus activity (steepest ascent)
% Find the best parameters for one subject and apply them to the other
%    subjects. 
% 

%% TASK2.1: Within step adjustment of EMG acticity due to natural angle variation (single subject)
% Find whether the EMG activity is passively adjusted to changes in angle
%    trajectories. 
% All steps and seperated step 
% nassarro
% two plot: inden burst activitet og til max burst activitet (all subject)
fprintf('script: TASK2.1  . . . '); tic

show_plot = true;
subject = 1; % 
protocol = CTL; 
before = 100; % [ms]
after = 0; % [ms]

pre = 1; search_area(pre,:) = [0.15; 0.60]; % predictor search window
dep = 2; search_area(dep,:) = [0.15; 0.70]; % depended search window

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
% load data from subject
data = total_data{1,1,subject};  
step_index = total_step{1,1,subject};
        
if show_plot
    figSize = [200 60 1000 700];
    figure('Position', figSize)
    sgtitle("Within step modulation due to natural angle variation. Subject "+ subject + newline + "{\color{blue}Graph} larger than avg. and {\color{red}Graph} less than avg.")
    
    clear win_avg_ang win_avg_sol win_avg_ta
    for step = 1:3
        % align data needed to be plotted
        clear temp_plot
        temp_plot = cell(3,7); 
        [temp_plot{protocol,:}] = func_align(step_index{protocol}, data{proto,[1:4,6:7]}, 'sec_before', msToSec(before), 'sec_after', msToSec(after), 'alignStep', align_with_obtions(step));
        x_axis = temp_plot{protocol, time};
          
        % align data needed to calculate window avg
        clear temp_data
        temp_data = cell(3,7); 
        [temp_data{protocol,:}] = func_align(step_index{protocol}, data{proto,[1:4,6:7]}, 'alignStep', align_with_obtions(step));
    
        % Defined window  
        for cc = [dep, pre]
            pct_start_index(cc) = ceil(size(temp_data{protocol,time},2)*search_area(cc,1));   % window begin in idx
            pct_end_index(cc)   = ceil(size(temp_data{protocol,time},2)*search_area(cc,2));   % window end in idx
            pct_start_sec(cc) = temp_data{protocol,time}(1,pct_start_index(cc));              % window begin in sec
            pct_end_sec(cc)   = temp_data{protocol,time}(1,pct_end_index(cc));                % window end in sec
        end 

        % Calculate average in window 
        win_avg_ang(step,:) = mean(temp_data{protocol,ANG}(:,pct_start_index(pre) : pct_end_index(pre)),2); % mean value of sweeps window   
        win_avg_sol(step,:) = mean(temp_data{protocol,SOL}(:,pct_start_index(dep):pct_end_index(dep)),2);   % mean value of sweeps window       
        win_avg_ta(step,:) = mean(temp_data{protocol,TA}(:,pct_start_index(dep):pct_end_index(dep)),2);     % mean value of sweeps window       
        
        % Sort data according 
        [win_avg_ang(step,:), win_avg_idx] = sort(win_avg_ang(step,:), 'ascend'); % sorts elements in ascending order.
        win_avg_sol(step,:) = win_avg_sol(step, win_avg_idx); % sort these with found order
        win_avg_ta(step,:) = win_avg_ta(step, win_avg_idx);   % sort these with found order

        clear lower_index mean_index upper_index
        % sort data in groups [low, mid, up] 
        dim = size(win_avg_ang,2); % num of sweeps
        lower_index(:) = 1:floor(dim/3); 
        mean_index(:) = floor(dim/3)+1:dim-floor(dim/3); 
        upper_index(:) = dim-floor(dim/3)+1:dim; 

        % mean plot properties
        color_mean = [150 152 158]/255;
        LineWidth_mean = 2; 
    
        % upper and lower plot properties
        color_upper = "blue";
        color_lower = "red";
        upper_style = "--";
        lower_style = "-.";
        LineWidth_UL = 1; 
    
        % patch properties
        y_pat = [-1000 -1000 2000 2000];
        patchcolor = [251 244 199]/255; 
        FaceAlpha = 1; 
        EdgeColor = [37 137 70]/255;
        lineWidth_patch = 2; 
        for cc = [dep, pre]
            x_pat(cc,:) = [pct_start_sec(cc), pct_end_sec(cc), pct_end_sec(cc), pct_start_sec(cc)]; 
        end 

        % plot figures 
        % ankel
        subplot(330+step); hold on; 
        title("Step " + steps_tested(step))
        subtitle("Ankel"); 
        ylabel(labels_sec(ANG))
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(lower_index),:),1), lower_style,'LineWidth', LineWidth_UL, "color", color_lower)
            YL = get(gca, 'YLim'); 
            ylim([YL(1) YL(2)]);
            patch(x_pat(pre,:),y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', EdgeColor, 'LineWidth', lineWidth_patch)
            set(gca, 'Layer', 'top')
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, ANG}(win_avg_idx(lower_index),:),1), lower_style,'LineWidth', LineWidth_UL, "color", color_lower)
    
    
        % soleus
        subplot(333+step); hold on; 
        subtitle(labels(SOL))
        ylabel("Soleus"+newline+"[normalized]")
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(lower_index),:),1), lower_style,'LineWidth', LineWidth_UL, "color", color_lower)
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat(dep,:),y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', EdgeColor, 'LineWidth', lineWidth_patch)
            set(gca, 'Layer', 'top')
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, SOL}(win_avg_idx(lower_index),:),1), lower_style,'LineWidth', LineWidth_UL, "color", color_lower)
    
        % tibialis
        subplot(336+step); hold on;
        subtitle(labels(TA))
        ylabel("Tibialis"+newline+"[normalized]")
        xlabel(labels_sec(time))
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(lower_index),:),1), lower_style,'LineWidth', LineWidth_UL, "color", color_lower)
            YL = get(gca, 'YLim'); ylim([YL(1) YL(2)]);
            patch(x_pat(dep,:),y_pat,patchcolor,'FaceAlpha',FaceAlpha, 'EdgeColor', EdgeColor, 'LineWidth', lineWidth_patch)
            set(gca, 'Layer', 'top')
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(mean_index),:),1), 'LineWidth', LineWidth_mean, "color", color_mean)
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(upper_index),:),1), upper_style, 'LineWidth', LineWidth_UL, "color", color_upper)
            plot(x_axis, mean(temp_plot{protocol, TA}(win_avg_idx(lower_index),:),1), lower_style, 'LineWidth', LineWidth_UL, "color", color_lower)
    end 
fprintf('done [ %4.2f sec ] \n', toc);
else 
fprintf('disable \n');
end   

%% task 2.2 Within step adjustment of EMG acticity due to natural angle variation. (single subject)
fprintf('script: TASK2.2  . . . '); tic

show_plt = true; 

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
if show_plt
    clear x_sol x_ta
    % Rearrange data to box plot data
    step2 = 1; step4 = 2; step6 = 3; 
    % soleus
    group_low_sol = [win_avg_sol(step2,lower_index) ; win_avg_sol(step4,lower_index) ; win_avg_sol(step6,lower_index)]';
    group_mid_sol = [win_avg_sol(step2,mean_index) ; win_avg_sol(step4,mean_index) ; win_avg_sol(step6,mean_index)]';
    group_up_sol = [win_avg_sol(step2,upper_index) ; win_avg_sol(step4,upper_index) ; win_avg_sol(step6,upper_index)]';
    x_sol = {group_low_sol, group_mid_sol, group_up_sol};
    % tibialis
    group_low_ta = [win_avg_ta(step2,lower_index) ; win_avg_ta(step4,lower_index) ; win_avg_ta(step6,lower_index)]';
    group_mid_ta = [win_avg_ta(step2,mean_index) ; win_avg_ta(step4,mean_index) ; win_avg_ta(step6,mean_index)]';
    group_up_ta = [win_avg_ta(step2,upper_index) ; win_avg_ta(step4,upper_index) ; win_avg_ta(step6,upper_index)]';
    x_ta = {group_low_ta, group_mid_ta, group_up_ta};

    % Plot data to get natural y limits
    figure(13); 
    boxplotGroup(x_sol);  
    YL_sol = get(gca, 'YLim'); 
    close 13

    figure(13); 
    boxplotGroup(x_ta);  
    YL_ta = get(gca, 'YLim'); 
    close 13
    
    % Plot boxplot for single subject
    figure; hold on; 
        blue = 	[0 0 1]; red = [1 0 0]; gray = color_mean; 

    subplot(211)
        grpLabels = {'Step 2', 'Step 4', 'Step 6'}; 
        sublabels = {'low', 'mid', 'up'};
        boxplotGroup(x_sol,'primaryLabels',sublabels,'SecondaryLabels',grpLabels, 'interGroupSpace',2,'GroupLines',true, 'Colors',[red; gray; blue],'GroupType','betweenGroups')
        ylim([YL_sol(1) YL_sol(2)])
        ylabel("Normalized muscle activity"+newline+"Soleus")
    

    % Plot boxplot for single subject
    subplot(212);hold on; 
        title("Within step adjustment of EMG acticity due to natural angle variation",'FontName','FixedWidth')
        subtitle("Single Subject. Num: " + subject ,'FontName','FixedWidth')
        boxplotGroup(x_ta,'primaryLabels',sublabels,'SecondaryLabels',grpLabels, 'interGroupSpace',2,'GroupLines',true, 'Colors',[red; gray; blue],'GroupType','betweenGroups')
        ylim([YL_ta(1) YL_ta(2)])
        ylabel("Normalized muscle activity"+newline+"Tibialis")

fprintf('done [ %4.2f sec ] \n', toc);
else 
fprintf('disable \n');
end  

%% task2.3 Within step adjustment of EMG acticity due to natural angle variation. (All subject)
fprintf('script: TASK2.3  . . . '); tic

show_plt = true; 

%  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
low = 1; mid = 2; up = 3; 
step2 = 1; step4 = 2; step6 = 3; 

if show_plt
    for sub = 1:numel(names) % loop through names 
        for step = 1:3 % loop through steps
            % align data to find window avg
            clear data step_index
            data = total_data{1,1,sub};  
            step_index = total_step{1,1,sub};
            
            clear temp_data
            temp_data = cell(3,7); 
            [temp_data{protocol,:}] = func_align(step_index{protocol}, data{proto,[1:4,6:7]}, 'alignStep', align_with_obtions(step));
        
            % Defined window  
            for cc = [dep, pre]
                pct_start_index(cc) = ceil(size(temp_data{protocol,time},2)*search_area(cc,1));   % window begin in idx
                pct_end_index(cc)   = ceil(size(temp_data{protocol,time},2)*search_area(cc,2));   % window end in idx
                pct_start_sec(cc) = temp_data{protocol,time}(1,pct_start_index(cc));              % window begin in sec
                pct_end_sec(cc)   = temp_data{protocol,time}(1,pct_end_index(cc));                % window end in sec
            end 
    
            clear win_avg_ang win_avg_sol win_avg_ta
            % Calculate average in window 
            win_avg_ang(step,:) = mean(temp_data{protocol,ANG}(:,pct_start_index(pre) : pct_end_index(pre)),2); % mean value of sweeps window   
            win_avg_sol(step,:) = mean(temp_data{protocol,SOL}(:,pct_start_index(dep):pct_end_index(dep)),2);   % mean value of sweeps window       
            win_avg_ta(step,:) = mean(temp_data{protocol,TA}(:,pct_start_index(dep):pct_end_index(dep)),2);     % mean value of sweeps window       
        
             % Sort data according
            [win_avg_ang(step,:), win_avg_idx] = sort(win_avg_ang(step,:), 'ascend'); % sorts elements in ascending order.
            win_avg_sol(step,:) = win_avg_sol(step, win_avg_idx); % sort these with found order
            win_avg_ta(step,:) = win_avg_ta(step, win_avg_idx);   % sort these with found order
        
            clear lower_index mean_index upper_index
            % Sort data in groups [low, mid, up] 
            dim = size(win_avg_ang,2); % num of sweeps
            lower_index(:) = 1:floor(dim/3); 
            mean_index(:) = floor(dim/3)+1:dim-floor(dim/3); 
            upper_index(:) = dim-floor(dim/3)+1:dim; 
        
            % Subject mean of each subgroup
            subgroup_sol(sub, step, low) = mean(win_avg_sol(step,lower_index));
            subgroup_sol(sub, step, mid) = mean(win_avg_sol(step,mean_index));
            subgroup_sol(sub, step, up)  = mean(win_avg_sol(step,upper_index));
            subgroup_ta(sub, step, low)  = mean(win_avg_ta(step,lower_index));
            subgroup_ta(sub, step, mid)  = mean(win_avg_ta(step,mean_index));
            subgroup_ta(sub, step, up)   = mean(win_avg_ta(step,upper_index));
        end % step
    end % sub 
    
    % Rearrange data to box plot data
    step2 = 1; step4 = 2; step6 = 3; 
    group_low_sol = [subgroup_sol(:,step2,low)' ; subgroup_sol(:,step4,low)' ; subgroup_sol(:,step6,low)' ; (subgroup_sol(:,step2,low)'+subgroup_sol(:,step4,low)'+subgroup_sol(:,step6,low)')./3]';
    group_mid_sol = [subgroup_sol(:,step2,mid)' ; subgroup_sol(:,step4,mid)' ; subgroup_sol(:,step6,mid)' ; (subgroup_sol(:,step2,mid)'+subgroup_sol(:,step4,mid)'+subgroup_sol(:,step6,mid)')./3]';
    group_up_sol =  [subgroup_sol(:,step2,up)' ; subgroup_sol(:,step4,up)' ; subgroup_sol(:,step6,up)' ; (subgroup_sol(:,step2,up)'+subgroup_sol(:,step4,up)'+subgroup_sol(:,step6,up)')./3]';
    x_sol = {group_low_sol, group_mid_sol, group_up_sol};
    
    group_low_ta = [subgroup_ta(:,step2,low)' ; subgroup_ta(:,step4,low)' ; subgroup_ta(:,step6,low)'; (subgroup_ta(:,step2,low)'+subgroup_ta(:,step4,low)'+subgroup_ta(:,step6,low)')./3]';
    group_mid_ta = [subgroup_ta(:,step2,mid)' ; subgroup_ta(:,step4,mid)' ; subgroup_ta(:,step6,mid)' ; (subgroup_ta(:,step2,mid)'+subgroup_ta(:,step4,mid)'+subgroup_ta(:,step6,mid)')./3]';
    group_up_ta =  [subgroup_ta(:,step2,up)' ; subgroup_ta(:,step4,up)' ; subgroup_ta(:,step6,up)'; (subgroup_ta(:,step2,up)'+subgroup_ta(:,step4,up)'+subgroup_ta(:,step6,up)')./3]';
    x_ta = {group_low_ta, group_mid_ta, group_up_ta};

    x_all = {(group_low_sol+group_low_ta)./2, (group_mid_sol+group_mid_ta)./2, (group_up_sol+group_up_ta)./2 } ;


    % Plot boxplot for single subject
    figure; hold on; 
    blue = 	[0 0 1]; red = [1 0 0]; gray = [150 152 158]/255;

    subplot(311)
        title("Within step adjustment - All Subject",'FontName','FixedWidth')
        grpLabels = {'Step 2', 'Step 4', 'Step 6', 'All steps'}; 
        sublabels = {'low', 'mid', 'up'};
        %plot([0 4],[100 200])
        boxplotGroup(x_sol,'primaryLabels',sublabels,'SecondaryLabels',grpLabels, 'interGroupSpace',2,'GroupLines',true, 'Colors',[red; gray; blue],'GroupType','betweenGroups')
        %ylim([YL_sol(1) YL_sol(2)])
        ylim([0 0.5])
        ylabel("Normalized muscle activity"+newline+"Soleus")
    

    % Plot boxplot for single subject
    subplot(312);hold on; 
        boxplotGroup(x_ta,'primaryLabels',sublabels,'SecondaryLabels',grpLabels, 'interGroupSpace',2,'GroupLines',true, 'Colors',[red; gray; blue],'GroupType','betweenGroups')
        %ylim([YL_ta(1) YL_ta(2)]
        ylim([0 0.5])
        ylabel("Normalized muscle activity"+newline+"Tibialis")

    subplot(313);hold on; 
        boxplotGroup(x_all,'primaryLabels',sublabels,'SecondaryLabels',grpLabels, 'interGroupSpace',2,'GroupLines',true, 'Colors',[red; gray; blue],'GroupType','betweenGroups')
        ylim([0 0.5])
        ylabel("Normalized muscle activity"+ newline +"for all muscles")


    % Verify code
    % y = [mean(group_low); mean(group_mid); mean(group_up)];
    % % Bar plot 
    % figure; 
    % gray = [226,226,226]/255;
    % bar(y, 'FaceColor', gray)

fprintf('done [ %4.2f sec ] \n', toc);
else 
fprintf('disable \n');
end %show_plt




%% Task3.1 horizontal perturbation. 
% plot difference i hastighed ift soleus aktivitet som regression plot. 
% 5 første stræk mod 5 sidste stræk 

show_plt = true; 
subject = 2; 
before = 500; 
after = 0; 
xlimit = [-100 200];



clear offset 
offset(1) = 30; 
offset(2) = 13; %38;
offset(3) = 19;
offset(4) = 0;
offset(5) = 10;  %31;
offset(6) = 36;
offset(7) = 40;
offset(8) = 38;


if show_plt
    data = total_data{1,1,subject};  
    step_index = total_step{1,1,subject};
    type = total_type{1,1,subject}; yes = type{3}; no = type{4}; 

    clear temp_plot
    temp_plot = cell(3,7); 
    [temp_plot{HOR,:}] = func_align(step_index{HOR}, data{HOR,[1:4,6:7]}, 'sec_before', msToSec(before), 'sec_after', msToSec(after), 'alignStep', "four_begin");
    x_axis = secToMs(temp_plot{HOR, time});

    %plot properties 
    no_color = [0.75, 0.75, 0.75]; 
    yes_color = "black";
    zero_color = "red"; 
    no_LineWidth = 3; 
    yes_LineWidth = 1; 
    zero_LineWidth = 1; 


    % patch properties
    y_pat = [-1000 -1000 2000 2000];
    patchcolor_slr = "blue" ; %[251 244 199]/255; 
    patchcolor_mlr = "black" ; %[251 244 199]/255; 

    FaceAlpha = 0.2; 
    EdgeColor = [37 137 70]/255;
    lineWidth_patch = 2; 
    x_pat_SLR = [offset(subject)+39 offset(subject)+59 offset(subject)+59 offset(subject)+39];
    x_pat_MLR = [offset(subject)+59 offset(subject)+79 offset(subject)+79 offset(subject)+59];


    % Check if a figure with the name 'TASK3' is open
    fig = findobj('Name', 'TASK3');
    % If a figure is found, close it
    if ~isempty(fig), close(fig); end

    % Create a new figure with the name 'TASK3'
   % figure('Name','TASK3')
    figSize = [200 60 1000 700];
    figure('Name','TASK3','Position', figSize)
    sgtitle("Horizontal. Subject " + subject)
  
    sensortype = ANG; % position 
    subplot(411); hold on; 
        % plot formalia (411)       
        title(['{\color{gray}Control sweeps [n=' num2str(length(no)) '].} Perturbation sweeps [n=' num2str(length(yes)) '].  {\color{red} Perturbation onset ' num2str(offset(subject)) ' [ms] }' ])
        subtitle(labels(sensortype))
        ylabel(labels_ms(sensortype))
        if ~isempty(xlimit), xlim(xlimit); end 
                
        % plot data (411)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no,:),1), 'LineWidth',no_LineWidth, 'Color',no_color)
                plot(x_axis , mean(temp_plot{HOR,sensortype}(no(1:5),:),1), 'LineWidth',no_LineWidth, 'Color',"blue")

        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes,:),1), 'LineWidth',yes_LineWidth, 'Color',yes_color)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes(1:5),:),1), 'LineWidth',yes_LineWidth, 'Color',"red")

        YL = get(gca, 'YLim'); ylim([YL(1) YL(2)])        
        plot([offset(subject) offset(subject)],[-100 100], 'lineWidth', zero_LineWidth, 'Color',zero_color)

    
    sensortype = VEL; % velocity
    subplot(412); hold on; 
        % plot formalia (412)
        subtitle(labels(sensortype))
        ylabel(labels_ms(sensortype))
        if ~isempty(xlimit), xlim(xlimit); end 
        grid on

        % plot data (412)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no,:),1), 'LineWidth',no_LineWidth, 'Color',no_color)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no(1:5),:),1), 'LineWidth',no_LineWidth, 'Color',"blue")

        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes,:),1), 'LineWidth',yes_LineWidth, 'Color',yes_color)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes(1:5),:),1), 'LineWidth',yes_LineWidth, 'Color',"red")

        YL = get(gca, 'YLim'); ylim([YL(1) YL(2)])
        plot([offset(subject) offset(subject)],[-100 100], 'lineWidth', zero_LineWidth, 'Color',zero_color)

    
    sensortype = SOL; % soleus 
    subplot(413); hold on; 
        % plot formalia (413)
        subtitle(labels(SOL))
        ylabel(labels(sensortype))
        
        % plot data (413)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no,:),1), 'LineWidth',no_LineWidth, 'Color',no_color)
                plot(x_axis , mean(temp_plot{HOR,sensortype}(no(1:5),:),1), 'LineWidth',no_LineWidth, 'Color',"Blue")

        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes,:),1), 'LineWidth',yes_LineWidth, 'Color',yes_color)
                plot(x_axis , mean(temp_plot{HOR,sensortype}(yes(1:5),:),1), 'LineWidth',yes_LineWidth, 'Color',"red")

        YL = get(gca, 'YLim'); ylim([YL(1) YL(2)])
        plot([offset(subject) offset(subject)],[-100 100], 'lineWidth', zero_LineWidth, 'Color',zero_color)
        if ~isempty(xlimit), xlim(xlimit); end 
        patch(x_pat_SLR,y_pat,patchcolor_slr,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
        patch(x_pat_MLR,y_pat,patchcolor_mlr,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")


    sensortype = TA; % tibialis 
    subplot(414); hold on;
        % plot formalia (414)
        subtitle(labels(sensortype))
        ylabel(labels(sensortype))
        xlabel(labels_ms(time))
    
        % plot data (414)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no,:),1), 'LineWidth',no_LineWidth, 'Color',no_color)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(no(1:5),:),1), 'LineWidth',no_LineWidth, 'Color',"blue")

        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes,:),1), 'LineWidth',yes_LineWidth, 'Color',yes_color)
        plot(x_axis , mean(temp_plot{HOR,sensortype}(yes(1:5),:),1), 'LineWidth',yes_LineWidth, 'Color',"red")

        YL = get(gca, 'YLim'); ylim([YL(1) YL(2)])
        plot([offset(subject) offset(subject)],[-100 100], 'lineWidth', zero_LineWidth, 'Color',zero_color)
        if ~isempty(xlimit), xlim(xlimit); end 
        patch(x_pat_SLR,y_pat,patchcolor_slr,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")
        patch(x_pat_MLR,y_pat,patchcolor_mlr,'FaceAlpha',FaceAlpha, 'EdgeColor', "none")

%     filename = "Horizontal. Subject "+subject+".png";
%     filepath = 'C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/png files/Task3 - Horizontal perturbation/';
%     fullpath = fullfile(filepath, filename);
%     saveas(gcf, fullpath, 'png');
end 


%% TASK3.2 Horizontal perturbation boxplot 
% Show boxplot where all subject can clearly be identified. 
show_plt = false; 

cnt = 0; 
SLR = 1; MLR = 2; 
if show_plt
    clear avg_CTL avg_CTL
    for sub = [1,2,3,5,6,7,8]
        cnt = cnt + 1; 
        data = total_data{1,1,sub};  
        step_index = total_step{1,1,sub};
        type = total_type{1,1,sub}; yes = type{3}; no = type{4}; 
        
        clear temp_data
        temp_data = cell(3,7); 
        [temp_data{HOR,:}] = func_align(step_index{HOR}, data{HOR,[1:4,6:7]}, 'alignStep', "four_begin");
        [temp_data{CTL,:}] = func_align(step_index{CTL}, data{CTL,[1:4,6:7]}, 'alignStep', "four_begin");

        SLR_begin = floor(msToSec(offset(sub)+39)*Fs);
        SLR_end = floor(msToSec(offset(sub)+59)*Fs);
        MLR_begin = floor(msToSec(offset(sub)+60)*Fs);
        MLR_end = floor(msToSec(offset(sub)+80)*Fs);

        avg_CTL(SLR,cnt,SOL) = mean(temp_data{CTL,SOL}(:,[SLR_begin : SLR_end]),[1:2]);
        avg_CTL(MLR,cnt,SOL) = mean(temp_data{CTL,TA}(:,[SLR_begin : SLR_end]),[1:2]);
        avg_CTL(SLR,cnt,TA) = mean(temp_data{CTL,SOL}(:,[MLR_begin : MLR_end]),[1:2]);
        avg_CTL(MLR,cnt,TA) = mean(temp_data{CTL,TA}(:,[MLR_begin : MLR_end]),[1:2]);

        avg_hor(SLR,cnt,SOL) = mean(temp_data{HOR,SOL}(:,[SLR_begin : SLR_end]),[1:2]);
        avg_hor(MLR,cnt,SOL) = mean(temp_data{HOR,TA}(:,[SLR_begin : SLR_end]),[1:2]);
        avg_hor(SLR,cnt,TA) = mean(temp_data{HOR,SOL}(:,[MLR_begin : MLR_end]),[1:2]);
        avg_hor(MLR,cnt,TA) = mean(temp_data{HOR,TA}(:,[MLR_begin : MLR_end]),[1:2]);
    end 

   group_control_sol = [avg_CTL(SLR,:,SOL)' ; avg_CTL(MLR,:,SOL)'];
   group_control_sol = [avg_CTL(SLR,:,SOL)' ; avg_CTL(MLR,:,SOL)'];

   group_control_ta = [avg_CTL(SLR,:,TA)' ; avg_CTL(MLR,:,TA)'];

%     group_mid_sol = [subgroup_sol(:,step2,mid)' ; subgroup_sol(:,step4,mid)' ; subgroup_sol(:,step6,mid)' ; (subgroup_sol(:,step2,mid)'+subgroup_sol(:,step4,mid)'+subgroup_sol(:,step6,mid)')./3]';
%     group_up_sol =  [subgroup_sol(:,step2,up)' ; subgroup_sol(:,step4,up)' ; subgroup_sol(:,step6,up)' ; (subgroup_sol(:,step2,up)'+subgroup_sol(:,step4,up)'+subgroup_sol(:,step6,up)')./3]';
%     x_sol = {group_low_sol, group_mid_sol, group_up_sol};

end 


%% TASK4 Pre-baseline vs Post-baseline 
% Does the spinal influence chance due to the experienced protocols. 

%% TASK5 Show individual Unload trials


%% TASK6 make foot movement graph





%% finsihed
fprintf('\n\n Processed finished \n') 
