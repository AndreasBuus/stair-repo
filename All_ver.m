clc 
clear all; 
close all; 

%% Folders 
folderpath_10sem = "C:/Users/BuusA/OneDrive - Aalborg Universitet/10. semester (Kandidat)/Matlab files/Indiv_data/"; 

load(folderpath_10sem + "Thomas.mat");  Thomas = data_indiv;                TotalData{:,:,1} = data_indiv; 
load(folderpath_10sem + "Thomas_Type.mat");  Thomas_type = type;            TotalType{:,:,1} = type; 
load(folderpath_10sem + "Andreas.mat"); Andreas = data_indiv;               TotalData{:,:,2} = data_indiv; 
load(folderpath_10sem + "Andreas_Type.mat");  Andreas_type = type;          TotalType{:,:,2} = type; 
load(folderpath_10sem + "Maria.mat"); Maria = data_indiv;                   TotalData{:,:,3} = data_indiv; 
load(folderpath_10sem + "Maria_Type.mat");  Maria_type = type;              TotalType{:,:,3} = type; 
load(folderpath_10sem + "Benedikte.mat"); Benedikte = data_indiv;           TotalData{:,:,4} = data_indiv; 
load(folderpath_10sem + "Benedikte_Type.mat");  Benedikte_type = type;      TotalType{:,:,4} = type; 
load(folderpath_10sem + "Trine2.mat"); Trine2 = data_indiv;                 TotalData{:,:,5} = data_indiv; 
load(folderpath_10sem + "Trine2_Type.mat");  Trine2_type = type;            TotalType{:,:,5} = type; 
load(folderpath_10sem + "Gritt.mat"); Gritt = data_indiv;                   TotalData{:,:,6} = data_indiv; 
load(folderpath_10sem + "Gritt_Type.mat");  Gritt_type = type;              TotalType{:,:,6} = type; 
load(folderpath_10sem + "Mia.mat"); Mia = data_indiv;                       TotalData{:,:,7} = data_indiv; 
load(folderpath_10sem + "Mia_Type.mat");  Mia_type = type;                  TotalType{:,:,7} = type; 
load(folderpath_10sem + "Trine.mat"); Trine = data_indiv;                   TotalData{:,:,8} = data_indiv; 
load(folderpath_10sem + "Trine_Type.mat");  Trine_type = type;              TotalType{:,:,8} = type; 
load(folderpath_10sem + "Andrew.mat"); Andrew = data_indiv;                 TotalData{:,:,9} = data_indiv; 
load(folderpath_10sem + "Andrew_Type.mat");  Andrew_type = type;            TotalType{:,:,9} = type; 

subject{1,1} = "Thomas";    subject{1,2} = "Thomas_type"; 
subject{2,1} = "Andreas";   subject{2,2} = "Andreas_type"; 
subject{3,1} = "Beneditke"; subject{3,2} = "Beneditke_type"; %
subject{4,1} = "Maria";     subject{4,2} = "Maria_type"; 
subject{5,1} = "Trine2";    subject{5,2} = "Trine2_type"; 
subject{6,1} = "Griit";     subject{6,2} = "Gritt_Type.mat";
subject{7,1} = "Mia";       subject{7,2} = "Mia_Type.mat"; % probabaly exclude
subject{8,1} = "Trine";     subject{8,2} = "Trine_Type.mat"; % probala
subject{9,1} = "Andrew";     subject{9,2} = "Andrew_Type.mat"; % probala

%Andrew group nr 2. 

NSub = size(subject,1); 

%type{1} = VER_yes; type{2} = VER_no; type{3} = HOR_yes; type{4} = HOR_no; 
clear data_indiv type folderpath_10sem
%% Abbreviation
    CTL = 1; VER = 2; HOR = 3; CTL2 = 4; time = 5;
    SOL = 1; TA = 2; ANG = 3; FSR = 4; FOO = 5; 
    VEL = 6; ACC = 7; test = 8; SOL_diff = 9; SOL_diff_no_delay = 10; 
    ProtoAll = [CTL, VER, HOR, CTL2];

%% Functions
    msToSec = @(x) x*10^-3; % ms to sec 
    secToMs = @(x) x*10^3; % ms to sec 

%% Acquisition Set-Up
    sweep_length = 10;              % Signal length in second
    Fs = 2000;                      % Samples per second
    dt = 1/Fs;                      % Seconds per sample
    pre_trig = 4;                   % Pre-trigger 
    N = Fs*sweep_length;            % Total number of samples per signal

%% Speed and aceleration
span = 10; 
for sub = 1:NSub 
    for PROTO = [VER]
        data = TotalData{1,1,sub}; 
        type = TotalType{1,1,sub};
        for i = 1:size(data{PROTO,1},1) %sweeps

            %data{PROTO,ANG}(i,:) = data{PROTO,ANG}(i,:)*3.4114+18.9260; 
            data{PROTO,ANG}(i,:) = [rescale(data{PROTO,ANG}(i,1:3801)), rescale(data{PROTO,ANG}(i,3802:end))]; 
            data{PROTO,VEL}(i,:) = diff(data{PROTO,ANG}(i,:))/(dt*10^3);     % [deg/ms]
            data{PROTO,ACC}(i,:) = diff(smooth(data{PROTO,VEL}(i,:), span, 'moving'))/(dt*10^3);     % [deg/ms^2]
        end
        TotalData{1,1,sub} = data; 
        TotalType{1,1,sub} = type;
    end 
end

%% Plot
SLR_x = [39 59 59 39];
MLR_x = [60 80 80 60];
LLR_x = [81 100 100 81];

proto = VER; 
x_range = [-1000 800]; 
subjecsInc = [1,2,3,4,5,6,7,8,9];% 
%subjecsInc = [1,2];
%subjecsInc = 1;

if true    
    figure;
    if proto == HOR;  no = type{4}; yes = type{3}; sgtitle("Horizontal Pertubation"); end
    if proto == VER;  no = type{2}; yes = type{1}; sgtitle("Vertical Pertubation"); end
    % set(gcf,'Position',[100 100 1400 500])
    
    subplot(511);hold on; ylabel(""), xlim(x_range)%, ylim([-5 10])
        y = [0 0 1 1];
        patch(SLR_x,y,[0.75, 0.75, 0.75], 'EdgeColor', 'none','FaceAlpha',0.3)
        patch(MLR_x,y,[0.60, 0.60, 0.60], 'EdgeColor', 'none','FaceAlpha',0.3)
        patch(LLR_x,y,[0.75, 0.75, 0.75], 'EdgeColor', 'none','FaceAlpha',0.3)
        for sub = subjecsInc
            data = TotalData{1,1,sub}; 
            type = TotalType{1,1,sub};
            yes_mean = mean(data{proto,SOL}(type{1},:),1); 
            no_mean  = mean(data{proto,SOL}(type{2},:),1); 
            
            normalized_soleus_perturbation{sub,:} = data{proto,SOL}(type{1},:); 
            normalized_angle_perturbation{sub,:}  = data{proto,ANG}(type{1},:); 
            normalized_velocity_perturbation{sub,:}  = data{proto,VEL}(type{1},:); 
            normalized_acceleration_perturbation{sub,:}  = data{proto,ACC}(type{1},:); 


            time_axis = data{proto,time};

            plot(data{proto,time}, no_mean, "color", [0.7, 0.7, 0.7], "LineWidth",1)
            plot(data{proto,time}, yes_mean, "color", "black")
            plot(data{proto,time}, rescale((yes_mean)./no_mean), "color", "red")
            ylabel("Normalized"+newline+"Soleus")
            legend(["","","","control", "Vertical", "SOL (% back)"])           
        end 

    subplot(512);hold on; ylabel(""),title("Cross-correlation "); xlim(x_range)%, ylim([-5 10])  
        for sub = subjecsInc
            data = TotalData{1,1,sub}; 
            type = TotalType{1,1,sub}; VER_yes = type{1}; VER_no = type{2}; 

            y_yes =  mean(data{proto,SOL}(VER_yes,:),1); 
            y_no =  mean(data{proto,SOL}(VER_no,:),1); 

            y_yesTA =  mean(data{proto,TA}(VER_yes,:),1); 
            y_noTA =  mean(data{proto,TA}(VER_no,:),1); 

            y_ang_no =  mean(data{proto,ANG}(VER_no,:),1); 
            template  =  mean(data{proto,SOL}(VER_no,400:1500),1); 
            [rx, lags] = xcorr(y_yes, template);

            % Cross-correlation peaks 
            N = numel(y_yes); 
            plot(0:N-1, y_yes, 'color', 'black')
            plot(lags, rescale(rx), 'color', 'blue');

            % Find Peaks
            [pks, locs] = findpeaks(rescale(rx), 'MinPeakDistance', 500);
            tmp = find(lags(locs) > 0);
            pks = pks(tmp);
            locs = locs(tmp);

            % New signal with delay
            delay = 400-lags(locs(1));
            y_no = [y_no(delay+1:end), zeros(1, delay)]; 
            y_noTA = [y_noTA(delay+1:end), zeros(1, delay)]; 

            y_ang_no = [y_ang_no(delay+1:end), zeros(1, delay)]; 

            signal_diff = ((y_yes-y_no)./y_no); 
            TF = isinf(signal_diff); 
            signal_diff(TF) = 0; 
    
            % Find Difference for each YES sweeps
            for i = 1:numel(VER_yes)
                signal_diff_s = (data{proto,SOL}(VER_yes(i),:)-y_no)./y_no; 
                TF = isinf(signal_diff_s); 
                signal_diff_s(TF) = 0; 
                data{PROTO,SOL_diff}(i,:) = signal_diff_s; 
            end 
            TotalData{1,1,sub} = data; % Save data 

            % Cross-correlation peaks 
            plot(lags(locs), pks, 'rx');
            plot(lags(locs(1)), pks(1), 'o');
            plot(lags(locs(1)):lags(locs(1))+numel(template)-1, template, 'color', [0.7, 0.7, 0.7], 'LineWidth', 1);

        subplot(513); hold on; xlim(x_range)
            plot(data{proto,time}, y_no, 'color', [0.7, 0.7, 0.7], 'LineWidth', 2)
            plot(data{proto,time}, y_yes,  'color', 'black') 
            plot(data{proto,time}, rescale(signal_diff), 'color', 'red' )  
            ylabel("Normalized"+newline+"soleus"+newline+"Re-aligned")
            %plot(data{proto,time}, rescale((data{PROTO,SOL_diff}(1,:))), "blue")


            zero_idx = 2001; 
            SLR = [39,59]; % ms
            MLR = [60 80]; % ms
            SLR_idx = [SLR(1)*Fs*10^-3+zero_idx, SLR(2)*Fs*10^-3+zero_idx]; 
            MLR_idx = [MLR(1)*Fs*10^-3+zero_idx, MLR(2)*Fs*10^-3+zero_idx];
    

            kage(sub,1) = mean(y_yes(SLR_idx(1):SLR_idx(2)), "all");
            kage(sub,2) = mean(y_no(SLR_idx(1):SLR_idx(2)), "all");
            kage(sub,3) = mean(y_yes(MLR_idx(1):MLR_idx(2)), "all");
            kage(sub,4) = mean(y_no(MLR_idx(1):MLR_idx(2)), "all");
            kage(sub,5) = mean(y_yes(MLR_idx(1):MLR_idx(2)), "all");
            kage(sub,6) = mean(y_no(MLR_idx(1):MLR_idx(2)), "all");

            maage(sub,1) = mean(y_yesTA(SLR_idx(1):SLR_idx(2)), "all");
            maage(sub,2) = mean(y_noTA(SLR_idx(1):SLR_idx(2)), "all");
            maage(sub,3) = mean(y_yesTA(MLR_idx(1):MLR_idx(2)), "all");
            maage(sub,4) = mean(y_noTA(MLR_idx(1):MLR_idx(2)), "all");
            maage(sub,5) = mean(y_yesTA(MLR_idx(1):MLR_idx(2)), "all");
            maage(sub,6) = mean(y_noTA(MLR_idx(1):MLR_idx(2)), "all");

        
        subplot(514); hold on; xlim(x_range)
            plot(data{proto,time}, y_noTA, 'color', [0.7, 0.7, 0.7], 'LineWidth', 2)
            plot(data{proto,time}, y_yesTA,  'color', 'black') 

        subplot(515); hold on; xlim(x_range)

           plot(data{proto,time}, y_ang_no, 'color', [0.7, 0.7, 0.7], 'LineWidth', 2)
           plot(data{proto,time}, mean(data{proto,ANG}(VER_yes,:),1),  'color', 'black') 
           ylabel("Position"+newline+"[Deg]")

            
        end
end



%% Linear Regression
% Hvad er soleus aktivitet indenfor 38 ms efter onset til 70 ms efter onset

viewBefore = 1000; % [ms]
viewAfter  = 1000; % [ms]
zero_idx = find(data{proto,time} == 0);     %[sample]

% SOLEUS diff : 44,-5 and 55,9

%25, 15

% Unit [ms] , OBS must align 0.5
delay =      44 %10:1:100;  %5:1:30 %12.5 ; %[10:5:50]; % X X X
offset =     -5% -40:1:30;   %27.5; % [-100:5:100]; %  Y Y Y 
depIdxDiff = 30; %[20, 30, 40];
preIdxDiff = 30; %[20, 30, 40];

for off = 1:numel(offset) 
for de = 1:numel(depIdxDiff)
for pr = 1:numel(preIdxDiff)
for del = 1:numel(delay)

testSubject = 1; 
plt           = 0; 
pltAllSubject = 1; 
surfplt       = 1; 

    % Predictor
preIdx1 = zero_idx + msToSec(offset(off))*Fs;
preIdx2 = preIdx1 + floor(Fs*msToSec(preIdxDiff(pr)));
preIdx = [preIdx1:preIdx2];

    % Dependent
depIdx1 = zero_idx + msToSec(offset(off))*Fs + msToSec(delay(del))*Fs;
depIdx2 = depIdx1 + floor(Fs*msToSec(depIdxDiff(de)));
depIdx = [depIdx1:depIdx2];

pre = []; dep = [];
pre1 = []; dep1 = []; 
for sub = subjecsInc
    data = TotalData{1,1,sub}; 
    type = TotalType{1,1,sub}; VER_yes = type{1}; VER_no = type{2}; 
    for PROTO = [VER]
        
        data_dep{SOL_diff} = []; 
        for i = 1:size(data{VER,SOL_diff},1)
            data_dep{SOL_diff}(i) = mean(data{VER,SOL_diff}(i,depIdx)); 
        end 

        data_dep{SOL} = [];
        data_pre{VEL} = [];
        for i = 1:size(data{PROTO,SOL},1) % Sweeps  
            % Dependent Soleus
            %data_dep{SOL}(i) = mean(maxk(data{PROTO,SOL}(i,depIdx),30)); 
            data_dep{SOL}(i) = mean(data{PROTO,SOL}(i,depIdx)); 

            % Predictors
            data_pre{VEL}(i) = mean(abs(data{PROTO,VEL}(i,preIdx))); 
            %data_pre{VEL}(i) = mean(abs(mink(data{PROTO,VEL}(i,preIdx),20))); 
        end
    end 
    predictor1{sub} = data_pre{VEL}(VER_yes) - mean(data_pre{VEL}(VER_yes));   % Velocity, removed mean 
    %predictor1{sub} = rescale(data_pre{VEL}(VER_yes));                        % Velocity, rescaled
  
    dependent1{sub} = data_dep{SOL_diff}(:) - mean(data_dep{SOL_diff});       % Soleus Diff, removed mean
    %dependent1{sub} = rescale(data_dep{SOL_diff}(:));                          % Soleus Diff, rescaled
    %dependent1{sub} = data_dep{SOL}(VER_yes) - mean(data_dep{SOL}(VER_yes)); 
        pre1 = [pre1, predictor1{sub}];  % HERE
        dep1 = [dep1, dependent1{sub}']; % HERE 

    predictor{sub} = data_pre{VEL}(VER_yes)  - mean(data_pre{VEL}(VER_yes));     
    dependent{sub} = data_dep{SOL}(VER_yes)  - mean(data_dep{SOL}(VER_yes)); 
        pre = [pre, predictor{sub}];
        dep = [dep, dependent{sub}]; 
        
end 


% plot variables (SOL,ANG, VEL)
if plt 
    sub = testSubject; 
    data = TotalData{1,1,sub}; 
    type = TotalType{1,1,sub}; VER_yes = type{1}; VER_no = type{2}; 
    figure; 
    sweep = 2; % perturbation

    arrDep = depIdx; 
    arrPre = preIdx; 
     subplot(511); hold on; title("Mean SOL. "+ (depIdx1-zero_idx)*dt*1000 + "[ms] - "+ (depIdx2-zero_idx)*dt*1000+"[ms]") 
        plot(data{proto,time}, mean(data{VER,SOL}(VER_yes(:),:),1), "color", [0.7,0.7,0.7], "LineWidth",1)
        plot(data{proto,time}(1,arrDep), mean(data{VER,SOL}(VER_yes(:),arrDep),1), "color", "black", "LineWidth",2); 
    subplot(512); hold on; title("Sweep SOL")
        plot(data{proto,time}, data{VER,SOL}(VER_yes(sweep),:), "color", [0.7,0.7,0.7], "LineWidth",2)
        plot(data{proto,time}(1,arrDep), data{VER,SOL}(VER_yes(sweep),arrDep), "color", "black"); 
    subplot(513); hold on; title("SOL diff")
        plot(data{proto,time},mean(data{VER,SOL_diff}(:,:),1), "color", [0.7,0.7,0.7], "LineWidth",2); 
        plot(data{proto,time}(arrDep), mean(data{VER,SOL_diff}(:,arrDep),1), "color", "black"); 

    subplot(514); hold on; title("mean Velocity. "+ (preIdx1-zero_idx)*dt*1000 + "[ms]. - "+ (preIdx2-zero_idx)*dt*1000+"[ms]")
        plot(data{proto,time}(1,2:end), mean(data{VER,VEL}(VER_yes(:),:),1), "color", [0.7,0.7,0.7], "LineWidth",1); 
        plot(data{proto,time}(arrPre), mean(data{VER,VEL}(VER_yes(:),arrPre),1), "color", "black","LineWidth",2); 
    subplot(515); hold on; xlim([0 4001]); title("mean ankle")
       % plot(mean(data{VER,ANG}(VER_no,:),1), "color", [0.5,0.5,0.5], "LineWidth",2); 
        plot(mean(data{VER,ANG}(VER_yes,:),1), "color", [0.7,0.7,0.7], "LineWidth",1); 
        plot(arrPre, mean(data{VER,ANG}(VER_yes,arrPre),1), "color", "black","LineWidth",2); 
end

sub = testSubject; 
%mdl = fitlm(predictor{sub},dependent{sub}); 
if plt
   % mdl = fitlm(predictor{sub},dependent{sub}) 
   % mdl = fitlm(rescale(predictor{sub}),rescale(dependent{sub}))
   mdl = fitlm(predictor1{sub}, dependent1{sub})
     figure; plot(mdl); ylabel("Soleus"); xlabel("predictor"); title("subject " + sub)
end

% Linear regression, all subject 
if pltAllSubject 
    disp("prediktor1")
    mdl = fitlm(pre1,dep1) 
    b = table2array(mdl.Coefficients(1,1)); 
    a = table2array(mdl.Coefficients(2,1)); 
    linearReg = @(x) x*a + b; 
    figure; hold on; title("pre1")
    for i = subjecsInc
        plot(predictor1{i}, dependent1{i},'o'), ylabel(["Soleus (%background)","Window: 40-69 ms"]); xlabel(["Velocity","Window: -15-15 ms"])
    end
    plot(pre1, linearReg(pre1)); ylabel("Soleus /(%background)"); xlabel("Predictor")
end

% if pltAllSubject 
%     disp("prediktor")
%     mdl = fitlm(pre,dep) 
%     b = table2array(mdl.Coefficients(1,1)); 
%     a = table2array(mdl.Coefficients(2,1)); 
%     linearReg = @(x) x*a + b; 
%     figure; hold on; title("pre")
%     for i = subjecsInc
%         plot(predictor{i}, dependent{i},'o')
%     end
%     plot(pre, linearReg(pre)); ylabel("Soleus"); xlabel("Predictor")
% end
% 
% 
% mdl = fitlm(pre1, dep1);
% resultSig(off,del) = table2array(mdl.Coefficients(2,4)); 
% resultSize(off,del)= table2array(mdl.Coefficients(2,1));
% resultRsquare(off,del)= mdl.Rsquared.Adjusted; 

end; end; end; end

%% SURF PLOT
if surfplt 

    [X,Y] = meshgrid(delay, offset);
    Z = delay + offset'
    [x,y] = find(Z > 36 & Z < 44)
    k = zeros(size(Z))
    for i = 1:numel(x)
        k(x(i),y(i)) = 1;
    end
    
    Zpower = resultSize;
    Zsig = resultSig*(-1); 
    ZRsquare = resultRsquare; 
    Zpos = k; 
    
    figure
        subplot(221); surf(X,Y,Zpower); xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - power ")
        subplot(222); surf(X,Y,Zsig);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - sig "); zlim([-0.05 0])
        subplot(223); surf(X,Y,ZRsquare);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - Rsquare ")
        subplot(224); surf(X,Y,Zpos);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - pos ")
    
    [X,Y] = meshgrid(delay, offset);
    Z = delay + offset'
    [x,y] = find(or(Z < 35, Z > 40))
    for i = 1:numel(x)
        Zpower(x(i),y(i)) = 0;
        Zsig(x(i),y(i)) = -1;
        ZRsquare(x(i),y(i)) = 0;
    end
    
    figure 
        subplot(221); surf(X,Y,Zpower); xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - power "); title("Coefficients x1")
        subplot(222); surf(X,Y,Zsig);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - pValue "); zlim([-0.05 0]); title("p Value")
        subplot(223); surf(X,Y,ZRsquare);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - R^2 "); title("Adjusted R-Squared:")
        subplot(224); surf(X,Y,Zpos);   xlabel("x - delay [ms]"); ylabel("y - offset [ms]"); zlabel("z - pos "); title("Variables")
end

%% Plot differences: 
figure; 
xrange = [-150 150];
    for sub = [1,2,3,4,5,6,8,9]
        data = TotalData{1,1,sub}; 
        type = TotalType{1,1,sub}; VER_yes = type{1}; VER_no = type{2}; 
test = 9
subplot(311); hold on; ylabel("Normalized"+newline+"average SOL"); xlim(xrange)
    if sub == test 
        plot(data{VER,time}, mean(data{VER,SOL}(VER_yes,:),1), "color", "red", "LineWidth",1)
    else 
        plot(data{VER,time}, mean(data{VER,SOL}(VER_yes,:),1), "color", "black", "LineWidth",1)
    end 

subplot(312); hold on; ylabel("Normalized"+newline+"average"+newline+"SOL (%background)");xlim(xrange)
    if sub == test
        plot(data{VER,time}, mean(data{VER,SOL_diff}(:,:),1), "color", "red", "LineWidth",1); 
    else 
        plot(data{VER,time}, mean(data{VER,SOL_diff}(:,:),1), "color", "black", "LineWidth",1); 
    end 

subplot(313); hold on; ylabel("Normalized"+newline+"ankle"); xlim(xrange)
    if sub == test
        plot(data{VER,time}, mean(data{VER,ANG}(VER_yes,:),1)-mean(data{VER,ANG}(VER_yes,:),1:2), "color", "red", "LineWidth",1); 
    else 
        plot(data{VER,time}, mean(data{VER,ANG}(VER_yes,:),1)-mean(data{VER,ANG}(VER_yes,:),1:2), "color", "black", "LineWidth",1); 
    end 
end






%% DUMMY

%% BOKS plot 

zero_idx = 2001; 
SLR = [39,59]; % ms
MLR = [60 80]; % ms

offsets = [17.5,16.5,31.5,32.5,35.5,25.0,25.0,25.0]; %ORIGINAL
subjecsInc = [1,2,3,4,5,8,9]; %ORIGINAL!

%offsets = [17.5, 31.5, 32.5, 35.5]; POWERPUFF
%subjecsInc = [1,3,4,5]; <-- POWERPUFF!

% offsets =    [21,35,20,30,20,29];
% subjecsInc = [ 1, 2, 3, 4, 5,8];

% for i = 1:numel(subjecsInc)
%     data = TotalData{1,1,subjecsInc(i)}; 
%     type = TotalType{1,1,subjecsInc(i)}; no = type{4}; yes = type{3};
% 
%     zero_aligned_idx = zero_idx + (offsets(i)*Fs*10^-3); 
%     SLR_idx = [SLR(1)*Fs*10^-3+zero_aligned_idx, SLR(2)*Fs*10^-3+zero_aligned_idx]; 
%     MLR_idx = [MLR(1)*Fs*10^-3+zero_aligned_idx, MLR(2)*Fs*10^-3+zero_aligned_idx];
%     
%     
%     kage(i,1) = mean(data{HOR,SOL}(yes,SLR_idx(1):SLR_idx(2)), "all");
%     kage(i,2) = mean(data{HOR,SOL}(no, SLR_idx(1):SLR_idx(2)), "all");
%     
%     kage(i,3) = mean(data{HOR,SOL}(yes,MLR_idx(1):MLR_idx(2)), "all");
%     kage(i,4) = mean(data{HOR,SOL}(no, MLR_idx(1):MLR_idx(2)), "all");
% 
%     kage(i,5) = mean(data{HOR,SOL}(yes,SLR_idx(1):MLR_idx(2)), "all");
%     kage(i,6) = mean(data{HOR,SOL}(no, SLR_idx(1):MLR_idx(2)), "all");
% 
%     maage(i,1) = mean(data{HOR,TA}(yes,SLR_idx(1):SLR_idx(2)), "all");
%     maage(i,2) = mean(data{HOR,TA}(no, SLR_idx(1):SLR_idx(2)), "all");
%     
%     maage(i,3) = mean(data{HOR,TA}(yes,MLR_idx(1):MLR_idx(2)), "all");
%     maage(i,4) = mean(data{HOR,TA}(no, MLR_idx(1):MLR_idx(2)), "all");
% 
%     maage(i,5) = mean(data{HOR,TA}(yes,SLR_idx(1):MLR_idx(2)), "all");
%     maage(i,6) = mean(data{HOR,TA}(no, SLR_idx(1):MLR_idx(2)), "all");
% 
% end

%[p,h] = signrank(kage(:,1),kage(:,2)) % not significant. SLR SOL
% [p,h] = signrank(kage(:,3),kage(:,4)) % not significant, MLR, SOL
%[p,h] = signrank(kage(:,5),kage(:,6)) % not significant, SLR+MLR SOL

%[p,h] = signrank(maage(:,1),maage(:,2)) % not significant. SLR SOL
%[p,h] = signrank(maage(:,3),maage(:,4)) % not significant, MLR, SOL
%[p,h] = signrank(maage(:,5),maage(:,6)) % not significant, SLR+MLR SOL


%%
HOR8 = ["Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation";"Vertical perturbation"]; 
CTL8 = ["Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control";"Vertical Control"]; 
base8 = ["baseline";"baseline";"baseline";"baseline";"baseline";"baseline";"baseline";"baseline";"baseline"];

SLR8 = ["SLR";"SLR";"SLR";"SLR";"SLR";"SLR";"SLR";"SLR";"SLR"];
MLR8 = ["MLR";"MLR";"MLR";"MLR";"MLR";"MLR";"MLR";"MLR";"MLR"]; 
SLRMLR8 = ["SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR";"SLR+MLR"];

HOR4=HOR8(1:4); 
CTL4=CTL8(1:4);
base4=base8(1:4);
SLR4 = SLR8(1:4);
MLR4= MLR8(1:4);
SLRMLR4 = SLRMLR8(1:4);

HOR7 = HOR8(1:7); 
CTL7 = CTL8(1:7);
base7= base8(1:7);
SLR7 = SLR8(1:7);
MLR7 = MLR8(1:7);
SLRMLR7 = SLRMLR8(1:7);


HOR6 = HOR8(1:6); 
CTL6 = CTL8(1:6);
base6= base8(1:6);
SLR6 = SLR8(1:6);
MLR6 = MLR8(1:6);
SLRMLR6 = SLRMLR8(1:6);


% type = [SLR4; SLR4; MLR4; MLR4; SLRMLR4; SLRMLR4];
% protocol = [HOR4; CTL4;HOR4; CTL4;HOR4; CTL4];

type = [SLR8; SLR8; MLR8; MLR8; SLRMLR8; SLRMLR8];
protocol = [HOR8; CTL8; HOR8; CTL8; HOR8; CTL8];

% type = [SLR7; SLR7; MLR7; MLR7; SLRMLR7; SLRMLR7];
% protocol = [HOR7; CTL7; HOR7; CTL7; HOR7; CTL7];

% type = [SLR6; SLR6; MLR6; MLR6; SLRMLR6; SLRMLR6];
% protocol = [HOR6; CTL6; HOR6; CTL6; HOR6; CTL6];

averagesSOL = [kage(:,1); kage(:,2);kage(:,3);kage(:,4);kage(:,5);kage(:,6)];
averagesTA = [maage(:,1); maage(:,2);maage(:,3);maage(:,4);maage(:,5);maage(:,6)];

kagen = table(type,protocol,averagesSOL,averagesTA);

monthOrder = {'SLR','MLR','SLR+MLR'};
kagen.type = categorical(kagen.type,monthOrder);


figure; sgtitle("All subject (N=9)")

subplot(211)
boxchart(kagen.type, kagen.averagesSOL, 'GroupByColor', kagen.protocol)
ylabel("Average normalized"+newline+"soleus activity")
legend

subplot(212)
boxchart(kagen.type, kagen.averagesTA, 'GroupByColor', kagen.protocol)
ylabel("Average normalized"+newline+"tibialis activity")

