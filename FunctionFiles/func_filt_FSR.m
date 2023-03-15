function [FSR] = func_filt_FSR(FSR, varargin)
    % func_filt_FSR - This function returns filtered FSR signal and removed
    % noise gaps. 
    %
    % FSR_filtered = func_filt_FSR(FSR)
    % 
    % func_filt_FSR(X,'gap_size2') changes the length of the gap noises
    % removel of filter 2 - default 300. 
    %
    % func_filt_FSR(X,'gap_size3') changes the length of the gap noises
    % removel of filter 3 - default 600. 
    %
    % func_filt_FSR(X,'limit_pct') changes the threshold to bypass 
    %
    % func_filt_FSR(X,'test') enable testing. 
    
    sweep_length = 10;              % Signal length in second
    Fs = 2000;                      % Samples per second
    dt = 1/Fs;                      % Seconds per sample
    N = Fs*sweep_length;            % Total number of samples per signal

    x_axis = linspace(-4, 6-dt, N); 
    Nsweep = size(FSR,1); 
    default_gap_size2 = 300; 
    default_gap_size3 = 600; 
    default_pct = 0.95; 
    default_test = 0; 

    p = inputParser;
    addParameter(p,'gap_size2', default_gap_size2)
    addParameter(p,'gap_size3', default_gap_size3)
    addParameter(p,'limit_pct', default_pct)
    addParameter(p,'test', default_test)

    parse(p,varargin{:}) 

    if p.Results.test == true
        fprintf("(1) Unique numbers before: " + numel(unique(FSR)) );
    end

    % Remove small fluxuations
    for n = 1:Nsweep
        for i = 1:size(FSR,2)
            if FSR(n,i) > 2.5
                FSR(n,i) = 5; 
            elseif FSR(n,i) < 2.5
                FSR(n,i) = 0; 
            end
        end
    end

    if p.Results.test == true    % Test code 
        fprintf(". Unique numbers after: " + numel(unique(FSR)) + "\n")
        figure(2); 
    end 


    %% (2) Remove noise spikes if gap is larger than 300 samples (0.15 sec) 

    Continue_loop = p.Results.test;
    correctInput = false; 
    prompt = "Continue, press >c<" + newline + "Quite, press >q<" + newline;   

    for sweep = 1:Nsweep % loop through sweeps
        edge_indexes = find(edge(FSR(sweep,:))); % Index of change in FSR
    
        for i = 1:numel(edge_indexes)-1 % loop though changes
            dur = edge_indexes(i+1) - edge_indexes(i); % find length between two changes 
            if dur < p.Results.gap_size2 % remove change if less than 
                index = edge_indexes(i);
                
                if FSR(sweep, index+2) < 2.5 % change if is FSR signal is positiv to the rigth side of change
                    if Continue_loop == true % display changes if enabled
                        sgtitle("(2) Sweep " + sweep + ". Sec: " + string(index*dt-4) )
                        subplot(211); 
                        title("before"); 
                        plot(x_axis, FSR(sweep,:)) 
                        ylim([0 6])

                    end 
                    FSR(sweep,edge_indexes(i)-1:edge_indexes(i+1)+1) = 5; % remove gap 
                    if Continue_loop == true
                        subplot(212); 
                        title("after");  
                        plot(x_axis, FSR(sweep,:)) 
                        ylim([0 6])

                        while correctInput == false
                            str = input(prompt, 's');
                            if strcmp(str,"q")
                                disp("Loop stopped")
                                Continue_loop = false; correctInput = true; 
                            elseif strcmp(str,"c")
                                correctInput = true; 
                            end 
                            if correctInput == false
                                warning("Input not accepted")
                            end
                        end                         
                        correctInput = false; 
                        clc
                    end
                end
            end
        end
    end 

    %% (3) Remove noise if gap deviate from from 75 pct of average

    avg_FSR = rescale(sum(FSR,1)); 

    if p.Results.test == true    % Test code 
        close 2
        figure(2)

        subplot(311); hold on 
        plot(x_axis, avg_FSR)
        plot([-4 6],[p.Results.limit_pct, p.Results.limit_pct])       
        hold off
    end 
  
    Continue_loop = p.Results.test; 
    correctInput = false; 
   
    for sweep = 1:Nsweep % loop through sweeps
        edge_indexes = find(edge(FSR(sweep,8000:end)))+8000; % Index of change in FSR    
        for i = 1:numel(edge_indexes)-1 % loop though changesc


            if FSR(sweep, edge_indexes(i)+10) < 2.5 % change if is FSR signal is 0 to the rigth side of index   

                dur = floor(numel(edge_indexes(i):edge_indexes(i+1))/2);
                limit = avg_FSR(edge_indexes(i) + dur) ;
                gap =  numel(edge_indexes(i):edge_indexes(i+1));
                if p.Results.limit_pct < limit && gap < p.Results.gap_size3
                    if Continue_loop == true % display changes if enabled
                        sgtitle("(3) Sweep " + sweep + ". Sec: " + string(edge_indexes(i)*dt-4) + newline + ". Limit: " + limit)
                        
                        subplot(312);
                        title("before");  
                        ylim([0 6])
                        plot(x_axis, FSR(sweep,:)) 
                        ylim([0 6])

                    end 
                    FSR(sweep,edge_indexes(i)-1:edge_indexes(i+1)+1) = 5; % remove gap 
                    if Continue_loop == true
                        subplot(313);
                        title("after"); 
                        plot(x_axis, FSR(sweep,:)) 
                        ylim([0 6])

                        while correctInput == false
                            str = input(prompt, 's');
                            if strcmp(str,"q")
                                disp("Loop stopped")
                                Continue_loop = false; correctInput = true; 
                            elseif strcmp(str,"c")
                                correctInput = true; 
                            end 
                            if correctInput == false
                                warning("Input not accepted")
                            end
                        end                         
                        correctInput = false; 
                        clc
                    end
                end
            end 
        end
    end 


%% (4) swipe through all 

    if p.Results.test == true    % Test code 
        close 2
        figure(2)
        subplot(211)
        plot(x_axis, avg_FSR)
        ylim([0 6])
        Continue_loop = p.Results.test; 
        correctInput = false; 

        for sweep = 1:Nsweep % loop through sweeps

            if Continue_loop == true % display changes if enabled
                sgtitle("(4) Sweep " + sweep)  
                subplot(212)
                plot(x_axis, FSR(sweep,:)) 
                ylim([0 6])
    
                while correctInput == false
                    str = input(prompt, 's');
                    if strcmp(str,"q")
                        disp("Loop stopped")
                        clear 2
                        Continue_loop = false; correctInput = true; 
                    elseif strcmp(str,"c")
                        correctInput = true; 
                    end 
                    if correctInput == false
                        warning("Input not accepted")
                    end
                end                         
                correctInput = false; 
                clc
            end
        end 
    end
end  