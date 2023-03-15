function [EMG_norm] = func_normalize_EMG(step_index, data, varargin)

% MYFUNCTION - This function returns the EMG signal normalized to specified
% step
%
% EMG_norm = func_normalize_EMG(step_index, data, 'normalize_to_step')
% EMG_norm = func_normalize_EMG(step_index, data, 'span')
% EMG_norm = func_normalize_EMG(step_index, data, 'protocols')

    % abbreviation
    CTL = 1; VER = 2; HOR = 3; CTL2 = 4; % protocols 
    SOL = 1; TA = 2; ANG = 3; FSR = 4; time = 5; VEL = 6; ACC = 7; % sensor types

    % default values and validiation 
    defaultStep = 2;
    defaultSpan = 15; 
    expectedStepsObtions = [0,2,4,6];
    expectedProtocolObtions = [CTL, VER, HOR, CTL2]; 
    ecpectedProtocol = 2;
    
    validationSpanFcn = @(x) assert(numel(x) <= 1, 'Assertion failed: x contains more than 1 element.');
    validationStepFcn = @(x) assert(ismember(x, expectedStepsObtions), 'Assertion failed.');
    validationProtoFcn = @(x) assert(all(ismember(x, expectedProtocolObtions)), 'Error: Invalid input options.');

    % variable input argument list
    p = inputParser;
    addParameter(p, 'normalize_to_step',defaultStep, validationStepFcn); 
    addParameter(p, 'span', defaultSpan, validationSpanFcn)
    addParameter(p, 'protocols',expectedProtocolObtions, validationProtoFcn)
    parse(p,varargin{:})

    % find rise of FSR sensor to specified step
    [rise, fall] = func_find_edge(p.Results.normalize_to_step); 
   
    for proto = p.Results.protocols % loop through protocols
        clear SOL_max TA_max
        for i = 1:size(data{proto,1},1) % loop through sweeps      
            index_array = step_index{proto}((i),rise) + ceil((step_index{proto}((i),fall)-step_index{proto}((i),rise))/2)  : step_index{proto}((i),fall); 
            SOL_max(i) = max(smooth(data{proto,SOL}((i), index_array),p.Results.span, 'moving'));
            TA_max(i)  = max(smooth(data{proto,TA} ((i), index_array),p.Results.span, 'moving'));
        end
    
        % divide every element with found average
        data{proto,SOL} = data{proto,SOL} / mean(SOL_max);
        data{proto,TA} = data{proto,TA} / mean(TA_max);

    end 
    EMG_norm = data; 
end 
