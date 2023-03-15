function [step_index] = func_step_index_corr(varargin)  

%function [step_index] = func_step_index_corr(FSR, str_direction, edge_num, move_num, sweep)  

inputNames = {'FSR', 'direction', 'edge_num', 'move_num', 'sweep'};
p = inputParser;

for i = 1:length(inputNames)
    addOptional(p, inputNames{i}, []);
end
parse(p, varargin{:});

FSR = p.Results.FSR;
str_direction = p.Results.direction;
edge_num = p.Results.edge_num;
move_num = p.Results.move_num;
sweep = p.Results.sweep;


    switch str_direction 
        case "right"
            dir = 1; 
         case "left"
            dir = -1; 
        otherwise
            msg = "Not valid direction. Possible inputs: left or right";
            error(msg)
    end

    offset = 7950; 
    edges_index_total = find(edge(FSR(sweep,offset:end)))+offset; %  Index pos of change in FSR
    edges_index_removed = edges_index_total([1:edge_num-1, (edge_num+move_num):end]);
    step_index = flip(edges_index_removed(1:9));
  
end
