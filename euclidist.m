function [distance, i] = euclidist(data, weight, som_row,...
                                som_col, m, n )
% This function finds the best matched code vector (Winning neuron) based
% on the input image
% 
    % Initialize matrix for storing the Euclidean distance between the input
    % vector and each neuron
    distance = zeros(som_row, som_col);

    i = randi([1 m]);
    
    for row = 1:som_row
        for col = 1:som_col
            sub = data(i,:) - reshape(weight(row,col,:),1,n);
            distance(row,col) = sqrt(sub * sub');
        end
    end

end
