function dist = computedist( somRow, somCol, row_winner, col_winner, width_var)
% This function compute the lateral distance between neurons i and the winner neurons 
   
    dist = zeros(somRow, somCol);
    
    for row = 1:somRow
       for col = 1:somCol
           if (row == row_winner) && (col == col_winner) % winner neuron
               
               dist(row,col) = 1;
           else
               
               distance = (row_winner - row)^2+(col_winner - col)^2;
               dist(row,col) = exp(-distance/(2*width_var)); %neighborhood function for other neurons
           end    
       end
    end
end
