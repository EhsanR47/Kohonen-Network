%Ehsan Rassekh

clear ;
close all;
clc;
 
file = fopen('iris.data');

%read the txt file using textscan function
textdata = textscan(file,'%f %f %f %f %s', 200, 'Delimiter',',');

%form the data matrix
data = cell2mat(textdata(:,1:4));
[m,n] = size(data);

%%%plot data ,only using first two features
figure(1);
plot(data(:,1),data(:,2),'k*','MarkerSize',5)
title ' Iris Data';
xlabel 'sepal length (cm)'; 
ylabel 'sepal Width  (cm)';

%################## Define the SOM Architecture ##################%
% Determine the number of rows and columns in the som
som_row = 20;
som_col = 20;

% Number of epoch
epoch = 500;
%##################initialization##################%
%%winner neuron neighbour
width_Initial = 8;

l_rate = 1 ;

a_max = 0.9 ;
a_min = 0.1 ;
t_max = epoch;
eta = a_max;

t_width = epoch/log(width_Initial);
%##################initialization weight matrix##################%
weight = zeros(som_row, som_col, n);

for i=1:som_row	
	for j=1:som_col
		weight(i, j, :) = rand(1,n);
	end
end

%########################################################################%
for iter=1:epoch
    
    
    width_var = (width_Initial * exp(-iter / t_width))^2;
	l_rate = (a_max - a_min)* ((t_max - iter)/(t_max - 1)) + a_min;

	[distance, index] = euclidist(data, weight, som_row,som_col, m, n );

	[minm,ind] = min(distance(:));
	[row_winner,col_winner] = ind2sub(size(distance),ind);

	dist = computedist( som_row, som_col, row_winner, col_winner, width_var);

	%##################update weight##################

	
    
    for row = 1: som_row
       for col = 1:som_col
           
           % Reshape the dimension of the current weight vector
           weight_vec = reshape(weight(row,col,:),1,n);
           
           % Update the weight vector for each neuron
           weight(row,col,:) = weight_vec + l_rate*dist(row,col)*(data(index,:)-weight_vec);
            
       end
    end
	
	
	%% ################## Illustrate The Updated Clustering Results  ##################
    % Weight vector of neuron
    dot = zeros(som_row*som_col, n);
    % Matrix for SOM plot grid
    matrix = zeros(som_row*som_col,1);
    % Matrix for SOM plot grid for deletion
    matrix_old = zeros(som_row*som_col,1);
    
   
    ind = 1;  
    hold on;
    f1 = figure(1);
    set(f1,'name',strcat('Iteration #',num2str(iter)),'numbertitle','off');

    % Retrieve the weight vector of neuron
    for r = 1:som_row
        for c = 1:som_col      
            dot(ind,:)=reshape(weight(r,c,:),1,n);
            ind = ind + 1;
        end
    end

    % Plot SOM
    for r = 1:som_row
        Row_1 = 1+som_row*(r-1);
        Row_2 = r*som_row;
        Col_1 = som_row*som_col;

        matrix(2*r-1,1) = plot(dot(Row_1:Row_2,1),dot(Row_1:Row_2,2),'--ro','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4);
        matrix(2*r,1) = plot(dot(r:som_col:Col_1,1),dot(r:som_col:Col_1,2),'--ro','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4);

        matrix_old(2*r-1,1) = matrix(2*r-1,1);
        matrix_old(2*r,1) = matrix(2*r,1);

    end

    % Delete the SOM plot from previous iteration
    if iter~=epoch  
        for r = 1:som_row
            delete(matrix_old(2*r-1,1));
            delete(matrix_old(2*r,1));
            drawnow;
        end
    end

    
    %% ################## The End Of Illustration  ##################
	

end

%%find winner neuron,No weight updates
top_map = zeros(som_row, som_col);

for i=1:m

		distance = zeros(som_row, som_col);

		
		
		for row = 1:som_row
			for col = 1:som_col
				sub = data(i,:) - reshape(weight(row,col,:),1,n);
				distance(row,col) = sqrt(sub * sub');
			end
		end
		
		
		
		
		[minm,ind] = min(distance(:));
		[row_winner,col_winner] = ind2sub(size(distance),ind);
		
        top_map(row_winner,col_winner) = i
		

end

display(top_map)
