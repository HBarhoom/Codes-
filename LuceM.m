function [Bias,Sim] = LuceM(confMat)

% the Luce's choice model (Luce, 1963) is used to calculate response biases and 
% similarities from a given confusion matrix. 
%This function returns the response biases vector (Bias) and similarities
%matrix (Sim) for a given confusion matrix (confMat).  

% Hatem Barhoom 
% 25/05/2020

matrix = confMat; %% cofMat is the confusion matrix of the raw data (for example: presented Vs responsded)
matrix(matrix==0)=0.005;%% to avoid zero values in the matrix. 
n = length(matrix);
in_matrix = ones (n); %% the starting matrix. 


for k = 1:1000
    
    % iterative proportional fitting 
    
       %adjust for rows
    for c = 1:n;
        for r = 1:n;
            d = ((in_matrix(r,c))/(sum(in_matrix(r,1:n)))*(sum(matrix(r,1:n))));
            matrix_fit(r,c)= d;
        end
    end
    
    matrix1 = matrix_fit;
    
    %adjust for columns
    for c = 1:n;
        for r = 1:n;
            v = ((matrix1(r,c))/(sum(matrix1(1:n,c)))*(sum(matrix(1:n,c))));
            matrix_fit(r,c)= v;
        end
    end
    
   matrix2 = matrix_fit;
   
   %adjust for similarities:
   for c = 1:n;
        for r = 1:n;
            v = ((matrix2(r,c))/((matrix2(r,c)+matrix2(c,r)))*((matrix(r,c)+matrix(c,r))));
            matrix_fit(r,c)= v;
        end
   end
    
   
   in_matrix = matrix_fit; %% the starting matrix for teh next cycle. 
   
end

MLEM = in_matrix;    %maximum liklihood estimates


% calculate similarities (matrix)
for r = 1:n;
    for c = 1:n;
        
conf(r,c) = sqrt (((MLEM(r,c))*(MLEM(c,r))/((MLEM(r,r))*(MLEM(c,c)))));

    end 
end 
  
Sim = round(conf,3);

%calculate resposne biases (vector)
s = zeros(1,n);
for c = 1:n;
for  k = 1:n;         
s(c) = s(c) + ((sqrt((MLEM(c,k)*MLEM(k,k))/((MLEM(k,c)*(MLEM(c,c)))))));  
end 
end 
Bias = ((1./s)); 
Bias = round(Bias,3);
end 
