function [v,V, unique_operations] = EditDistance_all_cpp(string1,string2)
% Given two strings s1 and s2 (DNA or protein), the edit distance between 
% them is the minimum number of operations required to convert s1 to s2. 
% Allowed operations are:
%   Replacing one character of string by another character.
%   Deleting a character from string
%   Adding a character to string
% Cost of each type of operation are all the same, 1


% compute edit distance matrix
m=length(string1);
n=length(string2);
% Call cpp to compute editdist
[v,V] = editDist_only(string1,string2,m,n); 

% Call cpp to backtrace
opts = unique(backtrace(V,string1,string2));

% Translate backtrace from cpp output
unique_operations = cell(v,1);
m = size(V,1);
n = size(V,2);
for k = 1:length(opts)
    ij = double(opts(k)); 
    i = mod(mod(ij,m*n),m);
    j = floor(mod(ij,m*n)/m);
    operation = floor(ij/(m*n));
    %fprintf('%d\t%d\t%d\n',i,j,operation);
    if operation == 1
        unique_operations{k} = sprintf('mutate position %d to %s',i, string2(j));
    elseif operation == 2
        unique_operations{k} = sprintf('delete position %d',i);
    elseif operation == 3
        unique_operations{k} = sprintf('insert %s after position %d',string2(j),i);
    end        
end
unique_operations = unique(unique_operations);