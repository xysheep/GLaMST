function [v,V] = EditDistance_only_cpp(string1,string2)
% Given two strings s1 and s2 (DNA or protein), the edit distance between 
% them is the minimum number of operations required to convert s1 to s2. 
% Allowed operations are:
%   Replacing one character of string by another character.
%   Deleting a character from string
%   Adding a character to string
% Cost of each type of operation are all the same, 1

%unused variables
% compute edit distance matrix
m=length(string1);
n=length(string2);
[v,V] = editDist_only(string1,string2,m,n); 