function [V,v] = EditDistance_only(string1,string2)
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
v=zeros(m+1,n+1);
edge_list = [];
for i=1:1:m
    v(i+1,1)=i;
end
for j=1:1:n
    v(1,j+1)=j;
end
for i=1:m
    for j=1:n
        if (string1(i) == string2(j))
            v(i+1,j+1)=v(i,j);
        else
            v(i+1,j+1)=1+min(min(v(i+1,j),v(i,j+1)),v(i,j));
        end
    end
end
V=v(m+1,n+1);

