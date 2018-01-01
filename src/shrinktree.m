function [adj] = shrinktree(adj)
queue = zeros(size(adj,1),1);
top = 1;
queue(top) = 1;
while top > 0
    head = queue(top);
    top = top - 1;
    current = head;
    nexts = find(adj(current,:));
    for nt = nexts
        next = nt;
        current = head;
        while length(next) == 1
            adj(current,next) = 0;
            adj(head,next) = adj(head,current) + 1;
            adj(head,current) = 0;
            current = next;
            next = find(adj(current,:));
        end
        queue(top+1) = current;
        top = top + 1;
    end
end
idx = intersect(find(sum(adj,1)==0), find(sum(adj,2)==0));
adj(idx,:) = [];
adj(:,idx) = [];