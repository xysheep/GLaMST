function d2leave = leavedist(adj)
    d2leave = zeros(size(adj,1),1);
    leaves = find(sum(adj,2)==0)';
    queue = leaves;
    while ~isempty(queue)
        top = queue(1);
        queue = queue(2:end);
        parent = find(adj(:,top));
        if ~isempty(parent) && ~ismember(parent, queue)
            queue = [queue parent];
        end
        children = find(adj(top,:));
        if isempty(children)
           d2leave(top) = 1;
        else
           d2leave(top) = max(1+d2leave(children));
        end
    end
end