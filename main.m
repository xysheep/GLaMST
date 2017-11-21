function main(fastafilename, varargin)
pnames = {'rewire'};
dflts = {1};
[rewire] = internal.stats.parseArgs(pnames,dflts,varargin{:});
[~, s] = fastaread(fastafilename);
reconstruct_tree_minimun_tree_size(s, rewire);