%% Prepare Environment
% Set working directory to the upzipped folder. If your matlab does not have
% C/C++ compiler setup, need to first install MinGW, instructions at 
% following link (Home->Add-Ons->search MinGW->install) 
% <https://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html>

restoredefaultpath;
addpath lib
addpath src
% If your matlab does not have C/C++ copmiler setup, also run 
mex src/editDist_only.cpp
mex src/backtrace.cpp


%% Simulate data
% Simulate a dataset that 
operation_probability = [0.99 0.005 0.005];% mutation, insertion, deletion chance
sequence_length = 300; % Length of the seqeunces
num_tree_nodes = 200; % Total number of the real tree
sample_size = 100; % Number of observed nodes
rng(42); %% Fix the random seed;
[observed_sequences, true_sequences, adj, is_selected] = simulate_data_v2(...
    sequence_length,num_tree_nodes,sample_size,operation_probability);
draw_hierarchy_tree(adj, is_selected');

%% Run GLaMST and visualize the results
% This section is a example of using GLaMST. You need to pass a variable
% "observed_sequences" that contain a cell array of known sequences to
% function "reconstruct_tree_minimun_tree_size". ** Make sure the first
% sequence of "observed_sequences" should be root(G.L.). **
% "reconstructed_nodes" is a cell array of sequences of the reconstructed
% tree. "reconstructed_is_selected" is the index of tree nodes that known
% before the reconstruction. "reconstructed_directed_adj" is the adjecent
% matrix that represent the hierarchical structure of the tree.
[reconstructed_nodes,mst_adj,reconstructed_is_selected, reconstructed_directed_adj] = ...
    reconstruct_tree_minimun_tree_size(observed_sequences);
draw_hierarchy_tree(reconstructed_directed_adj, reconstructed_is_selected);


%% Generate lineage tree from real data
% Given this data is not yet publicly avaliable, we have simulated a
% dataset that has exactly the same topology as used in our paper. This
% simulated data is avaliable here <>. 
fsa = fastaread('demodata/real.fasta');
observed_sequences = fas.Sequence;
[reconstructed_nodes,mst_adj,reconstructed_is_selected, reconstructed_directed_adj] = ...
    reconstruct_tree_minimun_tree_size(observed_sequences);
draw_hierarchy_tree(reconstructed_directed_adj, reconstructed_is_selected);