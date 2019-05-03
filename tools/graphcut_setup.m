function graph_cut=graphcut_setup(num_vertices)
%% test para
% clear all
% close all
% num_vertices=50;
%% 
graph_cut.num=num_vertices;
for i=1:num_vertices
    graph_cut.elt(i).rank=0;
    graph_cut.elt(i).size=1;
    graph_cut.elt(i).posi=i;
    graph_cut.elt(i).possi=0;
end