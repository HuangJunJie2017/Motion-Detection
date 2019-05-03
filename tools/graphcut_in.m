function [graph_cut,center]=graphcut_in(num_vertices,edge,c,min_size,min_output_number,max_edge_w,min_size_output)
%
% close all 
% clear all
% 
% w=576;
% h=720;
% num_vertices=w*h;
% data=load('data.mat');
% edge=data.edge_rgb;
% c=5000;
% min_size=500;
% min_output_number=2;
% max_edge_w=200;
%%
[~,ord]=sort(edge.w,'ascend');
edge.w=edge.w(ord,:);
edge.end=edge.end(ord,:);
graph_cut=graphcut_setup(num_vertices);
threshold=c*ones(num_vertices,1);
for i=1:edge.num
%%     a=graphcut_find(graph_cut,edge.end(i,1));
    if(graph_cut.num<=min_output_number)
        break;
    end
    
    x=edge.end(i,1);
    y=x;
    while(y~=graph_cut.elt(y).posi)
        y= graph_cut.elt(y).posi;
    end
    graph_cut.elt(x).posi=y;
    a=y;
%%     b=graphcut_find(graph_cut,edge.end(i,2));    
    x=edge.end(i,2);
    y=x;
    while(y~=graph_cut.elt(y).posi)
        y= graph_cut.elt(y).posi;
    end
    graph_cut.elt(x).posi=y;
    b=y;
    
    if(a~=b)
        if((edge.w(i,1)<=threshold(a))&&(edge.w(i,1)<=threshold(b)))
%%             graph_cut=graphcut_joint(graph_cut,a,b);
            if(graph_cut.elt(a).rank>graph_cut.elt(b).rank)
                graph_cut.elt(b).posi=a;
                graph_cut.elt(a).size=graph_cut.elt(a).size+graph_cut.elt(b).size;
            else
                graph_cut.elt(a).posi=b;
                graph_cut.elt(b).size=graph_cut.elt(b).size+graph_cut.elt(a).size;
                if(graph_cut.elt(a).rank==graph_cut.elt(a).rank)
                    graph_cut.elt(b).rank=graph_cut.elt(b).rank+1;
                end
            end
            graph_cut.num=graph_cut.num-1;
%%             a=graphcut_find(graph_cut,a);
            y=a;
            while(y~=graph_cut.elt(y).posi)
                y= graph_cut.elt(y).posi;
            end
            graph_cut.elt(a).posi=y;
            a=y;
            
            threshold(a)=edge.w(i,1)+c/graph_cut.elt(a).size;
        end
    end
end


for i=1:edge.num
%%     a=graphcut_find(graph_cut,edge.end(i,1));
    if(graph_cut.num<=min_output_number)
        break;
    end
    x=edge.end(i,1);
    y=x;
    while(y~=graph_cut.elt(y).posi)
        y= graph_cut.elt(y).posi;
    end
    graph_cut.elt(x).posi=y;
    a=y;
%%     b=graphcut_find(graph_cut,edge.end(i,2));    
    x=edge.end(i,2);
    y=x;
    while(y~=graph_cut.elt(y).posi)
        y= graph_cut.elt(y).posi;
    end
    graph_cut.elt(x).posi=y;
    b=y;
    
    if((a~=b)&&((graph_cut.elt(a).size<min_size)||(graph_cut.elt(b).size<min_size))&&(edge.w(i,1)<max_edge_w))
%%             graph_cut=graphcut_joint(graph_cut,a,b);
            if(graph_cut.elt(a).rank>graph_cut.elt(b).rank)
                graph_cut.elt(b).posi=a;
                graph_cut.elt(a).size=graph_cut.elt(a).size+graph_cut.elt(b).size;
            else
                graph_cut.elt(a).posi=b;
                graph_cut.elt(b).size=graph_cut.elt(b).size+graph_cut.elt(a).size;
                if(graph_cut.elt(a).rank==graph_cut.elt(a).rank)
                    graph_cut.elt(b).rank=graph_cut.elt(b).rank+1;
                end
            end
            graph_cut.num=graph_cut.num-1;
    end
end

center=zeros(graph_cut.num,2);
n=1;
%% 同一类标签归一化
for i=1:num_vertices
    y=i;
    while(y~=graph_cut.elt(y).posi)
        y= graph_cut.elt(y).posi;
    end
    graph_cut.elt(i).posi=y;
    if(size(find(center(:,1)==y),1)==0)
        center(n,1)=y;
        center(n,2)=graph_cut.elt(y).size;
        n=n+1;  
    end
end
[~,ord]=sort(center(:,2),'descend');
center=center(ord(1:min_output_number,:),:);
%%
if size(center,1)~=0
    min_size_output=center(1,2)/3;
end
%%
[idx_output,~]=find(center(:,2)>min_size_output);
center=center(idx_output,:);
%%
% color=zeros(num_vertices,3);
% for i=1:num_vertices
%     color(i,:)=floor(rand(1,3)*255);
% end
% im_out=zeros(w,h,3);
% for i=1:num_vertices
%     x=mod(i,w);
%     y=floor(i/w)+1;
%     if x==0
%         x=w;
%         y=y-1;
%     end
%     for j=1:3
%         im_out(x,y,j)=color(graph_cut.elt(i).posi,j);
%     end
% end
% figure
% imshow(uint8(im_out));