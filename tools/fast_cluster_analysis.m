function [center,output_number]=fast_cluster_analysis(point,similarity_para_flow,similarity_para_coor,min_desity,min_distanceflow,min_distancecoor,max_analyze_number)
%% analysis samples using fast cluster technology.
%% Example: density=fast_cluster_analysis(point,similarity_para)
%% Input: point             - samples ,number*channel mat
%%        similarity_para   - parameter used to define the similarity between samples' feature
%%        min_desity        - minimum desity required for selecting centers
%%        min_distance      - minimum distance required for selecting centers
%% Output: evolved image
%%

% clear all
% close all
% % point=rand(200,2);
% load test
% point=point_feature;
% similarity_para_flow=1;
% similarity_para_coor=1/30^2;
% min_desity=10;
% min_distanceflow=20;
% min_distancecoor=2000;
%%  显示前景采样点
% img_out=zeros(288,352);
% [N,c]=size(point);
% for i=1:N
%     img_out(point(i,3),point(i,4))=255;
% end
% figure(1)
% imshow(uint8(img_out));
%%
point=double(point);
[N,~]=size(point);
if N>max_analyze_number
    [idx,~]=find(rand(N,1)<(max_analyze_number/N));
    point=[point(idx,:),idx];
    min_desity=max_analyze_number/N*max_analyze_number/50;
else
    min_desity=max_analyze_number/50;
    point=[point,(1:N)'];
end

[N,~]=size(point);
distance_flow=zeros(N,N);
distance_coor=zeros(N,N);
for i=1:2
    feature_y=repmat(point(:,i),1,N);
    feature_x=repmat(point(:,i)',N,1);
    distance_flow=distance_flow+(feature_y-feature_x).^2;
end
for i=3:4
    feature_y=repmat(point(:,i),1,N);
    feature_x=repmat(point(:,i)',N,1);
    distance_coor=distance_coor+(feature_y-feature_x).^2;
end
distance=similarity_para_flow*distance_flow+similarity_para_coor*distance_coor;
% distance=similarity_para_coor*distance_coor;
dense=exp(-0.5*distance);
density=sum(dense,2);
%% 显示成分分析时各样本点的密度
ww=143;
hh=175;
maxd=max(density);
img_out=zeros(ww,hh);
for i=1:N
    img_out(point(i,3)/2,point(i,4)/2)=density(i)/maxd;%floor(density(i));%/maxd*255);
end
figure(21)
img_temp=ones(ww,hh);
[xx,yy]=find(img_temp);
xx=reshape(xx,ww,hh);
yy=reshape(yy,ww,hh);
surf(xx,yy,img_out);
colorbar
axis off
% ww=143;
% hh=175;
% img_out=zeros(ww,hh,3);
% img_out(:,:,3)=ones(ww,hh)*255*0.5;
% map=colormap;
% map=[0,0,0.5;map];
% maxdensity=max(max(density));
% for i=1:N
%     idx=floor(density(i)/maxdensity*64)+1;
%     img_out(point(i,3)/2,point(i,4)/2,:)=map(idx,:)*255;
% end
% figure(22)
% imshow(uint8(img_out));
% imwrite(uint8(img_out),'densitymap.png');
%%
maxdf=max(max(distance_flow));
maxdc=max(max(distance_coor));
[~,ord_density]=sort(density,'descend');
deltaflow=zeros(N,1);
deltacoor=zeros(N,1);
deltaflow(ord_density(1))=-1;
deltacoor(ord_density(1))=-1;
for i=2:N
    deltaflow(ord_density(i))=maxdf;
    deltacoor(ord_density(i))=maxdc;
    for j=1:i-1
        if(distance_flow(ord_density(i),ord_density(j))<deltaflow(ord_density(i)))
            deltaflow(ord_density(i))=distance_flow(ord_density(i),ord_density(j));
        end
        if(distance_coor(ord_density(i),ord_density(j))<deltacoor(ord_density(i)))
            deltacoor(ord_density(i))=distance_coor(ord_density(i),ord_density(j));
        end
    end
end

deltaflow(ord_density(1))=max(deltaflow(:));
deltacoor(ord_density(1))=max(deltacoor(:));

% 
% img_out=zeros(59,79);
% for i=1:N
%     img_out(point(i,3)/4,point(i,4)/4)=deltacoor(i);%floor(density(i));%/maxd*255);
% end
% figure(22)
% img_temp=ones(59,79);
% [xx,yy]=find(img_temp);
% xx=reshape(xx,59,79);
% yy=reshape(yy,59,79);
% surf(xx,yy,img_out);
% colorbar
% axis off


center=zeros(N,8);
output_number=0;
for i=1:N
    if density(i)>min_desity&&(deltaflow(i)>min_distanceflow||deltacoor(i)>min_distancecoor)
%     if density(i)>min_desity&&deltacoor(i)>min_distancecoor
        output_number=output_number+1;
        center(output_number,1)=point(i,5);
        center(output_number,2)=density(i);
        center(output_number,3)=deltaflow(i);
        center(output_number,4)=deltacoor(i);
        center(output_number,5:8)=point(i,1:4);
    end
end

if output_number==0
    if density(ord_density(1))>min_desity
        output_number=1;
        center(1,1)=point(ord_density(1),5);
        center(1,2)=density(ord_density(1));
        center(1,3)=deltaflow(ord_density(1));
        center(1,4)=deltacoor(ord_density(1));
        center(1,5:8)=point(ord_density(1),1:4);
    end
end
center=center(1:output_number,:);


