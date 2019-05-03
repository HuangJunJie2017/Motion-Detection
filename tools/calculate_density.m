function density=calculate_density(point,similarity_para)
%% analysis samples using fast cluster technology.
%% Example: density=fast_cluster_analysis(point,similarity_para)
%% Input: point             - samples ,number*channel mat
%%        similarity_para   - parameter used to define the similarity between samples' feature
%% Output: density
%%
% clear all
% close all
% point=rand(200,2);
% similarity_para=0.5;
%%
point=double(point);
[N,c]=size(point);
distance=zeros(N,N);
for i=1:c
    feature_y=repmat(point(:,i),1,N);
    feature_x=repmat(point(:,i)',N,1);
    distance=distance+(feature_y-feature_x).^2;
end
dense=exp(-similarity_para*distance);
density=sum(dense,2);
