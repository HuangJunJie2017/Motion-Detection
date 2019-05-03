addpath tools
addpath(fullfile(db_matlab_root_dir,'db_util'));
addpath(fullfile(db_matlab_root_dir,'measures'));

clear all

para_sample_rate=10;
number_sample=4;
%iter_h=50;
threshold_RANSAC=5;
para_a=0.5;
para_b=1.8;
ival =2;
Ftype = 'FN2';
% Ftype = 'FN2-css-ft-sd';
threshold_co=0.05;
%map =[1,1;2,1;3,1;4,1;1,2;2,2;3,2;4,2;1,3;2,3;3,3;4,3;1,4;2,4;3,4;4,4];


w=480;
h=854;
x=1:w;
y=1:h;
xx=repmat(x',1,h);
yy=repmat(y,w,1);
point=[reshape(xx,1,w*h);reshape(yy,1,w*h)];
point_extend = [point.^2;point(1,:).*point(2,:);point;ones(1,w*h)]; 

%[w_d,h_d,~]=size(flow_d);
w_d=47;
h_d=85;

img_temp=ones(w_d,h_d);
%         [point_d(:,1),point_d(:,2)]=find(img_temp);
[point_d1,point_d2]=find(img_temp);
point_d = [point_d1,point_d2]';
point_d=point_d*para_sample_rate;
point_d_extend=[point_d.^2;point_d(1,:).*point_d(2,:);point_d;ones(1,w_d*h_d)]; 
% Get the ids of all sequences
seq_ids = db_seqs();
%%
%parameter
%for ival =2
for para_a=2.85
for para_b = 0.33
for number_sample=16
    pieces_x = 4;%ceil(sqrt(number_sample*4*w/h));
    pieces_y = 8;%ceil(sqrt(number_sample*4*h/w));
    num_pieces = pieces_x*pieces_y;
    map = constructmap(pieces_x,pieces_y);
for iter_h = 5
    select_idx_record=zeros(iter_h,number_sample);
    num_record=zeros(iter_h,1);
for times =1
for speed =25   
% Name of a result technique
%%
%times
% result_id = ['Ours_ai_pa028_pb020_t' sprintf('%03d',int32(times))];
%result_id = ['Ours_FN2cssftsd' sprintf('%03d',int32(para_a*100)) sprintf('_pb%03d',int32(para_b*100))];
% result_id = 'Ours_FN2cssftsd_test';
result_id = 'Ours_FN2_test-iterh5';
% result_id = 'Ours_ai_test';

result_dir = fullfile(db_root_dir,'Results','Segmentations',db_imsize, result_id);
mkdir(result_dir);

%%
for s_id = 1:length(seq_ids)
    seq_id = seq_ids{s_id};
    mkdir(fullfile(result_dir,seq_id));
    
    frame_ids = db_frame_ids(seq_ids{s_id});
    fprintf('%s contains %d images: \n',seq_ids{s_id},length(frame_ids));
    scene_speed_tmp =0;
    for f_id = 1:length(frame_ids)-1
        flopath = fullfile(db_root_dir,'FLO',Ftype,db_imsize,sprintf('ival%01d',ival),seq_id,sprintf('%05d.flo',f_id));
        %flopath = 'G:\data\flow_final\13\flow_1\0432.flo';
        flow=readFlowFile(flopath);
%         imgcl=flowToColor(flow);
%         figure(2)
%         imshow(imgcl);
        
        
        %下采样
        flow_d=img_desample(flow,para_sample_rate);
%         [w_d,h_d,~]=size(flow_d);
        flow_d=reshape(flow_d,w_d*h_d,2);
        num_record_max=0;
        flo_diff_max =10^10;
        point_selet_idx_max = 0;
        for k=1:iter_h
            %with CRA
            area_idx=randperm(num_pieces,number_sample);
            area_idx_x=map(area_idx,1);
            area_idx_y=map(area_idx,2);
            area_idx_x=floor((area_idx_x'-1+rand(1,number_sample))*w_d/pieces_x)+1;
            area_idx_y=floor((area_idx_y'-1+rand(1,number_sample))*h_d/pieces_y)+1;
            point_selet_idx=((area_idx_y-1)*w_d+area_idx_x)';%尝试改这里

            select_idx_record(k,:)=point_selet_idx;
            
            
            point_bt=point_d(:,point_selet_idx);
            point_bt_extend = [point_bt.^2;point_bt(1,:).*point_bt(2,:);point_bt;ones(1,number_sample)]; 
            point_ft=[point_bt;ones(1,number_sample)]+[flow_d(point_selet_idx,:)';zeros(1,number_sample)];
            H = point_ft*point_bt_extend'*(point_bt_extend*point_bt_extend')^(-1);
            
            
            flow_d_temp=H*point_d_extend-[point_d;ones(1,w_d*h_d)];
            flow_diff=sum(abs(flow_d_temp(1:2,:)'-flow_d),2)*10;
            flow_diff = threshold_RANSAC-flow_diff;
            idx=find(uint8(flow_diff));
            sum_diff =sum(flow_diff(idx));
            num_record(k,1)=size(idx,1);
            if(size(idx,1)>num_record_max)
                point_selet_idx_max=k;
                num_record_max = size(idx,1);
                flo_diff_max=sum_diff;
            elseif(size(idx,1)==num_record_max) 
                if(flo_diff_max>sum_diff)
                    point_selet_idx_max=k;
                    num_record_max = size(idx,1);
                    flo_diff_max=sum_diff;
                end
            end
                
        end
        %[~,ord]=sort(num_record,1,'ascend');
        point_selet_idx=select_idx_record(point_selet_idx_max,:);
        
        point_bt=point_d(:,point_selet_idx);
        point_bt_extend = [point_bt.^2;point_bt(1,:).*point_bt(2,:);point_bt;ones(1,number_sample)]; 
        point_ft=[point_bt;ones(1,number_sample)]+[flow_d(point_selet_idx,:)';zeros(1,number_sample)];
        H = point_ft*point_bt_extend'*(point_bt_extend*point_bt_extend')^(-1);
            
            
           
            
%% 动态更新阈值
        %threshold_RANSAC=sqrt(H(1,3)^2+H(2,3)^2)*para_b+para_a;
        scene_speed =mean(sqrt(sum(flow_d(point_selet_idx,:).^2,2)));
        threshold_RANSAC=scene_speed*para_b+para_a;
        ival = floor(speed/scene_speed*ival)+1;
        if ival>5
            ival =5;
        end
%% 利用H和光流图求前景mask_ori
        
        flow_temp=H*point_extend-[point;ones(1,w*h)];
        flow_temp=reshape(flow_temp(1:2,:)',w,h,2);
        
%         imgcl=flowToColor(flow_temp);
%         figure(3)
%         imshow(imgcl);
        
        flow_diff=sum(abs(flow_temp-flow),3)*10;
        mask_ori=uint8(flow_diff-threshold_RANSAC*10)*255;
%         figure(1)
%         imshow(mask_ori);
        wpath = fullfile(result_dir, seq_id, [frame_ids{f_id} '.png']);
        imwrite(mask_ori,wpath);
    end
end
[eval, raw_eval] = eval_result(result_id,{'J'},'all');
result = mean(eval.J.mean);
end
end
end
end
end
end
%end


