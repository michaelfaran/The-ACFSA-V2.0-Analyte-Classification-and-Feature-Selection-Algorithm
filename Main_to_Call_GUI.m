function Main_to_Call_GUI(supertitle, input_vec,analyte_name_vec, sensor_name_vec, DAT, n_measure)
Data_type_name=supertitle;
Data_type=1;
%clear
close all
tic
%addpath 'C:\Users\admin\Documents\Run_Matlab_Fast_Folder\Gabriel\mat_files_all_data';
%addpath 'C:\Users\admin\Documents\Run_Matlab_Fast_Folder\Gabriel';
%addpath 'C:\Users\admin\Documents\Run_Matlab_Fast_Folder\Gabriel\elipse';
cdcd=cd;
cdcdcd=horzcat(cdcd,'\utilities');
addpath(cdcdcd);
Smart_Bin=1;%Smart chi-sqaured binning, based on labels
Inflate=input_vec(1); %Gaussian STD inflation for estimation error, 1 to inflate
STD_buff=0; %do not change this anymore
%STD_buff_new=floor(100*(10))./100;
STD_buff_new=input_vec(2);
Bogus=0; %This is just if we do not to rewrite the folder of recreating the same run conditions, but to do some figs etc, make it 1 if you want it
Stopping_condition=0; %0-maximum error per class, 1- weighted function
Decision_boundaries=input_vec(3);% 0- RQDA, 1- Voronoi
chi_squared_groups=input_vec(4);% 0- The default per two groups new chi sqaured metric, or 1- the previous one 
lambda_1=input_vec(5)/100;% Working point of error per sensor 1, maximum classifcation error per sensor
lambda_2=input_vec(6)/100;% Working point of error per sensor 2
if strcmp(supertitle,'Default')
lablesact=1;
elseif strcmp(supertitle,'Two-Sigma Inflation')
lablesact=1;
elseif strcmp(supertitle,'Dataset 1')
lablesact=1;
elseif strcmp(supertitle,'Two-Sigma Inflation Sweat Metabolomics')
lablesact=1;
elseif strcmp(supertitle,'Five-Sigma Sweat Data uFS')
lablesact=1;
else
lablesact=0;
end
mega_cda=cd;
mega_cda=horzcat(mega_cda,'\results\');
kk=1;
mm=1;
rand_flag=0;
faktor=1;

if Inflate==0
Inflate_name='_NI';
else
Inflate_name='_I';
end

if STD_buff_new==0
Buff_name='_Nbuff';
else
s = strrep(num2str(STD_buff_new,'%.15g'), '.', '_');   % 1.4 -> '1_4'
Buff_name=horzcat('_',s,'_buff');
end

Stop_name='';

% if Stopping_condition==0
% Stop_name=' ';
% else
% Stop_name=' ';
% end

if Decision_boundaries==0
Bond_name='_QDA';
else
Bond_name='_Vor';
end

if chi_squared_groups==0
chi_name='_WFS';
else
chi_name='_SFS';

end

chi_name2='';
% if Smart_Bin==0
% chi_name2='';
% else
% chi_name2='_SB';
% 
% end

if Bogus==0
config=horzcat('Config ',Data_type_name,Inflate_name,Buff_name,Stop_name,chi_name,chi_name2,Bond_name);
else
config='Bogus2';
end

all_sensors_length=length(sensor_name_vec);
Vektor_ARI=zeros(all_sensors_length,1);
latent_vec1=zeros(all_sensors_length,1);
latent_vec2=zeros(all_sensors_length,1);
avg_dist_vec=zeros(all_sensors_length,1);
all_analyte_length=length(analyte_name_vec);
Gaussian_error_val_mat=zeros(all_analyte_length,all_sensors_length);
classStatss=cell(1,all_sensors_length);
mega_cd=horzcat(mega_cda,config);
if ~exist(mega_cd,'dir'), mkdir(mega_cd); end
cd(mega_cd);
DAT=DAT-mean(DAT,1); %take off all the mean
if STD_buff_new~=0
[muMat, varMat, new_DAT] = analyte_stats_resample(DAT, n_measure, 20, 42, STD_buff_new);
DAT=new_DAT;
n_measure=20;
end

cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\Gabriel\mat_files_all_data');
dat_mean = squeeze(mean(reshape(DAT, n_measure, [], size(DAT,2)), 1, "omitnan"));
data_name=num2str(Data_type);
save(horzcat('dat_mean_Data_type_',data_name),'dat_mean','Data_type')
cd(mega_cd);

original_names_vec=sensor_name_vec;
names_vec=sensor_name_vec;
Active_indices=zeros(1,all_sensors_length);
eliminated_one= {};
for ll=1:1:(all_sensors_length-1)
   name_title=horzcat('PCA of ',config,' ',num2str(all_sensors_length-ll+1),' Sensors'); 
    if ll==1
    name_title_first=name_title;
    end
   counta=0;
   Marquis=['^','d','<','v','>','s'];
   colorit=[217 83 25;237 177 32;126 47 142;119 172 48;77 190 238; 100 100 3]/255;
   scorpion=zeros(2,1,size(DAT,1));
   if Smart_Bin==1
    [Binned, zoneInfo, mu_sorted, sigma_sorted] = smart_binning_1D(DAT, n_measure, all_analyte_length, Decision_boundaries);
    BinnedDat=Binned;
    end
   [coeff,score,latent] =pca(DAT);

prev_scoreX_signb4=sign(score(:,1));
prev_scoreY_signb4=sign(score(:,2));
if ll~=1
if sum(0.5*abs(prev_scoreX_sign-sign(score(:,1))))>0.5*length(score(:,1)) %updated if in one iterations some of the samples shift sign, and some not
score(:,1)=-score(:,1);
prev_scoreX_signb4=sign(score(:,1));
end
if sum(0.5*(abs(prev_scoreY_sign-sign(score(:,2)))))>0.5*length(score(:,2)) %updated if in one iterations some of the samples shift sign, and some not
score(:,2)=-score(:,2);
prev_scoreY_signb4=sign(score(:,2));
end
end
prev_scoreX_sign=prev_scoreX_signb4;
prev_scoreY_sign=prev_scoreY_signb4;   
latent_vec1(ll)=latent(1)./sum(latent);
latent_vec2(ll)=latent(2)./sum(latent);
[Mi,Ii]=sort(abs(coeff(:,1)));
[Mi2,Ii2]=sort(abs(coeff(:,2)));
rank_me=zeros(all_sensors_length-ll+1,1); 
for pol=1:1:(all_sensors_length-ll+1)
mop1=find(Ii==pol);
mop2=find(Ii2==pol);
inner_ranking=Mi(mop1)*latent(1)+Mi2(mop2)*latent(2);
rank_me(pol)=inner_ranking;
end
[vvv,uuu]=sort(rank_me); 

for lp=1:1:size(DAT,1)
[iii,mgmm]=max(abs(coeff(:,1)'.*DAT(lp,:)));
scorpion(kk,mm,lp)=mgmm;
end

figgg=figure; 
bbb=get(figgg,'Position');
new_width=5.9;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.75*new_width]);
mina=0;
maxa=0;
mmm2=size(DAT,1)/n_measure;
mmm=n_measure;
for j=1:1:mmm2
    index=(mmm*(j-1)+1):(mmm*(j-1)+mmm);
    kok=scatter(-score(index,1),score(index,2),30,'Marker',Marquis(j),'MarkerFaceColor',"none",'MarkerEdgeColor',colorit(j,:),'MarkerEdgeAlpha',0.5,'LineWidth',0.2); 
    hold on;
    mina=min(mina,min(-score(index,1)));
    maxa=max(maxa,max(-score(index,1)));
end
ylabel('PC2');
xlabel('PC1');
set(gca,'FontSize',6);
X2(:,1)=-score(:,1);
X2(:,2)=score(:,2);
idx=zeros(1,n_measure*all_analyte_length);
C2=zeros(all_analyte_length,2);
for rr=1:1:all_analyte_length
    C2(rr,:)=mean(X2((1:1:n_measure)+(rr-1)*n_measure,:));
    idx((1:1:n_measure)+(rr-1)*n_measure)=rr;
end
bs_ext=[-20 -20 20 20 ; 20 -20 -20 20]';%was before counter clockwise
bs_ext=[-20 -20 20 20 ; -20 20 20 -20]';
[Vv,CC,XY]=VoronoiLimit(C2(:,1),C2(:,2),'bs_ext',bs_ext,'figure','off');
oo2=zeros(size(XY,1),1);
for kdk=1:1:size(XY,1)
    I=find(XY(kdk,1)==C2(:,1));
    oo2(kdk)=I;    
end
[average_distance]= calc_avg_dist(C2);
avg_dist_vec(ll)=average_distance;
Gaussian_error_vec=zeros(1,all_analyte_length);
mister=zeros(1,size(X2,1));
nClasses = all_analyte_length;
classStats = repmat(struct('meanXY', [], 'covXY', []), nClasses, 1);

for tt=1:1:all_analyte_length
    A=X2(find(idx==tt),:);
    if Inflate==0 && STD_buff==0
        [handle_ellipse,meanXY1,covXY1] = plot_ellipse(A(:,1),A(:,2));
    elseif Inflate==1 && STD_buff==0
        nn=size(A,1);
        [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,0,Inflate);
    elseif Inflate==0 && STD_buff>0
         nn=size(A,1);
        [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,STD_buff,Inflate);
    else
         nn=size(A,1);
        [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,STD_buff,Inflate);   
end

A2=A;
X0=mean(A2(:,1));
Y0=mean(A2(:,2));
if Inflate==0 && STD_buff==0
    covXY=cov(A2(:,1),A2(:,2));
elseif Inflate==1 && STD_buff==0
    n=nn;
    if n <= 1
        error('Need at least 2 samples for a variance estimate.');
    end    
    nu = n;  % degrees of freedom
    if n > 3
        faktor = ((nu-1) / (nu - 2)) * ((n+1)/n) ; %factor for the varaince
    else
        q = chi2inv(0.5, nu);         % median of χ²_{nu}
        faktor = (nu *((n+1)/n)) / q;
    end
        covXY=cov(A2(:,1),A2(:,2)).*faktor;
elseif Inflate==0 && STD_buff>0
    faktor=STD_buff.^2;
    covXY=cov(A2(:,1),A2(:,2)).*faktor;
elseif Inflate==1 && STD_buff>0
    n=nn;
    if n <= 1
        error('Need at least 2 samples for a variance estimate.');
    end    
    nu = n;  % degrees of freedom
    if n > 2
        faktor = ((nu-1) / (nu - 2)) * ((n+1)/n) ; %factor for the varaince
    else
        q = chi2inv(0.5, nu);         % median of χ²_{nu}
        faktor = (nu *((n+1)/n)) / q;
    end
        faktor=faktor.*STD_buff.^2;
        covXY=cov(A2(:,1),A2(:,2)).*faktor;
end
meanXY=[X0 Y0];
handle_ellipse.Color=colorit(tt,:);  
handle_ellipse.LineStyle=':'; 
if Decision_boundaries==1
    [Gaussian_error] = Gaussian_Witch_Practice(meanXY,covXY,Vv(CC{(find(oo2==tt)),:},:));
    Gaussian_error_vec(tt)=Gaussian_error;
    A2_idx=Gaussian_Witch_Practice2(X2,Vv(CC{(find(oo2==tt)),:},:));
    mister(A2_idx==1)=tt;
end
classStats(tt).meanXY=meanXY;
classStats(tt).covXY=covXY;
end

if ll==1
    [x1min, x2min, x1max, x2max] = boundsFromClassStats(classStats);
    classStats_first=classStats;
end
classStatss{:,ll}=classStats;
if Decision_boundaries==0 %QDA
Gaussian_error_vec = Gaussian_Witch_Practice_QDA(classStats,'Res', 300, 'Trunc', 5);
mister = QDA_assign_points(X2, classStats, 'Res',300, 'Trunc',5);
mister =mister';
end
Gaussian_error_val_mat(:,ll)=Gaussian_error_vec;
K = numel(analyte_name_vec);
qw = cell(1,K);
    
% Normalize Marquis to a char vector
if isstring(Marquis) || ischar(Marquis)
    M = char(Marquis);
else
    error('Marquis must be a char vector or string of marker symbols.');
end

% Fallback marker sequence if Marquis is shorter than K
fallback = 'os^v>d<+xph*';

for k = 1:K
    if k <= numel(M)
        mk = M(k);
    else
        mk = fallback(min(k, numel(fallback))); % safe fallback
    end
    % Guard against short colorit
    cidx = min(k, size(colorit,1));
    qw{k} = scatter(nan, nan, 'Marker', mk, ...
        'MarkerEdgeColor', colorit(cidx,:),'LineWidth',0.2);
end

h = legend([qw{:}], analyte_name_vec, ...
    'Location','north', 'Orientation','horizontal', ...
    'NumColumns', numel(qw));   % force one-line legend
% Keep your axis tweak
ax = gca;
ax.Position = ax.Position - [0 0 0 0.1];
% Base legend box (your original numbers)
basePos = [0.2067 0.8655 0.6233 0.1261];   % [x y w h] in normalized units
baseK   = 6;                                % reference count
K       = numel(qw);
% Scale width proportional to item count, clamp to sane bounds
w = basePos(3) * (K / baseK);
w = min(max(w, 0.25), 0.95);               % clamp between 0.25 and 0.95
% Keep legend horizontally centered around original center
cx = basePos(1) + basePos(3)/2;
x  = cx - w/2;
% Apply new position (same y & height)
h.Units    = 'normalized';
h.Position = [x, basePos(2), w, basePos(4)];
h.ItemTokenSize(1)=3;
ylim([x2min x2max]);
xlim([x1min x1max]);
box on;
print (horzcat(name_title,''),'-dpng','-r600');
savefig(horzcat(name_title,''));

if Smart_Bin==0 && chi_squared_groups==0 % BinnedDat, Smart_Bin as well no here, but in the else
    [idx2,scores] = fscchi2_norm_UP(DAT,idx,classStats,'NumBins',all_analyte_length);
elseif Smart_Bin==0 && chi_squared_groups==1
    [idx2,scores] = fscchi2_norm(DAT,idx,'NumBins',all_analyte_length);
elseif Smart_Bin==1 && chi_squared_groups==0
    [idx2,scores] = fscchi2_norm_UP_smartBin(BinnedDat,idx,classStats,'NumBins',all_analyte_length);
else
    [idx2,scores] = fscchi2_norm_smartBin(BinnedDat,idx,'NumBins',all_analyte_length);
end

original_idx=idx2;
original_scores=scores;
Active_indices=original_idx;
if Decision_boundaries==1
    name_title_Voronoi=horzcat('Surviving_Sensors_',num2str(all_sensors_length-ll+1),'_Vor_Classifier'); 
    plotVoronoi_Tessellation_Marquis_single(analyte_name_vec, Marquis, classStats, name_title_Voronoi, colorit, x1min, x2min, x1max, x2max);
else
    name_title_QDA=horzcat('Surviving_Sensors_',num2str(all_sensors_length-ll+1),'_QDA_Classifier'); 
    plotQDA_Tessellation_Marquis_single(analyte_name_vec, Marquis, classStats, name_title_QDA, colorit, x1min, x2min, x1max, x2max);
end
  
DAT_ALL2=DAT;
DAT_ALL22=DAT_ALL2;

idx2_end=idx2(end);
DAT_ALL2(:,idx2_end) = [];
DAT=DAT_ALL2-mean(DAT_ALL2,1);

[RI,ARI,Dice,JD]=randindex(idx, mister);
Vektor_ARI(ll)=ARI; 
%remaining sensors
namexxxx=horzcat('Surviving_Sensors_',config,'_',num2str(all_sensors_length-ll),'.mat');
updated_num_sensors=all_sensors_length-ll;
char_array=char(original_names_vec);
num_matches = 0;
sensor_original_number_removed=0;
for i = 1:size(char_array,1)
     % Find occurrences in each cell element using strfind
     current_matches = strfind(char_array(i,:), char(names_vec(idx2(end))));
     % Update num_matches based on non-empty occurrences
     if ~isempty(current_matches)
            sensor_original_number_removed=i; 
     end
end
eliminated_one(ll)=names_vec(idx2(end));
names_vec(idx2(end)) = [];

if rand_flag==0
    save(namexxxx,'names_vec','sensor_original_number_removed','ARI','updated_num_sensors'); 
end
Active_indices(find(original_idx==sensor_original_number_removed))=nan;
sensor_original_number_removed_up_to_now=[];
end 

C2=zeros(all_analyte_length,1);
for rr=1:1:all_analyte_length
C2(rr,:)=mean(DAT_ALL2((1:1:n_measure)+(rr-1)*n_measure,:));
idx((1:1:n_measure)+(rr-1)*n_measure)=rr;
end
idxx=idx;
idxx2=idx;
[average_distance]= calc_avg_dist([C2 zeros(1,length(C2))']);
Idxx_idx=find(mister-idxx2~=0);
idxx2(Idxx_idx)=all_analyte_length+1;
Gaussian_error_vec=zeros(1,all_analyte_length);
X2=DAT_ALL2-mean(DAT_ALL2,1);
mean_vec=zeros(1,all_analyte_length);
var_vec=zeros(1,all_analyte_length);
for tt=1:1:all_analyte_length
    A=X2(find(idxx==tt),:);
    mean_vec(tt)=mean(A);
    try
        var_vec(tt)=faktor*var(A);
    catch
        var_vec(tt)=1*var(A);
    end
end
[mean_vec2,trt]=sort(mean_vec);
var_vec=var_vec(trt);

if Decision_boundaries==1
    mean_vec_diff=[-20 (mean_vec2(2:end)+mean_vec2(1:end-1))/2 20];
elseif Decision_boundaries==0
    mean_vec_diff=[-20 (mean_vec2(2:end)+mean_vec2(1:end-1))/2 20]; %just to initalize
    for pop=(1:1:length(mean_vec2)-1)
        x_star = qda_boundary_between_means(mean_vec2(pop), mean_vec2(pop+1), sqrt(var_vec(pop)), sqrt(var_vec(pop+1)));
        mean_vec_diff(pop+1)=x_star;
    end
end

Gaussian_error_vec=zeros(1,all_analyte_length);
miater=idx;
mister_UP=zeros(1,length(mister));
for tt=1:1:all_analyte_length
    A2=X2(find(mister==tt),:);
    covXY=var(A2).*faktor;
    MEAN= mean(A2);
    VAR= sqrt(covXY);
    bottom_lim=mean_vec_diff(find(trt==(tt)));
    upper_lim=mean_vec_diff( find(trt==(tt))+1);
    cp = normcdf([bottom_lim upper_lim], MEAN, VAR);
    Prob = cp(2) - cp(1);
    Gaussian_error_vec(tt)=1-Prob; 
    for momo=1:1:length(X2)
        if bottom_lim<X2(momo) && X2(momo)<upper_lim
            mister_UP(momo)=tt;
        end
    end
end
if Decision_boundaries==0
    idx_new = apply_trt_mapping(idx, trt);
    [avgErr, classErr, cuts, zones] = qda1d_gaussians(mean_vec2, var_vec);
    [mister_UP, ~] = qda_predict_1d_zones(X2, cuts, zones);
    idxx2=idx_new;
    Gaussian_error_vec=classErr;
end

[RI,ARI,Dice,JD]=randindex(idxx2, mister_UP);
Vektor_ARI(all_sensors_length)=ARI;
eliminated_one(all_sensors_length)=names_vec(1);
save('eliminated_one','eliminated_one'); 
Gaussian_error_val_mat(:,all_sensors_length)=Gaussian_error_vec;        
avg_dist_vec(all_sensors_length)=average_distance;
[ax11,ax22,money]=create3DAccuracyFig(config,all_sensors_length,1:1:all_sensors_length,avg_dist_vec,Vektor_ARI,Gaussian_error_val_mat,1,lambda_1,lambda_2);
close all;
h1=open(horzcat(name_title_first,'.fig'));
try
    name_title_UP=horzcat('PCA of ',config,' ',num2str(money),' Sensors'); 
    h2=open(horzcat(name_title_UP,'.fig'));
    classStats=classStatss{:,all_sensors_length+1-money};
catch
    name_title_UP=horzcat('PCA of ',config,' ',num2str(money)+1,' Sensors'); 
    h2=open(horzcat(name_title_UP,'.fig'));
    classStats=classStatss{:,all_sensors_length-money};
end

figOut = stackPCfigs_4x4_fromFIG(h1, h2, 'First_and_Last_Iterations_PCA');
name_title_QDA='QDA_Classifier_First_and_Last_Iteration';

if Decision_boundaries==0
    name_title_QDA='QDA_Classifier_First_and_Last_Iteration';
    plotQDA_Tessellation_Marquis_double_v2(analyte_name_vec, Marquis, classStats_first,classStats,name_title_QDA,colorit,x1min, x2min, x1max, x2max);
else
    name_title_QDA='Vor_Classifier_First_and_Last_Iteration';
    plotVoronoi_Tessellation_Marquis_double_v2(analyte_name_vec, Marquis, classStats_first, classStats, name_title_QDA, colorit, x1min, x2min, x1max, x2max);
end

figOut = collage1x3CopyPlus1Create([figure(3) figure(4)], config, all_sensors_length, 1:1:all_sensors_length,avg_dist_vec, Vektor_ARI, Gaussian_error_val_mat,supertitle,lablesact,lambda_1,lambda_2);
fname=supertitle;
sensor_idx = 1:all_sensors_length;
save(fname, 'config', 'all_sensors_length', 'sensor_idx', ...
     'avg_dist_vec', 'Vektor_ARI', 'Gaussian_error_val_mat', 'supertitle', '-v7.3');
print (horzcat(config,'_Summary_Fig'),'-dpng','-r600');
savefig(horzcat(config,'_Summary_Fig'));
toc

function [x1min,x2min,x1max,x2max] = boundsFromClassStats(classStats, varargin)
% boundsFromClassStats  Compute 2D plot bounds from Gaussian class models.
%
%   [x1min,x2min,x1max,x2max] = boundsFromClassStats(classStats)
%   assumes each classStats(k) has fields:
%       .meanXY  (1x2)   class mean
%       .covXY   (2x2)   class covariance
%
%   Name-value options:
%       'NSigma'         - padding in units of std along max-variance axis (default 3)
%       'SnapToInteger'  - snap bounds to floor/ceil integers (default true)
%
%   Example:
%       [x1min,x2min,x1max,x2max] = boundsFromClassStats(classStats,'NSigma',2.5);

    p = inputParser;
    p.addParameter('NSigma', 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter('SnapToInteger', true, @(x)islogical(x) && isscalar(x));
    p.parse(varargin{:});
    nSigma = p.Results.NSigma;
    snap   = p.Results.SnapToInteger;

    K = numel(classStats);
    assert(K>0, 'classStats is empty.');
    
    % Extract means and covariances (symmetrize cov, ensure PSD numerically)
    d = 2;
    mu    = zeros(K,d);
    evMax1= zeros(K,1); % std along principal axis (sqrt of max eigenvalue)
    evMax2= zeros(K,1); % std along principal axis (sqrt of max eigenvalue)  
    for k = 1:K
        mu(k,:) = classStats(k).meanXY(:).';
        evMax1(k)= sqrt(max(classStats(k).covXY(1,1)));  % std along the most variant direction
        evMax2(k)= sqrt(max(classStats(k).covXY(2,2)));  % std along the most variant direction  
    end
    evMax11=max(evMax1);
    evMax22=max(evMax2);
    pad1 = nSigma * evMax11; % per-class radial padding
    pad2 = nSigma * evMax22; % per-class radial padding    
    x1min = min(mu(:,1) - pad1);
    x1max = max(mu(:,1) + pad1);
    x2min = min(mu(:,2) - pad2);
    x2max = max(mu(:,2) + pad2);

    if snap
        x1min = floor(4*x1min)./4;
        x2min = floor(4*x2min)./4;
        x1max = ceil(4*x1max)./4;
        x2max = ceil(4*x2max)./4;
    end
end

function idx_trt = apply_trt_mapping(idx, trt)
% APPLY_TRT_MAPPING  Replace each label i in idx with trt(i).
%   idx_trt = apply_trt_mapping(idx, trt)
%   - idx : vector/array of integer class labels in 1..numel(trt)
%   - trt : 1xK (or Kx1) mapping, e.g. [4 2 5 1 3 6]
%
% Example:
%   idx = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6];
%   trt = [4 2 5 1 3 6];
%   idx_trt = apply_trt_mapping(idx, trt)
%   % -> [4 4 4 2 2 2 5 5 5 1 1 1 3 3 3 6 6 6]

    map = trt(:);
    K   = numel(map);

    if any(idx(:) < 1 | idx(:) > K | idx(:) ~= round(idx(:)))
        error('idx must contain integers in 1..%d.', K);
    end

    sz = size(idx);
    idx_trt = reshape(map(idx(:)), sz);
end

function [mister_UP, bin] = qda_predict_1d_zones(X2, cuts, zones)
% QDA_PREDICT_1D_ZONES  Classify 1-D samples by QDA decision intervals.
%   [mister_UP, bin] = qda_predict_1d_zones(X2, cuts, zones)
%   - X2   : vector of samples (Nx1 or 1xN)
%   - cuts : boundaries from qda1d_gaussians (column or row), MUST include -Inf and +Inf
%   - zones: interval labels for each open interval (length = numel(cuts)-1)
% Returns:
%   - mister_UP : predicted labels (same shape as X2)
%   - bin       : interval index for each sample (1..numel(cuts)-1)

    % sanity checks
    if numel(zones) ~= numel(cuts)-1
        error('zones must have length numel(cuts)-1.');
    end
    if ~isinf(cuts(1)) || ~isinf(cuts(end))
        error('cuts must include -Inf as first and +Inf as last boundary.');
    end

    shp = size(X2);
    x = X2(:);

    % histcounts assigns each x to an interval (left-closed, right-open; last bin right-closed)
    edges = cuts(:).';
    [~,~,bin] = histcounts(x, edges);

    if any(isnan(bin))
        error('Some samples fell outside the provided cuts (should not happen with ±Inf).');
    end

    pred = zones(bin);
    mister_UP = reshape(pred, shp);
end

end