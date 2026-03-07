function Main_to_Call_GUI(supertitle, input_vec, analyte_name_vec, sensor_name_vec, DAT, n_measure)
Data_type_name = supertitle;
Data_type      = 1;

close all
tic

cdcd   = cd;
cdcdcd = horzcat(cdcd,'\utilities');
addpath(cdcdcd);

% =========================
% NEW (minimal): QDA regularization knobs
% =========================
DoQDA_Reg     = true;    % only applied when Decision_boundaries==0
QDA_ridgeTiny = 1e-12;   % numeric safety
QDA_cFloor    = 1.0;     % eigen-floor strength

Smart_Bin = 1; % Smart chi-squared binning, based on labels

Inflate               = input_vec(1);
STD_buff              = 0; % legacy, remains fixed
STD_buff_new          = input_vec(2);
Decision_boundaries   = input_vec(3); % 0-QDA, 1-Voronoi
chi_squared_groups    = input_vec(4); % 0-weighted/new, 1-previous/standard
lambda_1              = input_vec(5)/100;
lambda_2              = input_vec(6)/100;

% NEW: optional threshold for flagging effective 1D dominance
if numel(input_vec) >= 7
    PC1VarThresh = input_vec(7)/100;
else
    PC1VarThresh = 0.95;
end

if strcmp(supertitle,'Default')
    lablesact = 1;
elseif strcmp(supertitle,'Two-Sigma Inflation')
    lablesact = 1;
elseif strcmp(supertitle,'Dataset 1')
    lablesact = 1;
elseif strcmp(supertitle,'Two-Sigma Inflation Sweat Metabolomics')
    lablesact = 1;
elseif strcmp(supertitle,'Five-Sigma Sweat Data uFS')
    lablesact = 1;
else
    lablesact = 0;
end

mega_cda = horzcat(cd, '\results\');
kk = 1;
mm = 1;
rand_flag = 0;
faktor = 1;

if Inflate==0
    Inflate_name='_NI';
else
    Inflate_name='_I';
end

if STD_buff_new==0
    Buff_name='_Nbuff';
else
    s = strrep(num2str(STD_buff_new,'%.15g'), '.', '_');
    Buff_name=horzcat('_',s,'_buff');
end

Stop_name='';

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

if DoQDA_Reg && Decision_boundaries==0
    Reg_name = '_RegQDA';
else
    Reg_name = '';
end

config=horzcat('Config ',Data_type_name,Inflate_name,Buff_name,Stop_name,chi_name,chi_name2,Bond_name,Reg_name);

all_sensors_length = length(sensor_name_vec);
Vektor_ARI         = zeros(all_sensors_length,1);
latent_vec1        = zeros(all_sensors_length,1);
latent_vec2        = zeros(all_sensors_length,1);
avg_dist_vec       = zeros(all_sensors_length,1);
pc1_is_1d_flag     = false(all_sensors_length,1); % NEW
all_analyte_length = length(analyte_name_vec);
Gaussian_error_val_mat = zeros(all_analyte_length,all_sensors_length);
classStatss = cell(1,all_sensors_length);

mega_cd=horzcat(mega_cda,config);
if ~exist(mega_cd,'dir'), mkdir(mega_cd); end
cd(mega_cd);

DAT=DAT-mean(DAT,1); % take off all the mean
if STD_buff_new~=0
    [muMat, varMat, new_DAT] = analyte_stats_resample(DAT, n_measure, 20, 42, STD_buff_new); %#ok<ASGLU>
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

    counta=0; %#ok<NASGU>
    Marquis=['^','d','<','v','>','s'];
    colorit=[217 83 25;237 177 32;126 47 142;119 172 48;77 190 238; 100 100 3]/255;
    scorpion=zeros(2,1,size(DAT,1)); %#ok<NASGU>

    % =========================
    % SMART BINNING UPDATE:
    % regularized smart binning only in QDA mode
    % =========================
    if Smart_Bin==1
        if (Decision_boundaries==0) && DoQDA_Reg
            [Binned, zoneInfo, mu_sorted, sigma_sorted] = smart_binning_1D_reg(DAT, n_measure, all_analyte_length, Inflate); %#ok<ASGLU>
        else
            [Binned, zoneInfo, mu_sorted, sigma_sorted] = smart_binning_1D(DAT, n_measure, all_analyte_length, Decision_boundaries); %#ok<ASGLU>
        end
        BinnedDat=Binned;
    end

    [coeff,score,latent] = pca(DAT);

    prev_scoreX_signb4=sign(score(:,1));
    prev_scoreY_signb4=sign(score(:,2));
    if ll~=1
        if sum(0.5*abs(prev_scoreX_sign-sign(score(:,1))))>0.5*length(score(:,1))
            score(:,1)=-score(:,1);
            prev_scoreX_signb4=sign(score(:,1));
        end
        if sum(0.5*(abs(prev_scoreY_sign-sign(score(:,2)))))>0.5*length(score(:,2))
            score(:,2)=-score(:,2);
            prev_scoreY_signb4=sign(score(:,2));
        end
    end
    prev_scoreX_sign=prev_scoreX_signb4;
    prev_scoreY_sign=prev_scoreY_signb4;

    latent_vec1(ll)=latent(1)./sum(latent);
    latent_vec2(ll)=latent(2)./sum(latent);

    % NEW: flag 1D-dominant subsets
    pc1_is_1d_flag(ll) = latent_vec1(ll) >= PC1VarThresh;

    [Mi,Ii]=sort(abs(coeff(:,1)));
    [Mi2,Ii2]=sort(abs(coeff(:,2)));
    rank_me=zeros(all_sensors_length-ll+1,1);
    for pol=1:1:(all_sensors_length-ll+1)
        mop1=find(Ii==pol);
        mop2=find(Ii2==pol);
        inner_ranking=Mi(mop1)*latent(1)+Mi2(mop2)*latent(2);
        rank_me(pol)=inner_ranking;
    end
    [vvv,uuu]=sort(rank_me); %#ok<ASGLU>

    for lp=1:1:size(DAT,1)
        [iii,mgmm]=max(abs(coeff(:,1)'.*DAT(lp,:))); %#ok<ASGLU>
        scorpion(kk,mm,lp)=mgmm;
    end

    figgg=figure;
    bbb=get(figgg,'Position'); %#ok<NASGU>
    new_width=5.9;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.75*new_width]);
    mina=0;
    maxa=0;

    mmm2=size(DAT,1)/n_measure;
    mmm=n_measure;
    for j=1:1:mmm2
        index=(mmm*(j-1)+1):(mmm*(j-1)+mmm);
        scatter(-score(index,1),score(index,2),30,'Marker',Marquis(j), ...
            'MarkerFaceColor',"none",'MarkerEdgeColor',colorit(j,:), ...
            'MarkerEdgeAlpha',0.5,'LineWidth',0.2);
        hold on;
        mina=min(mina,min(-score(index,1))); %#ok<NASGU>
        maxa=max(maxa,max(-score(index,1))); %#ok<NASGU>
    end

    xlabel(sprintf('PC1 (%.1f%%)',100*latent_vec1(ll)));
    ylabel(sprintf('PC2 (%.1f%%)',100*latent_vec2(ll)));
    set(gca,'FontSize',6);

    X2(:,1)=-score(:,1);
    X2(:,2)= score(:,2);

    idx=zeros(1,n_measure*all_analyte_length);
    C2=zeros(all_analyte_length,2);
    for rr=1:1:all_analyte_length
        C2(rr,:)=mean(X2((1:1:n_measure)+(rr-1)*n_measure,:));
        idx((1:1:n_measure)+(rr-1)*n_measure)=rr;
    end

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
    nClasses = all_analyte_length; %#ok<NASGU>
    classStats = repmat(struct('meanXY', [], 'covXY', []), all_analyte_length, 1);

    % =========================
    % NEW: pooled target for regularized QDA
    % =========================
    if (Decision_boundaries==0) && DoQDA_Reg
        Sw = zeros(2,2);
        nTot = 0;
        for kk2 = 1:all_analyte_length
            A_pool = X2(find(idx==kk2),:); %#ok<FNDSB>
            nk2 = size(A_pool,1);
            if nk2 < 2, continue; end
            mu_pool = mean(A_pool,1);
            Xc_pool = A_pool - mu_pool;
            Sw = Sw + (Xc_pool' * Xc_pool);
            nTot = nTot + nk2;
        end
        if nTot <= 0
            S_pool = eye(2);
        else
            S_pool = Sw / nTot; % MLE pooled target
        end
        S_pool = (S_pool + S_pool')/2;
        avgVar = trace(S_pool)/2;
        avgVar = max(avgVar, QDA_ridgeTiny);
    else
        S_pool = [];
        avgVar = [];
    end

    for tt=1:1:all_analyte_length
        A=X2(find(idx==tt),:); %#ok<FNDSB>

        if Inflate==0 && STD_buff==0
            [handle_ellipse,meanXY1,covXY1] = plot_ellipse(A(:,1),A(:,2)); %#ok<ASGLU>
        elseif Inflate==1 && STD_buff==0
            nn=size(A,1);
            [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,0,Inflate); %#ok<ASGLU>
        elseif Inflate==0 && STD_buff>0
            nn=size(A,1);
            [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,STD_buff,Inflate); %#ok<ASGLU>
        else
            nn=size(A,1);
            [handle_ellipse,meanXY1,covXY1] = plot_ellipse_inflate(A(:,1),A(:,2),nn,STD_buff,Inflate); %#ok<ASGLU>
        end

        A2=A;
        X0=mean(A2(:,1));
        Y0=mean(A2(:,2));
        meanXY=[X0 Y0];

        % keep original inflation scalar logic
        factor_total = 1;
        if Inflate==0 && STD_buff==0
            factor_total = 1;
        elseif Inflate==1 && STD_buff==0
            n=nn;
            if n <= 1
                error('Need at least 2 samples for a variance estimate.');
            end
            nu = n;
            if n > 2
                faktor = ((nu-1) / (nu - 2)) * ((n+1)/n);
            else
                q = chi2inv(0.5, nu);
                faktor = (nu *((n+1)/n)) / q;
            end
            factor_total = faktor;
        elseif Inflate==0 && STD_buff>0
            faktor = STD_buff.^2;
            factor_total = faktor;
        elseif Inflate==1 && STD_buff>0
            n=nn;
            if n <= 1
                error('Need at least 2 samples for a variance estimate.');
            end
            nu = n;
            if n > 2
                faktor = ((nu-1) / (nu - 2)) * ((n+1)/n);
            else
                q = chi2inv(0.5, nu);
                faktor = (nu *((n+1)/n)) / q;
            end
            faktor = faktor.*STD_buff.^2;
            factor_total = faktor;
        end

        % =========================
        % NEW: regularized covariance in QDA mode
        % =========================
        if (Decision_boundaries==0) && DoQDA_Reg
            nk = size(A2,1);
            if nk < 2
                error('Class %d has <2 samples; cannot estimate covariance.', tt);
            end

            Xc = A2 - mean(A2,1);
            S  = (Xc' * Xc) / nk;  % MLE covariance
            S  = (S + S')/2;

            % 1) inflate
            S = S .* factor_total;

            % 2) LW shrinkage to pooled target
            [Ssh, ~] = ledoit_wolf_target_2d(Xc, S, S_pool);

            % 3) eigen-floor ridge
            eigMin = min(eig((Ssh + Ssh')/2));
            tau   = QDA_cFloor * avgVar * (1/max(nk-1,1));
            gamma = max(0, tau - eigMin);

            covXY = (Ssh + Ssh')/2 + gamma*eye(2) + QDA_ridgeTiny*eye(2);
            covXY = (covXY + covXY')/2;
        else
            covXY_raw = cov(A2(:,1),A2(:,2));
            covXY = covXY_raw .* factor_total;
        end

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

    classStats_reg = classStats;

    if ll==1
        [x1min, x2min, x1max, x2max] = boundsFromClassStats(classStats);
        classStats_first=classStats;
    end
    classStatss{:,ll}=classStats_reg;

    if Decision_boundaries==0
        Gaussian_error_vec = Gaussian_Witch_Practice_QDA(classStats_reg,'Res', 300, 'Trunc', 5);
        mister = QDA_assign_points(X2, classStats_reg, 'Res',300, 'Trunc',5);
        mister = mister';
    end

    Gaussian_error_val_mat(:,ll)=Gaussian_error_vec;
    K = numel(analyte_name_vec);
    qw = cell(1,K);

    if isstring(Marquis) || ischar(Marquis)
        M = char(Marquis);
    else
        error('Marquis must be a char vector or string of marker symbols.');
    end

    fallback = 'os^v>d<+xph*';
    for k = 1:K
        if k <= numel(M)
            mk = M(k);
        else
            mk = fallback(min(k, numel(fallback)));
        end
        cidx = min(k, size(colorit,1));
        qw{k} = scatter(nan, nan, 'Marker', mk, ...
            'MarkerEdgeColor', colorit(cidx,:),'LineWidth',0.2);
    end

    h = legend([qw{:}], analyte_name_vec, ...
        'Location','north', 'Orientation','horizontal', ...
        'NumColumns', numel(qw));

    ax = gca;
    ax.Position = ax.Position - [0 0 0 0.1];
    basePos = [0.2067 0.8655 0.6233 0.1261];
    baseK   = 6;
    K       = numel(qw);

    w = basePos(3) * (K / baseK);
    w = min(max(w, 0.25), 0.95);
    cx = basePos(1) + basePos(3)/2;
    x  = cx - w/2;

    h.Units    = 'normalized';
    h.Position = [x, basePos(2), w, basePos(4)];
    h.ItemTokenSize(1)=3;

    ylim([x2min x2max]);
    xlim([x1min x1max]);
    box on;
    print (horzcat(name_title,''),'-dpng','-r600');
    savefig(horzcat(name_title,''));

    % =========================
    % FEATURE REDUCTION UPDATE:
    % use regularized classStats in QDA-UP branches
    % =========================
    if Smart_Bin==0 && chi_squared_groups==0
        [idx2,scores] = fscchi2_norm_UP(DAT,idx,classStats_reg,'NumBins',all_analyte_length);
    elseif Smart_Bin==0 && chi_squared_groups==1
        [idx2,scores] = fscchi2_norm(DAT,idx,'NumBins',all_analyte_length);
    elseif Smart_Bin==1 && chi_squared_groups==0
        [idx2,scores] = fscchi2_norm_UP_smartBin(BinnedDat,idx,classStats_reg,'NumBins',all_analyte_length);
    else
        [idx2,scores] = fscchi2_norm_smartBin(BinnedDat,idx,'NumBins',all_analyte_length);
    end

    original_idx=idx2;
    original_scores=scores; %#ok<NASGU>
    Active_indices=original_idx;

    if Decision_boundaries==1
        name_title_Voronoi=horzcat('Surviving_Sensors_',num2str(all_sensors_length-ll+1),'_Vor_Classifier');
        plotVoronoi_Tessellation_Marquis_single(analyte_name_vec, Marquis, classStats, name_title_Voronoi, colorit, x1min, x2min, x1max, x2max);
    else
        name_title_QDA=horzcat('Surviving_Sensors_',num2str(all_sensors_length-ll+1),'_QDA_Classifier');
        plotQDA_Tessellation_Marquis_single(analyte_name_vec, Marquis, classStats_reg, name_title_QDA, colorit, x1min, x2min, x1max, x2max);
    end

    DAT_ALL2=DAT;
    DAT_ALL22=DAT_ALL2; %#ok<NASGU>

    idx2_end=idx2(end);
    DAT_ALL2(:,idx2_end) = [];
    DAT=DAT_ALL2-mean(DAT_ALL2,1);

    [RI,ARI,Dice,JD]=randindex(idx, mister); %#ok<ASGLU>
    Vektor_ARI(ll)=ARI;

    namexxxx=horzcat('Surviving_Sensors_',config,'_',num2str(all_sensors_length-ll),'.mat');
    updated_num_sensors=all_sensors_length-ll;
    char_array=char(original_names_vec);
    num_matches = 0; %#ok<NASGU>
    sensor_original_number_removed=0;
    for i = 1:size(char_array,1)
        current_matches = strfind(char_array(i,:), char(names_vec(idx2(end))));
        if ~isempty(current_matches)
            sensor_original_number_removed=i;
        end
    end
    eliminated_one(ll)=names_vec(idx2(end));
    names_vec(idx2(end)) = [];

    if rand_flag==0
        save(namexxxx,'names_vec','sensor_original_number_removed','ARI','updated_num_sensors');
    end
    Active_indices(find(original_idx==sensor_original_number_removed))=nan; %#ok<FNDSB>
    sensor_original_number_removed_up_to_now=[]; %#ok<NASGU>
end

% =========================
% Last iteration, 1 sensor
% =========================
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
    A=X2(find(idxx==tt),:); %#ok<FNDSB>
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
    mean_vec_diff=[-20 (mean_vec2(2:end)+mean_vec2(1:end-1))/2 20];
    for pop=(1:1:length(mean_vec2)-1)
        x_star = qda_boundary_between_means(mean_vec2(pop), mean_vec2(pop+1), sqrt(var_vec(pop)), sqrt(var_vec(pop+1)));
        mean_vec_diff(pop+1)=x_star;
    end
end

Gaussian_error_vec=zeros(1,all_analyte_length);
miater=idx; %#ok<NASGU>
mister_UP=zeros(1,length(mister));
for tt=1:1:all_analyte_length
    A2=X2(find(mister==tt),:); %#ok<FNDSB>
    covXY=var(A2).*faktor;
    MEAN= mean(A2);
    VAR= sqrt(covXY);
    bottom_lim=mean_vec_diff(find(trt==(tt)));
    upper_lim=mean_vec_diff(find(trt==(tt))+1);
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
    [avgErr, classErr, cuts, zones] = qda1d_gaussians(mean_vec2, var_vec); %#ok<ASGLU>
    [mister_UP, ~] = qda_predict_1d_zones(X2, cuts, zones);
    idxx2=idx_new;
    Gaussian_error_vec=classErr;
end

[RI,ARI,Dice,JD]=randindex(idxx2, mister_UP); %#ok<ASGLU>
Vektor_ARI(all_sensors_length)=ARI;
pc1_is_1d_flag(all_sensors_length)=true; % one remaining sensor is inherently 1D

eliminated_one(all_sensors_length)=names_vec(1);
save('eliminated_one','eliminated_one');
Gaussian_error_val_mat(:,all_sensors_length)=Gaussian_error_vec;
avg_dist_vec(all_sensors_length)=average_distance;

[ax11,ax22,money]=create3DAccuracyFig(config,all_sensors_length,1:1:all_sensors_length,avg_dist_vec,Vektor_ARI,Gaussian_error_val_mat,1,lambda_1,lambda_2); %#ok<ASGLU>
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

% figOut = stackPCfigs_4x4_fromFIG(h1, h2, 'First_and_Last_Iterations_PCA');
figOut = stackPCfigs_4x4_fromFIG_fusePCA(h1, h2, 'First_and_Last_Iterations_PCA');

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
     'avg_dist_vec', 'Vektor_ARI', 'Gaussian_error_val_mat', 'pc1_is_1d_flag', 'PC1VarThresh', 'supertitle', '-v7.3');
print (horzcat(config,'_Summary_Fig'),'-dpng','-r600');
savefig(horzcat(config,'_Summary_Fig'));
toc

% -------------------------------------------------------------------------
function [x1min,x2min,x1max,x2max] = boundsFromClassStats(classStats, varargin)
    p = inputParser;
    p.addParameter('NSigma', 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter('SnapToInteger', true, @(x)islogical(x) && isscalar(x));
    p.parse(varargin{:});
    nSigma = p.Results.NSigma;
    snap   = p.Results.SnapToInteger;

    K = numel(classStats);
    assert(K>0, 'classStats is empty.');

    mu    = zeros(K,2);
    evMax1= zeros(K,1);
    evMax2= zeros(K,1);
    for k = 1:K
        mu(k,:) = classStats(k).meanXY(:).';
        evMax1(k)= sqrt(max(classStats(k).covXY(1,1)));
        evMax2(k)= sqrt(max(classStats(k).covXY(2,2)));
    end
    evMax11=max(evMax1);
    evMax22=max(evMax2);
    pad1 = nSigma * evMax11;
    pad2 = nSigma * evMax22;
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
    map = trt(:);
    K   = numel(map);

    if any(idx(:) < 1 | idx(:) > K | idx(:) ~= round(idx(:)))
        error('idx must contain integers in 1..%d.', K);
    end

    sz = size(idx);
    idx_trt = reshape(map(idx(:)), sz);
end

function [mister_UP, bin] = qda_predict_1d_zones(X2, cuts, zones)
    if numel(zones) ~= numel(cuts)-1
        error('zones must have length numel(cuts)-1.');
    end
    if ~isinf(cuts(1)) || ~isinf(cuts(end))
        error('cuts must include -Inf as first and +Inf as last boundary.');
    end

    shp = size(X2);
    x = X2(:);
    edges = cuts(:).';
    [~,~,bin] = histcounts(x, edges);

    if any(isnan(bin))
        error('Some samples fell outside the provided cuts (should not happen with ±Inf).');
    end

    pred = zones(bin);
    mister_UP = reshape(pred, shp);
end

% =========================
% NEW helper: LW shrinkage to pooled target (2D)
% =========================
function [Ssh, lambda] = ledoit_wolf_target_2d(Xc, S, T)
    n = size(Xc,1);

    S = (S+S')/2;
    T = (T+T')/2;

    D = S - T;
    delta = sum(D(:).^2);

    if delta <= eps
        lambda = 1;
        Ssh = T;
        return;
    end

    beta = 0;
    for i = 1:n
        xi = Xc(i,:).';
        Pi = xi*xi.';
        E  = Pi - S;
        beta = beta + sum(E(:).^2);
    end
    beta = beta / (n^2);

    lambda = min(1, beta / (delta + eps));
    Ssh = (1-lambda)*S + lambda*T;
    Ssh = (Ssh + Ssh')/2;
end

% =========================
% NEW helper: QDA-regularized smart binning
% =========================
function [Binned, zoneInfo, mu_sorted_all, sigma_sorted_all] = smart_binning_1D_reg(DAT, n_measure, K, Inflate)
    [N,P] = size(DAT);
    Binned = zeros(N,P,'int32');

    zoneInfo = [];
    mu_sorted_all = [];
    sigma_sorted_all = [];

    epsVar = 1e-12;
    nuPrior = 1;

    for j = 1:P
        mu = zeros(K,1);
        v  = zeros(K,1);

        for k = 1:K
            ix = (1:n_measure) + (k-1)*n_measure;
            xk = DAT(ix, j);
            xk = xk(~isnan(xk));
            nk = numel(xk);

            mu(k) = mean(xk,'omitnan');

            if nk >= 2
                xc = xk - mu(k);
                v(k) = (xc'*xc) / nk; % MLE variance
            else
                v(k) = epsVar;
            end

            if Inflate==1 && nk>=2
                v(k) = v(k) * inflate_factor_var(nk);
            end

            v(k) = max(v(k), epsVar);
        end

        v_pool = mean(v);
        v_pool = max(v_pool, epsVar);

        for k = 1:K
            nk = n_measure;
            lamk = nuPrior / (nuPrior + max(nk-1,1));
            v(k) = (1-lamk)*v(k) + lamk*v_pool;

            tauk = v_pool / max(nk-1,1);
            v(k) = max(v(k), tauk);
        end

        sg = sqrt(v);

        [mu_sorted, ord] = sort(mu,'ascend');
        sg_sorted = sg(ord);

        b = zeros(K-1,1);
        for t = 1:(K-1)
            b(t) = qda_boundary_between_means(mu_sorted(t), mu_sorted(t+1), sg_sorted(t), sg_sorted(t+1));
        end

        x = DAT(:,j);
        z = int32(sum(x >= b.', 2) + 1);
        z = max(1, min(K, z));
        Binned(:,j) = int32(ord(z));

        mu_sorted_all(:,j) = mu_sorted;
        sigma_sorted_all(:,j) = sg_sorted;
    end
end

function f = inflate_factor_var(n)
    if n <= 1
        f = 1;
        return;
    end
    nu = n;
    if n > 2
        f = ((nu-1) / (nu - 2)) * ((n+1)/n);
    else
        q = chi2inv(0.5, nu);
        f = (nu * ((n+1)/n)) / q;
    end
end

end