function [ax1,ax2,money]=create3DAccuracyFig(config,all_sensors_length,number_of_sensors,avg_dist_vec,Vektor_ARI,Gaussian_error_val,no_man,lambda_1,lambda_2)

% Create some data to work with
x = flip(number_of_sensors,2); 
y1 = avg_dist_vec; 
y2 = Vektor_ARI(1:all_sensors_length); 
y3 = 100*mean(Gaussian_error_val); 

% Plot on the left and right y axes
if no_man==1
figgg=figure; 
bbb=get(figgg,'Position');
new_width=5.9;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.75*new_width]);
ax1 = axes;
ax1.Position=ax1.Position- [-0.05 -0.1-0.03 +0.1 0.2+0.1];
else
    ax1=no_man;
end

yyaxis left   
% see [1]


plot(flip(number_of_sensors,2),avg_dist_vec(1:all_sensors_length),'.','Marker','o','MarkerSize',3);
xlabel('Remaining Sensors');
% ylabel('Mean Clusters Distance [a.u.]');
% ylim([0 ceil(max(y1))]);
ylim([0 3]);
pause(0.1)                  % see [3]
% set the y(left) and x tick values, make them permanent 
% This is the tricky part and shoudl receive a lot of thought when 
% you adapt this to your code...
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
xlim([0 all_sensors_length]);
grid(ax1,'off')
ytick = ax1.YTick;  
yyaxis right                % see [1]
ff=plot(ax1,x,y3,'.');
 yyy=ylim;
ylim([-2 100]);

hold on;
A =scatter(nan,nan,10,'Marker','o','MarkerEdgeColor',[0 0.4470 0.7410]); hold on;
B= scatter(nan,nan,10,'Marker','diamond','MarkerEdgeColor','k');hold on;
h = legend([A B ff], {'$\langle D \rangle$', 'ARI', 'Classifier Error [\%]'}, 'Location', 'north', 'Orientation', 'horizontal', 'Interpreter', 'latex');
h.Position=h.Position+[0 0.2 0 0];
      set(gca,'FontSize',6);
% create 2nd, transparent axes
ax2 = axes('position', ax1.Position);
plot(x,y2,'.k','Marker','diamond','MarkerSize',3);
ylim([0 1]);
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
grid(ax2, 'off')
% Horizontally scale the y axis to alight the grid (again, be careful!)
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ax2.YTick = round(linspace(yl(1), yl(2), length(ytick))*100)/100;      % see [2]
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'    '}); 
set ( ax1, 'xdir', 'reverse' );
set ( ax2, 'xdir', 'reverse' );
h.ItemTokenSize(1)=6;
h.Position(3)=0.6;
h.Position= [0.2067    0.8776-0.05    0.6233    0.1261];
set(gca,'FontSize',6);
if no_man==1
Namexs=horzcat('Accuracy_vs_sensor_',config,'_3_YAxis');
print (Namexs,'-dpng','-r300');
savefig(Namexs);
else
h.Position=[0.633    0.4681    0.2    0.0461];
end

%add two different stopping conditions:
%mean error overall
hold on;
[M,I1]=min(y3./100+lambda_1*flip(number_of_sensors));
xline(length(number_of_sensors)-I1+1,'r--');
hold on;
[M,I2]=min(y3./100+lambda_2*flip(number_of_sensors));
money=length(number_of_sensors)-I2+1;
xline(length(number_of_sensors)-I2+1,'r-');
%mean error+number of sensorsclose all
% Map indices -> remaining sensor counts (abscissa)
N_wp1 = length(number_of_sensors) - I1 + 1;
N_wp2 = length(number_of_sensors) - I2 + 1;

% Mean classification errors at those points (already in %)
err_wp1 = y3(I1);
err_wp2 = y3(I2);

% ======= Write output.txt ===============================================
try
    reportFile = fullfile(pwd,'output.txt');
    fid = fopen(reportFile,'w'); 
    
    if fid==-1
        warning('Could not open output.txt for writing.');
    else
        % Header
        fprintf(fid,'ACFSA V2.0 run summary\n');
        fprintf(fid,'Config: %s\n', config);
        fprintf(fid,'\nInterpretation: Two working points are shown as vertical lines in %s_Summary_Fig.png (dashed = WP1, solid = WP2).\n\n', config);

        % WP1 (Alternative; dashed)
        fprintf(fid,'[WP1 – Alternative] lambda_1 = %.3f%%\n', 100*lambda_1);
        fprintf(fid,'  Minimal sensor set size: %d\n', N_wp1);
        fprintf(fid,'  Mean classification error at WP1: %.4f%%\n', err_wp1);
        names_wp1 = localGetNamesForN(config,N_wp1);
        fprintf(fid,'  Remaining sensors (%d): %s\n\n', numel(names_wp1), localJoinNames(names_wp1));

        % WP2 (Default; solid)
        fprintf(fid,'[WP2 – Default]     lambda_2 = %.3f%%\n', 100*lambda_2);
        fprintf(fid,'  Minimal sensor set size: %d\n', N_wp2);
        fprintf(fid,'  Mean classification error at WP2: %.4f%%\n', err_wp2);
        names_wp2 = localGetNamesForN(config,N_wp2);
        fprintf(fid,'  Remaining sensors (%d): %s\n', numel(names_wp2), localJoinNames(names_wp2));

        fclose(fid);
        fprintf('Wrote %s\n', reportFile);
    end
catch ME
    warning('Failed to write output.txt: %s', ME.message);
end

end % ======= function end =======


% ---------- helpers (nested) -------------------------------------------
function names = localGetNamesForN(cfg, N)
% Return remaining sensor names for a working point with N remaining sensors.
% If N==1, read the last entry from eliminated_one.mat (the survivor).
    if N > 1
        f = sprintf('Surviving_Sensors_%s_%d.mat', cfg, N);
        if exist(f,'file')
            S = load(f,'names_vec');
            if isfield(S,'names_vec')
                names = S.names_vec(:).';
                return
            end
        end
        names = {};
    else
        % N == 1
        f = 'eliminated_one.mat';
        if exist(f,'file')
            E = load(f,'eliminated_one');
            if isfield(E,'eliminated_one') && ~isempty(E.eliminated_one)
                survivor = E.eliminated_one{end};
                names = {survivor};
                return
            end
        end
        names = {};
    end
end

function s = localJoinNames(c)
% Join cell array of names into a readable comma-separated string.
    if isempty(c)
        s = 'N/A';
    else
        try
            s = strjoin(string(c), ', ');
        catch
            % fallback for older MATLAB
            tmp = cellfun(@char, c, 'UniformOutput', false);
            s = strjoin(tmp, ', ');
        end
    end
end