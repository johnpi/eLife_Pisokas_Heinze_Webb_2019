%%
%% Plots a minimal version of connectivity weights
%%
%% Usage example
%% 
%% The Drosophila connectivity strengths matrix
%% w_vec_drosophila = [6.07203896674184	47.9173822984772	-35.6772759590409	-19.3358579258735	99.9680554605268] % Good
%% f1 = figure();
%% plot_con_matrix_minimal(w_vec_drosophila, 'Connectivity Matrix of Drosophila', 'Input', 'Output');
%% filename = './Reports/Report-6/connectivity-matrix-simplified-drosophila';
%% title('');
%% set(gca, 'FontSize', 16);
%% set(findall(gcf, 'type', 'text'), 'FontSize', 16);
%% fig = gcf;
%% fig.PaperUnits = 'inches';
%% fig.PaperPosition = [0 0 7 6]; % Outer Image size in inches
%% print(filename, '-dpng', '-r300');        % Image resolution
%% print(filename, '-dpdf', '-r300');        % Image resolution
%% 
%% The Locust connectivity strengths matrix
%% w_vec_locust = [10.8341743712825	30.7047329203947	-42.1014161760853	-9.78202338418053	20.2648078601249] % Good gradual (selected?)
%% f2 = figure();
%% plot_con_matrix_minimal(w_vec_locust, 'Connectivity Matrix of Locust', 'Input', 'Output');
%% filename = './Reports/Report-6/connectivity-matrix-simplified-locust';
%% title('');
%% set(gca, 'FontSize', 16);
%% set(findall(gcf, 'type', 'text'), 'FontSize', 16);
%% fig = gcf;
%% fig.PaperUnits = 'inches';
%% fig.PaperPosition = [0 0 7 6]; % Outer Image size in inches
%% print(filename, '-dpng', '-r300');        % Image resolution
%% print(filename, '-dpdf', '-r300');        % Image resolution
%%
%% OR
%%
%% For saving them in the format PLOS requires:
%% f=getframe(f1);
%% imwrite(f.cdata, '/Users/john/Documents/MATLAB/Exp_Recurrent_PB_EB/Data/Plots/connectivity-matrix-simplified-drosophila.tif', 'Resolution', 300);
%% f=getframe(f2);
%% imwrite(f.cdata, '/Users/john/Documents/MATLAB/Exp_Recurrent_PB_EB/Data/Plots/connectivity-matrix-simplified-locust.tif', 'Resolution', 300);
%% 

function plot_con_matrix_minimal(w_vec, title_str, xlabel_str, ylabel_str, z_range_min, z_range_max)
    if length(w_vec) == 5
        cols_rows_vec = {'P-EN', 'P-EG', 'E-PG', 'Delta7'};
        %cols_rows_vec = ['P-EN' 'P-EG' 'E-PG' 'P_{intr}'];
        con_matrix = [
        [0           0           w_vec(5)    w_vec(3)];
        [0           0           w_vec(5)    w_vec(3)];
        [w_vec(1)    w_vec(1)    0           0       ];
        [0           0           w_vec(2)    w_vec(4)];
        ];
    end
    con_matrix
    ticks = [1:size(con_matrix, 1)];
    %[x, y] = meshgrid([1:length(con_matrix)], [1:length(con_matrix)]);
    %figure();
    % Plot the connectivity weight matrix
    %p = pcolor(x, y, con_matrix);
    %pcolor(con_matrix);
    imagesc(con_matrix);
    xticks(ticks)
    xticklabels(cols_rows_vec)
    yticks(ticks)
    yticklabels(cols_rows_vec)
    grid off
    %view(0, 90); % azimuth, elevation
    colorbar;     % show scale colorbar
    %caxis([-20 20]);
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title(title_str);
    % Set the range of values (colors)
    if ~exist('z_range_min', 'var')
        z_range_min = -100;
    end
    if ~exist('z_range_max', 'var')
        z_range_max = 100;
    end
    caxis([z_range_min z_range_max]);
    
    %% Set the zero value to correspond to black color in the colormap
    %cmap = colormap; % Get the colormap values array
    %cmap_portion = z_range_max / (z_range_max - z_range_min); % in [0, 1]
    %indx = round(cmap_portion * size(cmap, 1)); % Index into the colormap corresponding to data value 0
    %cmap(indx, :) = [0.5 0.5 0.5];
    %cmap(indx-1, :) = [0.5 0.5 0.5];
    %cmap(indx+1, :) = [0.5 0.5 0.5];
    %colormap(cmap); % set the new array as the colormap
    
    % Blue-Black-Red scale
     m=zeros(101,3);
     b_1to50=linspace(1,0.2,50)';
     r_51to100 = linspace(0.2,1,50)';
     m(1:50,3) = b_1to50;
     m(52:101,1) = r_51to100;
     colormap(m)

    % Blue-White-Red scale
     m=ones(101,3);
     b_1to50=linspace(0.2,1,50)';
     r_51to100 = linspace(1,0.2,50)';
     m(1:50,1) = b_1to50;
     m(1:50,2) = b_1to50;
     m(52:101,2) = r_51to100;
     m(52:101,3) = r_51to100;
     colormap(m)
     
     m = brewermap(241, 'RdBu');
     m(100:119, :) = [];
     m(102:121, :) = [];
     m = flipud(m);
     colormap(m);

end
