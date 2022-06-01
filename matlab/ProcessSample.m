
function ProcessSample(image_paths, dest_dir, config)
    % disp(append("working on: ", images_dir));

    config = struct(config)

    mkdir(dest_dir);
    cd(dest_dir);

    % read in the images
    images = LoadImages(image_paths);

    % config.maxSpeed = config.scale * 40;


    % configure the configuration
    if strcmp(config.bac, 'first')
        config.bac = images(1);
    end

    if strcmp(config.threshold, 'auto')
        config.threshold = CalculateThreshold(images);
    end


    % run the front tracking routine and write results
    fprintf("Starting to track fronts (threshold=%d)....\n", config.threshold);
    fronts_gdf = fullfile(dest_dir, 'fronts.gdf');
    try
        fronts = TheFrontTrackerArrays(images, ...
                                 config.threshold, ...
                                 config.step, ...
                                 config.minarea, ...
                                 '', ...
                                 config.span, ...
                                 config.ROI, ...
                                 config.framerange, ...
                                 config.scale, ...
                                 config.frrate, ...
                                 config.invert, ...
                                 fronts_gdf, ...
                                 config.noisy, ...
                                 config.thick, ...
                                 config.bac, ...
                                 config.maxSpeed, ...
                                 config.maxPx);
        save('fronts.mat', 'fronts');
        index_gdf(fronts_gdf);
        % ConvertImagesToVideo(fullfile(dest_dir, 'frontsmovie'), fullfile(dest_dir, 'frontsmovie.avi'));
        disp("done.");

        %run the visualization routine and write results
        try
            [wx,wy,vx,vy,x,y,t,fi,th] = VF(fronts_gdf);
            saveas(gcf, 'vf', 'png');
            saveas(gcf, 'vf', 'pdf');
            saveas(gcf, 'vf', 'fig');
            writetable(table(t,fi,x,y,vx,vy,wx,wy,th), 'vf.csv');
        catch ERR
            % sometimes chokes on plotting, run without plot so we can get
            % the tabluar data out!
            [wx,wy,vx,vy,x,y,t,fi,th] = VF(fronts_gdf, [-inf inf], [0 0 inf inf], 0);
            writetable(table(t,fi,x,y,vx,vy,wx,wy,th), 'vf.csv');
            rethrow(ERR);
        end
    catch ERR
        disp("Error running front tracking:");
        disp(ERR);
        rethrow(ERR);
    end
end

% function config = CreateDefaultTrackingConfig()
%     config = containers.Map();
%     config('threshold') = 50; % int or 'auto'
%     config('step') = 1;
%     config('minarea') = 25;
%     config('span') = 50;
%     config('ROI') = []; % [x y width height], in pixels
%     config('framerange') = [-inf inf];
%     config('frrate') = 1/60;
%     config('invert') = 0;
%     config('noisy') = 0;
%     config('thick') = 1;
%     config('maxSpeed') = NaN;
%     config('maxPx') = 40;
%     config('bac') = 0; % number, 'first' (first image in series), or another image array
    
%     config('use_im_reg') = 0; % Use image registration? 1=Yes; 0=No
% end


function threshold = CalculateThreshold(images)
    [counts, binLocations] = imhist(images);
    %[knee, kneeI] = knee_pt(binLocations, counts);
    kneeI = triangle_th(counts, length(binLocations));
    threshold = binLocations(kneeI*length(binLocations));
end


function ConvertImagesToVideo(imdir, dest)
    images = dir(fullfile(imdir, '*.png'));
    outputVideo = VideoWriter(dest, 'Motion JPEG AVI');
    outputVideo.FrameRate = 2;
    outputVideo.Quality = 100;
    open(outputVideo);
    for ii = 1:length(images)
        img = imread(fullfile(images(ii).folder, images(ii).name));
        writeVideo(outputVideo, img);
    end
    close(outputVideo)
end


function images = LoadImages(image_paths)
    nfiles = length(image_paths);    % Number of files found
    %images = uint8.empty(img_height, img_width, 0);
    for ii=1:nfiles
        currentimage = imread(string(image_paths(ii)));
        if size(currentimage, 3) == 3
            currentimage = rgb2gray(currentimage);
        end
        images(:,:,ii) = currentimage;
    end
end
