%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for analyzing Morphometrics-processed pombe images.
% The general outline for identifying cells and tracking phylogeny is to
% answer some basic questions:
%
% 1. Is the identified 'cell', i.e. object with a CellID, an actual cell?
%    This is based on some basic histograms to apply limits, including
%    cell width, pole-pole length, and size.
% 2. Are there neighbor cells, and are they related to the given cell? 
%    This requires some cuts based on center of mass, parallel pole-pole 
%    vectors.
% 3. If there are 'related' neighbor cells are they (i) the same cell
%    growing or (ii) the cell dividing? If there are no neighbor cells, then
%    either phylogeny is lost in the next frame, or the one before.
%
% The outputs will be placed in a phylogenetic tree to answer questions
% about cell shape and growth rates. Another script will analyze and make
% plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Known constraints
% - Based on a good segmentation of identified 8300 cell IDs, the mean
%   eccentricity (pole-pole/width) is 2.3638 +/- 1.1388. The mean area
%   is 1912.2 +/- 448.7585. Eccentricity seems to be the best measure.
ecc_mean = 2.3638;
ecc_std = 1.1388;
ecc_min = ecc_mean-1.5*ecc_std;
ecc_max = ecc_mean+1.5*ecc_std;

% Cutoff for same-cell distance
cell_cut = 7; % Number of pixels (5-10 is generally ok)
div_cut = 30; % Daughter cells are usually less than 60 pixels

% Calculate phylogeny?
phylogeny = true;

% Number of frames
nframes = length(frame);

% Create intermediate 'pombe' tree for analyzing phylogeny more easily
   
for k=1:nframes
    for j=1:length(frame(k).object)
        frame(k).object(j).frame_num = k;
    end
end





        
% Now do phylogeny analysis based on cuts between frames


    
if phylogeny
    
    for k=1:nframes-1 % Do each frame
        for i=1:numel(frame(k).object) % Cell to be tested
            disp(strcat('Processing frame:',num2str(k));

            diffcent = [];
            index = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                for j=1:numel(frame(k+1).object)
                    try
                        cent1 = [frame(k).object(i).Xcent_cont frame(k).object(i).Ycent_cont];
                        cent2 = [frame(k+1).object(j).Xcent_cont frame(k+1).object(j).Ycent_cont];

                        diffcent(j) = norm(cent1 - cent2);
                        index(j) = j;
                    end


                end
            end
            % Sort differences with indices
            difference = [diffcent' index'];
            [D, I] = sort(difference(:,1));


            % This is the analysis vector
            sortdiff = difference(I,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Is it the same cell?
            if sortdiff(1,1)<cell_cut
                frame(k).object(i).relation = sortdiff(1,2);

            else
                frame(k).object(i).relation = [];
            

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This provides empty relation if there are no daughter cells or 
            % same-cells



        end
    end


        
end
        
        

        
        
        
        
        
        
        
        