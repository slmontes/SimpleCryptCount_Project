%code - originally from Daniel - uses the simulation output files find the
%x and y positions at each time for all cells
%tempx and tempy are positions for all cells

function [tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,ghostx,ghosty] = collectData_4CTs(multirun)
%%% Loading names in order to get the x and y positions of each of the
%%% proliferating organoid ring cells


if nargin > 0
    ages_Name = strcat('cellages_', num2str(multirun), '.dat'); % Name of cellages file
    types_Name = strcat('results_', num2str(multirun), '.vizcelltypes'); % Name of cell types file
else
    ages_Name = strcat('cellages_1.dat'); % Name of cellages file
    types_Name = strcat('results_1.vizcelltypes'); % Name of cell types file
end
%%

%trial = num2str(i);


ages_Import = LoadNonConstantLengthData(ages_Name); % Cell Ages data load
types_Import = LoadNonConstantLengthData(types_Name); % Cell types data load (mutant, proliferating, differentiated)

% ages_Import = LoadNonConstantLengthData('TestForMatlab_2MAR20/MultipleOrganoid3CTData_0.98_0.88_4.50/cellages_1.dat');
% types_Import = LoadNonConstantLengthData('TestForMatlab_2MAR20/MultipleOrganoid3CTData_0.98_0.88_4.50/results_1.vizcelltypes');

%%% Loop over the timecourse to save to x and y positions of all cells in
%%% the organoid ring changing dynamically over the simulation
for i = 1:length(ages_Import) % Runs a loop across the time course for each simulation
    cell_counter = 1; % Counter for the number of cells in the organoid ring
    cell_counter2 = 1;
    %%% Allocate temporary vectors containing each cell from the imported
    % agesa and types files
    temp_CellAges = cell2mat(ages_Import(1,i)); % Temp vector of cell ages at time step i
    temp_Type = cell2mat(types_Import(1,i)); % Temp vector of cell ages at time step i

     % How many cells in the simulation
    cell_count = (length(temp_CellAges)-1)/4;%(length(temp_CellAges)-1)/4;

    %If there are no paneth or TA(TA) cells
    stemy(i,cell_counter) = 0; % y location of each cell at time step i
    stemx(i,cell_counter) = 0; % x location of each cell at time step i
    panethy(i,cell_counter) = 0; %Paneth cells
    panethx(i,cell_counter) = 0;
    TAy(i,cell_counter) = 0;     %Transit Amplifying cells           
    TAx(i,cell_counter) = 0;
    ECy(i,cell_counter) = 0;     %Enterocyte cells
    ECx(i,cell_counter) = 0;

    if ismember(6,temp_Type)

        for j = 1:cell_count % Loop over all cells to find the proliferating (organoid ring cells)
    % cell types:
    % (0 = paneth cells, 1 = stem cells, 2 = matrigel nodes,...
    % 3 = TA cells, 5 = ghost nodes, 6 = EC (enterocytes))
%             if j == size(temp_Type,2)
%                 pause
%             end

            tempt(j) = temp_Type(j+1); 
            if tempt(j) == 0 || tempt(j) == 1 || tempt(j) == 3 || tempt(j) == 6
                tempy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                tempx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 0
                panethy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                panethx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 1
                stemy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                stemx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 3
                TAy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                TAx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 6
                ECy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                ECx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 2
                ghosty(i,cell_counter2) = temp_CellAges(j*4); % y location of each cell at time step i
                ghostx(i,cell_counter2) = temp_CellAges((j*4)-1);
                cell_counter2 = cell_counter2 + 1;
            end
        end

      elseif ~ismember(6,temp_Type) && ismember(3,temp_Type)

         for j = 1:cell_count % Loop over all cells to find the proliferating (organoid ring cells)

            tempt(j) = temp_Type(j+1); % cell types (0 = paneth cells (stiffer cells), 1 = stem cells, 2 = matrigel cells, 3 = TA cells, 5 = ghost cells)
            if tempt(j) == 0 || tempt(j) == 1 || tempt(j) == 3
                tempy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                tempx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 0
                panethy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                panethx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 1
                stemy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                stemx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 3
                TAy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                TAx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 2
                ghosty(i,cell_counter2) = temp_CellAges(j*4); % y location of each cell at time step i
                ghostx(i,cell_counter2) = temp_CellAges((j*4)-1);
                cell_counter2 = cell_counter2 + 1;
            end
        end

    else % Only stem and paneth cells

         for j = 1:cell_count % Loop over all cells to find the proliferating (organoid ring cells)

            tempt(j) = temp_Type(j+1); % cell types (0 = paneth cells (stiffer cells), 1 = stem cells, 2 = matrigel cells, 5 = ghost cells)

            if tempt(j) == 0 || tempt(j) == 1
                tempy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                tempx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 0
                panethy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                panethx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 1
                stemy(i,cell_counter) = temp_CellAges(j*4); % y location of each cell at time step i
                stemx(i,cell_counter) = temp_CellAges((j*4)-1);
                cell_counter = cell_counter + 1;
            end
            if tempt(j) == 2
                ghosty(i,cell_counter2) = temp_CellAges(j*4); % y location of each cell at time step i
                ghostx(i,cell_counter2) = temp_CellAges((j*4)-1);
                cell_counter2 = cell_counter2 + 1;
            end
            
         end
    end

end

end
