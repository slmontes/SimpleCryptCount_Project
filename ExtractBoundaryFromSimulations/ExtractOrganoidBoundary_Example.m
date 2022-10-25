%% Extracting the simulated images from set proportions
% clear all
% close all
tic
simulationName = {'MultOrg4CT_SetProp_ChangeAtD5'};
      
cnt_name = 1;

 %Add the names of all files to be analysed
for simulationAnalysed = 1:length(simulationName)
    addpath(simulationName{simulationAnalysed});
    numOfSimulations = 0;
    fileslist = dir(simulationName{simulationAnalysed});
    for i = 1:size(fileslist,1) %finds how many simulations were run to generate data
        if  ~isempty(strfind(fileslist(i).name,'cellages'))
            numOfSimulations = numOfSimulations + 1;
        end
    end

    for jTime = 1:3  %As we had 3 experiment dates
    
        simulated_organoids = [];
    
        for i = 1:numOfSimulations   

            [tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,ghostx,ghosty] = collectData_4CTs(i); %extracts data from files into usable format
    
            timePoint = [145, 241, 337]; % Day 3 = 145   %Day 5 = 241   %Day 7 = 337 or %[48:48:337]; -> Every 24 simulated hours    

            [sim_org] = ExtractOrganoidBoundary_4CTs(ghostx,ghosty,tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,timePoint(jTime));
        
            simulated_organoids=[simulated_organoids; sim_org];
            
        end
        
        names = {'simulated_organoids_example_D3'...
                'simulated_organoids_example_D5'...
                'simulated_organoids_example_D7'};
        
        save(names{cnt_name},'simulated_organoids','-mat')

        cnt_name = cnt_name + 1;
    end
end