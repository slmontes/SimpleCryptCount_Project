function [sim_org] = ExtractOrganoidBoundary_4CTs(ghostx,ghosty,tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,timePoint)

    [epithelialCells, nonGhostsx, nonGhostsy,panethCells,stemCells,TACells, ECcells] = Organoid_Plots_4CTs(ghostx,ghosty,tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,timePoint);
            
    epithelialCells = [tempx(timePoint,:)',tempy(timePoint,:)'];
    epithelialCells(~any(epithelialCells,2),:)=[];

%     figure
%     plot(epithelialCells(:,1),epithelialCells(:,2),'*')
    
    popSize = size(epithelialCells,1);
    
    if mod(popSize,4)~=0
    popSize = popSize + (4-mod(popSize,4));
    end
    
    userConfig = struct('showProg',false,'showResult',false, 'xy',epithelialCells, 'popSize',popSize, 'numIter',1e4);
    resultStruct = tsp_nn(userConfig);
    
    testx=resultStruct.xy(resultStruct.optSolution,1)*10;
    testy=resultStruct.xy(resultStruct.optSolution,2)*10;

%     figure
%     plot(testx,testy)
    
    bw = poly2mask(testx,testy,350,350);

%     figure
%     imshow(bw)
    
    Ibw_bound = bwboundaries(bw);
    
    if size(Ibw_bound,1)~=1
        tcnt=1;
        while size(Ibw_bound,1)~=1
            userConfig = struct('showProg',true,'showResult',false, 'xy',epithelialCells, 'popSize',popSize, 'numIter',power(10,4+tcnt));
            resultStruct = tsp_nn(userConfig);
            
            testx=resultStruct.xy(resultStruct.optSolution,1)*10;
            testy=resultStruct.xy(resultStruct.optSolution,2)*10;
            
            bw = poly2mask(testx,testy,350,350);
            Ibw_bound = bwboundaries(bw);
    
            tcnt=tcnt+1;
            
            if tcnt>5
                error('Error: Counter is larger than 8, need to fix simulated_organoids vector')
            end
        end
    end
    
    sim_org= [Ibw_bound];

end

