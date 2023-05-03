%% Objective function to find optimal parameters: 
% Objective function will have the counting crypts function plus the 
% Hand counted crypts vector plus the error calculation. 
% It will output the error value, it will adjust the parameter values to
% obtain the minimum error value possible.

function y = simple_objective_manualSeg_D3(x)

Input_min_area = x(1);
Input_max_area = x(2);
Input_min_arcLength = x(3);
fourier_harmonic_term = x(4);

load('Day3_stackedim_serie_01.mat');

    for w=1 % If your file has several masks you can add them to this loop
       
        NumCrypts_perMask = [];
        Ibw=binaryCombined;
        Ibw_bound = bwboundaries(Ibw);

        for v=1:length(Ibw_bound)
    
                NumCrypts_perBoundary = [];
    
                % change the boundary as chain code
                [cc]= chaincode(Ibw_bound{v});
    
                %Obtain fourier (x) with new code
                n=1000;
                factor = 1;
                org = v;
                test = fourier_approx((cc.code)', fourier_harmonic_term, n, 1);   %NOTE: Add The fourier term as input!
                                                                                  % 9 for Day3, 15 for Day5 and 20 for Day7 
    
                xs{org}=test*100;
                x=xs{org}(:, 1);
                y= xs{org}(:, 2);
    
                polyin = polyshape(test);
                perim = perimeter(polyin);  
    
                k=LineCurvature2D(xs{org});
    
                N=LineNormals2D(xs{org});
                plotx=[xs{org}(:,1), xs{org}(:,1)+k.*N(:,1)]';
                ploty=[xs{org}(:,2), xs{org}(:,2)+k.*N(:,2)]';
                [in, on] = inpolygon(plotx(2, :),ploty(2, :),x,y);
    
                xseg_in=[];
                yseg_in=[];
                section=[];
                midrow_section = [];
                Norm_crypt_area = [];
                CryptSection = [];
                Curv_sum = [];
                ArcLength_Crypts = [];
                cnt=1;
    
                for i=1:length(in)-1
    
                    if in(i)==1
                        if (in(i)==in(i+1))==true
                            xseg_in = [xseg_in, x(i)];
                            yseg_in = [yseg_in, y(i)];
                        else
                           xseg_in = [xseg_in, x(i)];
                           yseg_in = [yseg_in, y(i)];
                           section{cnt} = [xseg_in', yseg_in'];
                           xseg_in=[];
                           yseg_in=[];
                           cnt = cnt+1;
                        end      
                    end
                end
    
    
                for nn = 1:length(section)
                    midrow_section(nn,:) = section{nn}(ceil(end/2), :);
                end
    
    
                for g=1:size(midrow_section,1)
                    if g==length(midrow_section) || size(midrow_section,1)<2
                        [Norm_crypt_area(end+1), CryptSection{g}, Curv_sum(end+1)] = CalculateSectionArea(midrow_section(g,1),midrow_section(g,2),...
                                                                                            midrow_section(1,1),midrow_section(1,2),x,y,k,N);
                    else
                        [Norm_crypt_area(end+1), CryptSection{g}, Curv_sum(end+1)] = CalculateSectionArea(midrow_section(g,1),midrow_section(g,2),...
                                                                                        midrow_section(g+1,1),midrow_section(g+1,2),x,y,k,N);
                    end
                end
    
                for gg=1:length(CryptSection)
    
                    ArcLength_Crypts(end+1) = arclength(CryptSection{gg});
    
                end
    
                Org_perim = perimeter(polyin);
                Org_area = area(polyin);
                Org_circularity = (4*pi*Org_area)/(Org_perim^2);
    
                TotalArcLength = arclength([x,y]);
                Norm_arcLength = ArcLength_Crypts./TotalArcLength;

            NumCrypts_perBoundary = sum(1*(Input_max_area > Norm_crypt_area & Norm_crypt_area >= Input_min_area & Norm_arcLength>= Input_min_arcLength));
    
            NumCrypts_perMask(end+1) = NumCrypts_perBoundary;

        end
    
    NumCrypts_All{w}=NumCrypts_perMask;
    
    end
    
HandContedCrypts = {[3,1,1,1]};

%Error calculation
    TrialMinusTest = cellfun(@minus,HandContedCrypts,NumCrypts_All,'UniformOutput',false);

    Div_TrialTest_wNaN = cellfun(@(x,y) x./y, TrialMinusTest, HandContedCrypts, 'UniformOutput',false);

    Abs_value_wNaN = cellfun(@abs,Div_TrialTest_wNaN,'UniformOutput',false);

    Pct_Error_wNaN = cellfun(@(x,y) x.*100, Abs_value_wNaN, 'UniformOutput',false);
    
    %Omit NaN values for the percentage error mean
    pctError_IndvMean = cellfun(@(m) mean(m, 'omitnan'), Pct_Error_wNaN, 'UniformOutput', false);

    pctError_TotalMean = mean(cell2mat(pctError_IndvMean));    
    
y = pctError_TotalMean;

end