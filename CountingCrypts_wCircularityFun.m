function [NumCrypts Circularity] = CountingCrypts_wCircularityFun (Type, MaskSet_name,fourier_harmonic_term,crypt_parameters)

    Input_min_area = crypt_parameters(1);
    Input_max_area = crypt_parameters(2);
    Input_min_arcLength = crypt_parameters(3);
    
    load(MaskSet_name);
    
    cnt =1;
    
    if isequal(Type,'In vitro') 
        disp('Invitro data')
    
        if exist ('hiddenProps','var')
            Ibw=hiddenProps.bwfinal;
            Ibw_bound = bwboundaries(Ibw);
        elseif exist ('binaryCombined','var') 
            Ibw=binaryCombined;
            Ibw_bound = bwboundaries(Ibw);
        elseif exist ('binaryImage','var') 
            Ibw=binaryImage;
            Ibw_bound = bwboundaries(Ibw);
        else
            Ibw_bound = boundary_set;
        end
    
    %     for w=1%:length(Ibw_bound)
    
            NumCrypts_perMask = [];
            Circularity_perMask = [];
    
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
    
    
                for g=1:length(midrow_section)
                    if g==length(midrow_section) 
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
    
                Circularity_perMask(end+1) = Org_circularity;
    
            end
        
            NumCrypts=NumCrypts_perMask;
            Circularity=Circularity_perMask;
        
    elseif isequal(Type,'In silico')
        disp('Insilico data')
        
        if exist((MaskSet_name))~=1
            simulations = simulated_organoids;
        else
            simulations = eval(MaskSet_name);
        end
    
        NumCrypts_InSilico_All = [];
        Circularity_InSilico_All = [];
       
        for v=1:length(simulations)
        
            NumCrypts_perBoundary = [];
            
            % we change the boundary as chain code
            [cc]= chaincode(simulations{v});
    
            %Obtain fourier (x) with new code
            n=10000;
            factor = 1;
            org = v;
            test = fourier_approx((cc.code)', fourier_harmonic_term, n, 1);  %For all In Silico days, it is recommended 
                                                                              %to use 25 as fourier harmonic term                                                                       
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
            
            if length(midrow_section) == 2
                 [Norm_crypt_area(end+1), CryptSection{1}, Curv_sum(end+1)] = CalculateSectionArea(midrow_section(1,1),midrow_section(1,2),...
                                                                                        midrow_section(1,1),midrow_section(1,2),x,y,k,N);
            else                                                                           
                for g=1:length(midrow_section)
                    if g==length(midrow_section)
                        [Norm_crypt_area(end+1), CryptSection{g}, Curv_sum(end+1)] = CalculateSectionArea(midrow_section(g,1),midrow_section(g,2),...
                                                                                            midrow_section(1,1),midrow_section(1,2),x,y,k,N);
                    else
                        [Norm_crypt_area(end+1), CryptSection{g}, Curv_sum(end+1)] = CalculateSectionArea(midrow_section(g,1),midrow_section(g,2),...
                                                                                        midrow_section(g+1,1),midrow_section(g+1,2),x,y,k,N);
                    end
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
    
            if Org_circularity>0.98
                NumCrypts_perBoundary = 0;
            else
                NumCrypts_perBoundary = sum(1*(Input_max_area > Norm_crypt_area & Norm_crypt_area >= Input_min_area & Norm_arcLength>= Input_min_arcLength));
            end
            
            NumCrypts_InSilico_All(end+1) = NumCrypts_perBoundary;
            Circularity_InSilico_All(end+1,1) = Org_circularity;
            
        end
        NumCrypts=NumCrypts_InSilico_All;
        Circularity=Circularity_InSilico_All;
    
    else
        error('Error: Type must be "In vitro" or "In silico" ')
    end
        
    end