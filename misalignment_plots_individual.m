%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-05-13
% All rights reserved.
%========================================


% filename = 'E:\stimDustData\backscatter_and_scope_data\C11.mat'
% load(filename)
% for m = pulsesToUse
%     figure(1); hold off
%     plot(dataStore(m).scope.xData(2,:), dataStore(m).scope.yData(2,:))
%     ylim([-1 4]); 
%     drawnow
%     ROI = dataStore(m).scope.yData(2,1:5000);
%     misalignment(m) = dataStore(m).stepperCurrOffsetMM;
%     pkpk(m) = max(ROI) - min(ROI);
%     figure(2); hold on
%     plot(misalignment, pkpk); 
%     ylim([0 5]); 
%     drawnow
%     m
% end
clearvars -except dataStore filename
clear fn_acousticMeasurementsHydrophone
clear fn_format_and_save_figure
clear fn_gaussfiltcoef
clear fn_getSinusoidAmplitudeIQ
clear fn_getSinusoidAmplitudeSmooth
clear fn_plot_fft
clear fn_scopeFilter
clear fn_stimBackscatterProcess
clear fn_stimdustBenchtopDataPrep

PLOTTING = 0;
figuresFoldername = 'figures_out_bmod_norm';
figuresFoldername = '2018-06-28_tissuedistance';
figuresFoldername = '2018-08-06_bmod';
figuresFoldername = '2019-02-17_new_txdr_depth';
figuresFoldername = '2019-04-25_new_test';
figuresFoldername = '2019-05-12_new_transducer_characterization_large';
figuresFoldername = '2019-05-12_new_transducer_characterization_small';
figuresFoldername = '2019-05-12_vtx';

% SWEEPTYPE = 'angle';
% SWEEPTYPE = 'xtrans';
% SWEEPTYPE = 'ztrans';
% SWEEPTYPE = 'bmod';
% SWEEPTYPE = 'hydrophone';
% SWEEPTYPE = 'na';
% SWEEPTYPE = 'transHydrophone'
% SWEEPTYPE = 'xtransHydrophone'
% SWEEPTYPE = 'ztransHydrophoneLarge';
% SWEEPTYPE = 'ztransHydrophoneSmall';
% SWEEPTYPE = 'vTxHydrophoneLarge';
SWEEPTYPE = 'vTxHydrophoneSmall';
DOALLSWEEPS = 0;

scopeChannels.pz1 = 2;
scopeChannels.pz2 = 1;
scopeChannels.vdd = 3;
scopeChannels.stim = 4;


%======== prep
addpath('C:\Users\david\dkpStore\UCBSF\R\Writing\Resources\figure_colorschemes\DrosteEffect-CubeHelix-00335f8\DrosteEffect-CubeHelix-00335f8')
myColormap = cubehelix(128, 1.72, 1.31, 1.95, 1.30, [0.15, 0.82], [0.15, 0.82]);
colormap(myColormap);

numColors = 5;
for n = 1:numColors
    cmapIndex = (floor((n - 1).*((length(myColormap)-1)./(numColors-1)))+1);
    color5{n} = myColormap(cmapIndex, :);
%     figure(31)
%     plot((1+n):(10+n), 'Color', color5{n}, 'LineWidth', 1.5); hold on
end

currColormap = parula(256);
currColormap = currColormap(floor(.2.*length(currColormap)):end, :);
colormap(currColormap);










if DOALLSWEEPS; SWEEPTYPE = 'angle'; end
%======================================================================
% ANGLE AND XTRANS 2D SWEEP
%======================================================================
if strcmp(SWEEPTYPE, 'angle')
    filenames = {'C18', 55
                 'C19', 45
                 'C20', 38
                 'C21', 34
                 'C22', 30
                 'C23', 24
                 'C24', 17
                 'C25', 11
                 'C26', 6
                 'C27', -1
                 'C28', -8
                 'C29', -13
                 'C30', -20
                 'C31', -26
                 'C32', -32
                 'C33', -39
                 'C34', -43
                 'C35', -51
                 'C37', -58};
     % correct for angle centering
     for n = 1:length(filenames)
         filenames{n,2} = filenames{n,2} - 6
     end
     
     
     [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 1:5000);

     
    % correct extra elements in Misalignment and AngleOffset
    Misalignment = Misalignment(:,1:37);
    AngleOffset = AngleOffset(:,1:37);
    Pkpk = Pkpk(:,1:37);
    Ampvdd = Ampvdd(:,1:37);
    Stimval = Stimval(:,1:37);
    Back

%     figure
%     surf(Misalignment, AngleOffset, Pkpk)
%     figure
%     surf(Misalignment, AngleOffset, Ampvdd)
%     figure
%     surf(Misalignment, AngleOffset, Stimval)
%     figure
%     surf(Misalignment, AngleOffset, BackscatterL1normratio)
%     
%     
%     figure(40); clf;
%     stimBit = Stimval > 1e-4;
%     surf(Misalignment, AngleOffset, stimBit.*1.00, 'LineStyle', 'none')
% %     view(0,90)
%     xlabel('Translational misalignment (mm)')
%     ylabel('Angular misalignment (degrees)')
%     title('Stim dust misalignment operational range')
%     xlim([min(min(Misalignment)), max(max(Misalignment))])
%     ylim([min(min(AngleOffset)), max(max(AngleOffset))])
%     fn_format_and_save_figure(40, [pwd filesep 'operational_range_1'], [800 800])
    
    
%     figure(41); clf; hold on
%     stimBit = Stimval > 1e-4;
%     indicesOn = find(Stimval > 1e-4);
%     indicesOff = find(Stimval <= 1e-4);
%     plot3(Misalignment(indicesOff), AngleOffset(indicesOff), zeros(size(indicesOff)), 'r.', 'MarkerSize', 15)
%     plot3(Misalignment(indicesOn), AngleOffset(indicesOn), ones(size(indicesOn)), 'b.', 'MarkerSize', 15)
%     view(0,90)
%     xlabel('Translational misalignment (mm)')
%     ylabel('Angular misalignment (degrees)')
%     title('Stim dust misalignment operational range')




    currColormap = parula(256);
    currColormap = currColormap(floor(.2.*length(currColormap)):end, :);

    %======== Vpiezo pkpk 2D plot
    figure(201); clf; hold on
    [C,h] = contourf(Misalignment, AngleOffset, Pkpk, ((round(min(min(Pkpk)), 1) - .1):.4:(round(max(max(Pkpk)), 1) + .1)), 'EdgeColor', [.1 .1 .1]);
    ylim([-52 49])
    myColorbar = colorbar;
    myColorbar.Box = 'off';
    myColorbar.Label.String = 'Voltage (V)';
    plot(xlim, [0 0], 'color', [1 1 1], 'LineWidth', .15)
    plot([0 0], ylim, 'color', [1 1 1], 'LineWidth', .15)
    contourLevelsToUse = h.LevelList((length(h.LevelList) + 1 - floor(1:2:end)));
    clabel(C, h, contourLevelsToUse, 'BackgroundColor', 'none', 'Margin', 5, 'LabelSpacing', 10000);
    colormap(currColormap);
    xlabel('translational offset (mm)')
    ylabel('angular offset (mm)')
    title('piezo peak-peak voltage (V)')
    fn_format_and_save_figure(201, [pwd filesep figuresFoldername filesep 'mis_2D_Pkpk'], [700 600])


    
    %======== Vdd 2D plot
    figure(202); clf; hold on
    [C,h] = contourf(Misalignment, AngleOffset, Ampvdd, ((round(min(min(Ampvdd)), 1) - .1):.2:(round(max(max(Ampvdd)), 1) + .1)), 'EdgeColor', [.1 .1 .1]);
%     ylim([-min(abs(ylim)), min(abs(ylim))])
    ylim([-52 49])
    myColorbar = colorbar;
    myColorbar.Box = 'off';
    myColorbar.Label.String = 'Voltage (V)';
    plot(xlim, [0 0], 'color', [1 1 1], 'LineWidth', .15)
    plot([0 0], ylim, 'color', [1 1 1], 'LineWidth', .15)
%     clabel(C,h, round(h.LevelList, 1))
%     h_text = clabel(C,h, h.LevelList);
    contourLevelsToUse = h.LevelList((length(h.LevelList) + 1 - floor(2:2:end)));
    h_contourInfo = clabel(C, h, contourLevelsToUse, 'BackgroundColor', 'none', 'Margin', 5, 'LabelSpacing', 10000);
%     h_contourInfo = clabel(C, 'manual', 'BackgroundColor', [1 1 1]);
%     h_Text = findobj(h_contourInfo, 'Type', 'Text');
%     h_Line = findobj(h_contourInfo, 'Type', 'Line');
%     numericValues = [];
%     for i = 1:length(h_Text)
%         numericValueCurr = round(str2num(h_Text(i).String), 2);
%         set(h_Text(i), 'String', num2str(numericValueCurr));
%         textExtentWidthMargin = h_Text(i).Extent(3) .* .10;
%         % bottom left corner, top left corner, top right corner, bottom right corner
%         Xpatch = [h_Text(i).Extent(1) - textExtentWidthMargin, h_Text(i).Extent(1) - textExtentWidthMargin, h_Text(i).Extent(1) + h_Text(i).Extent(3) + textExtentWidthMargin, h_Text(i).Extent(1) + h_Text(i).Extent(3) + textExtentWidthMargin];
%         Ypatch = [h_Text(i).Extent(2), h_Text(i).Extent(2) + h_Text(i).Extent(4), h_Text(i).Extent(2) + h_Text(i).Extent(4), h_Text(i).Extent(2)];
%         p = patch('XData', Xpatch, 'YData', Ypatch);
%         set(p, 'FaceColor', [1 1 1], 'FaceAlpha', .6, 'EdgeColor', 'none');
%         uistack(h_Text(i), 'top')
%         set(h_Text(i), 'BackgroundColor', 'none');
%         % get rid of duplicates
% %         if any(numericValueCurr == numericValues)
% %             delete(p)
% %             delete(h_Text(i))
% %             delete(h_Line(i))
% %         end
%         numericValues(i) = numericValueCurr;
%     end
    colormap(currColormap);
    xlabel('translational offset (mm)')
    ylabel('angular offset (mm)')
    title('chip Vdd2p5 (V)')    
    fn_format_and_save_figure(202, [pwd filesep figuresFoldername filesep 'mis_2D_Vdd'], [700 600]);

    
    
    %======== stim on 2D plot
    figure(203); clf; hold on
    StimBit = Stimval > 1e-4;
    [C,h] = contourf(Misalignment, AngleOffset, StimBit, 1, 'EdgeColor', [.1 .1 .1]);
    ylim([-52 49])
    contourLevelsToUse = h.LevelList((length(h.LevelList) + 1 - floor(2:2:end)));
    clabel(C, h, contourLevelsToUse, 'BackgroundColor', 'none', 'Margin', 5, 'LabelSpacing', 10000);
    plot(xlim, [0 0], 'color', [1 1 1], 'LineWidth', .15)
    plot([0 0], ylim, 'color', [1 1 1], 'LineWidth', .15)
    plot(xlim, [-25 -25], 'color', [1 1 1], 'LineWidth', .15)
    plot(xlim, [25 25], 'color', [1 1 1], 'LineWidth', .15)
    plot([-1 -1], ylim, 'color', [1 1 1], 'LineWidth', .15)    
    plot([1 1], ylim, 'color', [1 1 1], 'LineWidth', .15)
    myColorbar = colorbar;
    myColorbar.Box = 'off';
    myColorbar.Label.String = 'Voltage (V)';
    colormap(currColormap);
    xlabel('translational offset (mm)')
    ylabel('angular offset (mm)')
    title('device able to stimulate')
    fn_format_and_save_figure(203, [pwd filesep figuresFoldername filesep 'mis_2D_stimOn'], [700 600])
    
     
end









if DOALLSWEEPS; SWEEPTYPE = 'xtrans'; end
%======================================================================
% XTRANS
%======================================================================
if strcmp(SWEEPTYPE, 'xtrans')
    filenames = {'C11', 0
                 'C12', 0};
             
    filenames = {'C66', 0};
end









if DOALLSWEEPS; SWEEPTYPE = 'ztrans'; end
%======================================================================
% ZTRANS
%======================================================================
if strcmp(SWEEPTYPE, 'ztrans')
    filenames = {'C48', 0}  % C48 also C49
    
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 1:5000);
    
    % correct for the starting height of transducer and direction of travel
    Misalignment = 68 - Misalignment;

    %======== plotting
    figure(230); clf; hold on
    h_pkpk = plot(Misalignment(1,:), Pkpk(1,:), 'LineWidth', 2)
    h_ampvdd = plot(Misalignment(1,:), Ampvdd(1,:), 'LineWidth', 2)
    StimBit = Stimval > 1e-4;
%     h_stimBit1 = plot(Misalignment(1,:), StimBit(1,:), 'LineWidth', 2)
    h_stimBit = stairs(Misalignment(1,:) - .5*(Misalignment(1,2) - Misalignment(1,1)), StimBit(1,:), 'LineWidth', 2)
    h_line1 = plot(xlim, [0 0], 'color', [1 1 1], 'LineWidth', .15)
    uistack(h_line1, 'bottom')

        
    set(h_pkpk, 'color', color5{5});  % color5 goes from dark to light
%     set(hbackscatter_2, 'color', color5{4});
    set(h_ampvdd, 'color', color5{3});
%     set(hvstim_2, 'color', color5{2});
    set(h_stimBit, 'color', color5{1});
    
    xlim([33 69])
    ylim([-.2 6.7])
    xlabel('distance external transducer to mote (mm)')
    ylabel('voltage (V)')
    title('z-translation')
    fn_format_and_save_figure(230, [pwd filesep figuresFoldername filesep 'ztrans'], [600 420])  % [width height]

end










if DOALLSWEEPS; SWEEPTYPE = 'bmod'; end
%======================================================================
% BMOD
%======================================================================
if strcmp(SWEEPTYPE, 'bmod')
    filenames = {'C8', 0
                 'C6', 1};  % name, bmod bit  C6 (1), C8 (0), C9 (1)  
%     filenames = {'C9', 1};
             
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 1:5000);
    
    bsL1R = [];
    bsL1S = [];
    stimvalue = [];
    for rowNum = 1:size(IndicesRealData,1)
        bsL1R = [bsL1R BackscatterL1normRecharge(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))] ;
        bsL1S = [bsL1S BackscatterL1normStim(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))];
        stimvalue = [stimvalue Stimval(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))];
    end

    stimbit = stimvalue > 1e-4;
    
    %======== backscatter relative difference
    beginningSkip = 0;
    endSkip = 0;
    bsRD = (bsL1S(beginningSkip+1:(end-endSkip)) - bsL1R(beginningSkip+1:(end-endSkip))) ./ bsL1R(beginningSkip+1:(end-endSkip));
    stimbit = stimbit((1+beginningSkip):(end-endSkip));
    
    popOnes = bsRD(stimbit);
    popZeros = bsRD(~stimbit);
    
    
    movingAverageWindow = 1;
    popZeros = conv(popZeros, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    popZeros = popZeros((movingAverageWindow+0):(end - movingAverageWindow));
    popOnes = conv(popOnes, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    popOnes = popOnes((movingAverageWindow+0):(end - movingAverageWindow));
    
    figure(231); hold on
    histogram(popZeros, 'BinWidth', .00335, 'EdgeColor', [0.15 0.15 0.15], 'FaceColor', color5{2}, 'FaceAlpha', 1.0, 'LineWidth', 0.3)
    histogram(popOnes, 'BinWidth', .00335, 'EdgeColor', [0.15 0.15 0.15], 'FaceColor', color5{end-1}, 'FaceAlpha', 1.0, 'LineWidth', 0.3)
    xLimUse = xlim();
    
    distributionX = ((xLimUse(1) - 1.4):.0001:(xLimUse(2) + 1.4));
    distributionRRadd = normpdf(distributionX, mean(popZeros), std(popZeros));
    distributionRSadd = normpdf(distributionX, mean(popOnes), std(popOnes));  %length(bsRRadd).*
%     distributionRR = normpdf(distributionX, mean(bsRRadd), std(bsRRadd));
%     distributionRS = normpdf(distributionX, mean(bsRSadd), std(bsRSadd));
%     distributionRR = normpdf(distributionX, 0, 1);
%     distributionRS = normpdf(distributionX, 0, 1);
    
%     plot(distributionX, distributionRRadd)
%     plot(distributionX, distributionRSadd)
    xlim([-.125, 0.055])
    title('In-vitro backscatter state differentiation')
    legend('device non-operating',' device operating')
    legend boxoff
    xlabel('backscatter relative difference')
    ylabel(sprintf('count out of N = %u', length(stimbit)))
    fn_format_and_save_figure(231, [pwd filesep figuresFoldername filesep 'histogramBenchtop_L2norm_3mov'], [900*.75 550*.75])
    
    
%     
%     distributionX = ((xLimUse(1) - 1.4):.001:(xLimUse(2) + 1.4));
%     distributionZeros = normpdf(distributionX, mean(BackscatterL1normratioVect_Zeros), std(BackscatterL1normratioVect_Zeros));
%     distributionOnes = normpdf(distributionX, mean(BackscatterL1normratioVect_Ones), std(BackscatterL1normratioVect_Ones));  %length(bsRR).*
% %     distributionRR = normpdf(distributionX, mean(bsRR), std(bsRR));
% %     distributionRS = normpdf(distributionX, mean(bsRS), std(bsRS));
% %     distributionRR = normpdf(distributionX, 0, 1);
% %     distributionRS = normpdf(distributionX, 0, 1);
%     
%     plot(distributionX, distributionZeros)
%     plot(distributionX, distributionOnes)
             


thresh = -.013899;
% falseNegative = 1 - normcdf(thresh, mean(popOnes), std(popOnes));
% falsePositive = normcdf(thresh, mean(popZeros), std(popZeros));
popZeros = popZeros(21:end-19)
thresh = -.03015;

thresh = -.0303075;
fprintf('thresh %d, falseNeg %d, falsePos %d\n', thresh, (1 - normcdf(thresh, mean(popOnes), std(popOnes))), (normcdf(thresh, mean(popZeros), std(popZeros))))
fprintf('ones mean %d, std %d; zeros mean %d, std %d\n', mean(popOnes), std(popOnes), mean(popZeros), std(popZeros))

end












%************************************************************************
%************************************************************************
%************************************************************************
% hydrophone
%************************************************************************

clear scopeChannels
scopeChannels.pz1 = 4;
scopeChannels.pz2 = 4;
scopeChannels.vdd = 4;
scopeChannels.stim = 4;

hydrophoneConversionFactor = 12774.34; % (W/m^2) per (V)  ***** NOTE ***** use fn_acousticMeasurementsHydrophone in instead of this


if DOALLSWEEPS; SWEEPTYPE = 'ztransHydrophoneLarge'; end
%======================================================================
% HYDROPHONE Z TRANS
%======================================================================
if strcmp(SWEEPTYPE, 'ztransHydrophoneLarge')
    filenames = {'C69_C70', 0}  % small | ?? volt    | z trans 0 to +24; start at 27.5 mm and increase distance   C66 for xtrans; C69 for ztrans
%     filenames = {'C73', 0}
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 12000:15000);

    % 2019-05-11 hack just for ztransHydrophone:
    if(1)
        vAmpHydrophoneAll = PzPlusAmp;
        Misalignment = Misalignment;
    else
        for rowNum = 1:size(IndicesRealData,1)
            vAmpHydrophoneAll = PzPlusAmp(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
            Misalignment = Misalignment(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
        end
    end
    
    % average values for common misalignments
    

%     vAmpHydrophoneAll = vAmpHydrophoneAll .* 3080 ./ 2525;  % correct for 25 for alignment data versus 32 for distance data
%     vAmpHydrophoneAll = vAmpHydrophoneAll.*.92;
%     vAmpHydrophoneAll = vAmpHydrophoneAll./.92;
%     vAmpHydrophoneAll = vAmpHydrophoneAll .* .9228;
    vAmpHydrophoneAll = 0.9441 .* vAmpHydrophoneAll;  % correct for longer pulses being lower amplitude than shorter pulses (60 us vs 24 us)
    
%     vAmpHydrophoneAll = 1.225178866 .* vAmpHydrophoneAll;  % scale for max ISPTA 0p3 derated of 720 mW/cm^2
    
    Distance_water = (32.6 - Misalignment) ./ 1000; % in meters
%     Distance_water = Misalignment; % in meters       % just for x-trans C66
    for n = 1:length(vAmpHydrophoneAll)
        [intensity_timeAverage(n), p_peakRare(n), MI(n)] = fn_acousticMeasurementsHydrophone(vAmpHydrophoneAll(n),  Distance_water(n));  % Distance_water(n)   (32.6 + 13.8)./1000
    end
    
    figure(234); clf; hold on
    plot(Distance_water.*1000, intensity_timeAverage ./ 10, 'LineWidth', 2, 'color', color5{2})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Derated acoustic intensity time average (mW/cm^2)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'ztransHydrophoneLarge32Intensity'], [900*.75 550*.75.*.75])
    fprintf('\nmax ispta %d', max(intensity_timeAverage ./ 10));
    clf

    plot(Distance_water.*1000, MI, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Mechanical Index')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'ztransHydrophoneLarge32MI'], [900*.75 550*.75.*.75])
    fprintf('\nmax MI %d', max(MI));
    clf

    plot(Distance_water.*1000, p_peakRare, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Peak rarefactional pressure (Pa)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'ztransHydrophoneLargep'], [900*.75 550*.75.*.75])
    fprintf('\nmax p_peakRare %d', max(p_peakRare));
    
end


%========
% 
if(0)
    vAmpHydrophoneAllOrig = vAmpHydrophoneAll;
    vAmpHydrophone = (36.8/25.6) .* vAmpHydrophoneAllOrig;  % correct for max Tx voltage relative to measure voltage; assume everything linear
    Distance_medium = ((32.6 - Misalignment)) ./ 1000; % in meters   32.6 ; the minus 2 for C71 (not used for C69 & C70)
    
    for n = 1:length(vAmpHydrophone)
        [intensity_timeAverage(n), p_peakRare(n), MI(n)] = fn_acousticMeasurementsHydrophone(vAmpHydrophone(n),  Distance_medium(n));  % Distance_water(n)   (32.6 + 13.8)./1000
    end
    
    figure(2341); clf; hold on

    plot(Distance_medium.*1000, intensity_timeAverage ./ 10, 'LineWidth', 2, 'color', color5{2})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Derated acoustic intensity time average (mW/cm^2)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    myxlim = xlim;
    xlim([30 70])
    myylim = ylim;
    ylim([0 myylim(2)])
    fn_format_and_save_figure(2341, [pwd filesep figuresFoldername filesep 'ztransHydrophoneLarge32Intensity_newmaterials28.9'], [900*.75 550*.75.*.75])
    fprintf('\nmax ispta %d', max(intensity_timeAverage ./ 10));
end





if DOALLSWEEPS; SWEEPTYPE = 'ztransHydrophoneSmall'; end
%======================================================================
% HYDROPHONE Z TRANS
%======================================================================
if strcmp(SWEEPTYPE, 'ztransHydrophoneSmall')
    filenames = {'C95', 0}  % small | ?? volt    | z trans 0 to +24; start at 27.5 mm and increase distance
    
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 12000:15000);

    % 2019-05-11 hack just for ztransHydrophone:
    if(1)
        vAmpHydrophoneAll = PzPlusAmp;
        Misalignment = Misalignment;
    else
        for rowNum = 1:size(IndicesRealData,1)
            vAmpHydrophoneAll = PzPlusAmp(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
            Misalignment = Misalignment(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
        end
    end
    
    
%     vAmpHydrophoneAll = vAmpHydrophoneAll .* 10.5 ./ 25.6; 
    vAmpHydrophoneAll = vAmpHydrophoneAll .* 24.8 ./ 25.6;  % correct for in-vivo voltage relative to tested voltage
%     vAmpHydrophoneAll = vAmpHydrophoneAll .*  (3080 ./ 2525) .* 25.6 ./ 24.8;  % correct for in-vivo voltage 32
    vAmpHydrophoneAll = vAmpHydrophoneAll .* 0.2127 / 0.21669;  % correct for long pulse amplitude (60 us) versus short pulse amplitude (24 us)
    
    Distance_water = ((-1.*Misalignment)) ./ 1000; % in meters
    Distance_actual = Distance_water + 0.01;
    for n = 1:length(vAmpHydrophoneAll)
        [intensity_timeAverage(n), p_peakRare(n), MI(n)] = fn_acousticMeasurementsHydrophone(vAmpHydrophoneAll(n), Distance_water(n));   % Distance_water(n)
    end
    
    
    
    figure(237); clf; hold on

    plot(Distance_actual.*1000, intensity_timeAverage ./ 10, 'LineWidth', 2, 'color', color5{2})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Derated acoustic intensity time average (mW/cm^2)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(237, [pwd filesep figuresFoldername filesep 'ztransHydrophoneSmall24p8Intensity'], [900*.75 550*.75.*.75])
    fprintf('\nmax ispta %d', max(intensity_timeAverage ./ 10));

    clf

    plot(Distance_actual.*1000, MI, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Mechanical Index')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(237, [pwd filesep figuresFoldername filesep 'ztransHydrophoneSmall24p8MI'], [900*.75 550*.75.*.75])
    fprintf('\nmax MI %d', max(MI));

    clf
    
    plot(Distance_actual.*1000, p_peakRare, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Peak rarefactional pressure (Pa)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(237, [pwd filesep figuresFoldername filesep 'ztransHydrophoneSmall24p8_p'], [900*.75 550*.75.*.75])
    fprintf('\nmax p_peakRare %d', max(p_peakRare));

end









if DOALLSWEEPS; SWEEPTYPE = 'xtransHydrophone'; end
%======================================================================
% HYDROPHONE X TRANS
%======================================================================
if strcmp(SWEEPTYPE, 'xtransHydrophone')
    %C66 is 32.2 vtx large xtrans; C68 is 25.6 vtx large xtrans; C97 is 32.5 vtx xtrans small; C98 is 25.6 vtx xtrans small
    filenames = {'C68', 0}  % small | ?? volt    | z trans 0 to +24; start at 27.5 mm and increase distance   C66 for xtrans; C69 for ztrans
    
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 12000:15000);

    if(1)
        vAmpHydrophoneAll = PzPlusAmp;
        Misalignment = Misalignment;
    else
        for rowNum = 1:size(IndicesRealData,1)
            vAmpHydrophoneAll = PzPlusAmp(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
            Misalignment = Misalignment(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
        end
    end
    
    
    
%     vAmpHydrophoneAll = vAmpHydrophoneAll .* 3080 ./ 2525;  % correct for 25 for alignment data versus 32 for distance data



    %== for small transducer
%     vAmpHydrophoneAll = vAmpHydrophoneAll .* 24.8 ./ 25.6;  % correct for in-vivo voltage relative to tested voltage
%     vAmpHydrophoneAll = vAmpHydrophoneAll .* 0.2127 / 0.21669;  % correct for long pulse amplitude (60 us) versus short pulse amplitude (24 us)
%     Distance_water = (19.632 ./ 1000) .* ones(size(Misalignment))  %  (32.6 - Misalignment) ./ 1000; % in meters

    %== for large transducer
    vAmpHydrophoneAll = 0.9441 .* vAmpHydrophoneAll;  % correct for longer pulses being lower amplitude than shorter pulses (60 us vs 24 us)
    Distance_water = (48 ./ 1000) .* ones(size(Misalignment))  %  (32.6 - Misalignment) ./ 1000; % in meters
    


    
    %     Distance_water = Misalignment; % in meters       % just for x-trans C66
    for n = 1:length(vAmpHydrophoneAll)
        [intensity_timeAverage(n), p_peakRare(n), MI(n)] = fn_acousticMeasurementsHydrophone(vAmpHydrophoneAll(n),  Distance_water(n));  % Distance_water(n)   (32.6 + 13.8)./1000
    end
    
    
    [M,I] = max(vAmpHydrophoneAll);
    xPosition = Misalignment - Misalignment(I);
    
    figure(235); clf; hold on

    plot(xPosition.*1000, intensity_timeAverage ./ 10, 'LineWidth', 2, 'color', color5{2})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Derated acoustic intensity time average (mW/cm^2)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(235, [pwd filesep figuresFoldername filesep 'xtransHydrophoneLarge25Intensity'], [900*.75 550*.75.*.75])
    fprintf('\nmax ispta %d', max(intensity_timeAverage ./ 10));
    clf

    plot(xPosition.*1000, MI, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Mechanical Index')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(235, [pwd filesep figuresFoldername filesep 'xtransHydrophoneLarge25MI'], [900*.75 550*.75.*.75])
    fprintf('\nmax MI %d', max(MI));
    
    clf

    plot(xPosition.*1000, p_peakRare, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Peak rarefactional pressure (Pa)')
    title('Transducer acoustic output vs z-distance')
    xlabel('Distance (mm)')
    fn_format_and_save_figure(235, [pwd filesep figuresFoldername filesep 'xtransHydrophoneLarge25p'], [900*.75 550*.75.*.75])
    fprintf('\nmax p_peakRare %d', max(p_peakRare));

    
end












if DOALLSWEEPS; SWEEPTYPE = 'transHydrophone'; end
%======================================================================
% HYDROPHONE X TRANS
%======================================================================
if strcmp(SWEEPTYPE, 'transHydrophone')
                                       % large transducer, 31.5 to 32.2 Vtx
    filenames = {
                 'C48', 0              % 1 starting 68 mm distance and moving 34.05 mm closer to be at 33.95 distance
                 'C49', 0              % 2 starting at 33.95 mm distance and moving 33 mm up to be at 66.95 mm distance
                 'C66', 0              % 3 x trans
                 'C67', 0              % 4 x trans coarse
                 'C93', 0              % 5 small transducer 32 vtx start at 9.632 mm from latex and move 10 mm closer to 0 mm from latex
                 'C94', 0              % 6 same as above
                 'C95', 0              % 7 z trans at 25.6 Vtx
                 'C96', 0              % 8 same as above
                 'C97', 0              % 9 x trans 32.5 Vtx
                 'C98', 0              % 10 x trans 25.6
                 
                 'C66', 0   % large | 32.2 Vpptx | x trans
                 'C67', 0   % large | 32.2 Vpptx | x trans coarse
                 'C68', 0   % large | 25.6 Vpptx | x trans
                 'C69', 0   % large | ?? volt    | z trans 0 to +24
                 'C70', 0   %                      z trans +24 to +28
                 'C71', 0   %                      z sweep (start at height minus 2 in z -scale)
                 
                 
                 };
             
    filenames = {
                 'C48', 0              % 1 starting 68 mm distance and moving 34.05 mm closer to be at 33.95 distance
                 'C95', 0              % 7 z trans at 25.6 Vtx
                };

    filenames = {'C48', 0};

    
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 13000:19000);
    
    for rowNum = 1:size(IndicesRealData,1)
        vAmpHydrophoneAll = PzPlusAmp(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
        Misalignment = PzPlusAmp(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2));
        
        if(rowNum == 1)
        end
            
        
        bsL1R = [bsL1R BackscatterL1normRecharge(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))] ;
        bsL1S = [bsL1S BackscatterL1normStim(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))];
        stimvalue = [stimvalue Stimval(rowNum, IndicesRealData(rowNum,1):IndicesRealData(rowNum,2))];
    end
    
    vTx = AngleOffset(:,1);
    vAmpHydrophone = hydrophoneConversionFactor.*mean(PzPlusAmp(:,2:4), 2);
    
    
    % sort
    [vTx,I] = sort(vTx);
    vAmpHydrophone = vAmpHydrophone(I);

    [intensity_timeAverage, p_peakRare, MI] = fn_acousticMeasurementsHydrophone(vAmpHydrophone, distance)
    
    figure(234); clf; hold on
    plot(vTx, vAmpHydrophone, 'LineWidth', 2, 'color', color5{2})
    xlim([0, 33])
    title('max acoustic pressure vs. V Tx (large transducer)')
    xlabel('Voltage Tx (V)')
    ylabel('Acoustic intensity (W/m^2)')
    fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'Vtx_hydrophone_small'], [900*.75 550*.75])

end









if DOALLSWEEPS; SWEEPTYPE = 'vTxHydrophoneLarge'; end
%======================================================================
% HYDROPHONE VTx SWEEP
%======================================================================
if strcmp(SWEEPTYPE, 'vTxHydrophoneLarge')
    filenames = {'C51', 5
                 'C52', 7
                 'C53', 9
                 'C54', 11
                 'C55', 13
                 'C56', 15
                 'C57', 17
                 'C58', 19
                 'C59', 21
                 'C60', 23
                 'C61', 25
                 'C62', 27
                 'C63', 29
                 'C64', 31
                 'C65', 32};  % name; Vtx
%     filenames = {'C65', 32}; 
             
             
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 1000:2000); %13000:19000);
    
    vTx = AngleOffset(:,1);
    %vAmpHydrophone = hydrophoneConversionFactor.*mean(PzPlusAmp(:,2:4), 2);
    vAmpHydrophone = hydrophoneConversionFactor.*(PzPlusAmp);
    
    % sort
    [vTx,I] = sort(vTx);
    vAmpHydrophone = vAmpHydrophone(I);

    distance = .048; % m
    [intensity_timeAverage, p_peakRare, MI] = fn_acousticMeasurementsHydrophone(vAmpHydrophone, distance);
    

    figure(234); clf; hold on
    plot(vTx, p_peakRare, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    %plot(vTx, intensity_timeAverage ./ 10, 'LineWidth', 2, 'color', color5{4})  % distance data in m, plot in mm; intensity in W/m^2 plot in mW/cm^2
    ylabel('Peak rarefactional pressure (Pa)')
    title('Max acoustic pressure vs. V Tx (25.4 mm diameter transducer)')
    xlabel('Transmit voltage (V)')
    xlim([0, 33])
%     fn_format_and_save_figure(235, [pwd filesep figuresFoldername filesep 'Vtx_hydrophone_large'], [900*.75 550*.75.*.75])
    fprintf('\nmax p_peakRare %d', max(p_peakRare));
    
    
    vq1 = interp1(vTx, intensity_timeAverage, 25.0)
    vq2 = interp1(vTx, intensity_timeAverage, 31.5)
    
    %== linear regression 
    x_regress = vTx(vTx<=21);
    x_regress = [ones(length(x_regress),1), x_regress];
    y_regress = p_peakRare(vTx<=21);
    b = mldivide(x_regress, y_regress);
    x_plot = 5:.01:31;
    y_plot = b(1) + b(2) * x_plot;
    y_normalize = y_plot((22-5)*100)*25/22;
    y_plot = y_plot./y_normalize;

    figure(234); clf; hold on
    plot(vTx, p_peakRare./y_normalize, 'LineWidth', 2, 'color', color5{2})
    plot(x_plot, y_plot, 'Color', [.6, .6, .6])
    xlim([0, 33])
    ylim([0, 1.3])
    title('normalized acoustic pressure amplitude at spatial peak versus transmit voltage (large transducer)')
    xlabel('Voltage Tx (V)')
    ylabel('normalized pressure amplitude)')
%     fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'Vtx_hydrophone_large_normalized_amp'], [900*.75 550*.75])    
end









if DOALLSWEEPS; SWEEPTYPE = 'vTxHydrophoneSmall'; end
%======================================================================
% HYDROPHONE VTx SWEEP
%======================================================================
if strcmp(SWEEPTYPE, 'vTxHydrophoneSmall')
    filenames = {'C76', 5
                 'C77', 7
                 'C78', 9
                 'C79', 11
                 'C80', 13
                 'C81', 15
                 'C82', 17
                 'C83', 19
                 'C84', 21
                 'C85', 23
                 'C86', 25
                 'C87', 27
                 'C88', 29
                 'C89', 31
                 'C90', 32
                 'C91', 25.6
                 'C92', 24.8};  % name; Vtx
%     filenames = {'C65', 32}; 
             
             
    [Misalignment, AngleOffset, PzPlusAmp, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, 13000:19000);
    
    vTx = AngleOffset(:,1);
    %vAmpHydrophone = hydrophoneConversionFactor.*mean(PzPlusAmp(:,2:4), 2);
    vAmpHydrophone = hydrophoneConversionFactor.*PzPlusAmp;
    
    
    % get rid of 25.6 and 24.8
    vTx = vTx(1:end-2);
    vAmpHydrophone = vAmpHydrophone(1:end-2);
    
    % sort
    [vTx,I] = sort(vTx);
    vAmpHydrophone = vAmpHydrophone(I);
    
    distance = .048; % m
    [intensity_timeAverage, p_peakRare, MI] = fn_acousticMeasurementsHydrophone(vAmpHydrophone, distance);
    
    %== linear regression 
    x_regress = vTx(vTx<=25);
    x_regress = [ones(length(x_regress),1), x_regress];
    y_regress = p_peakRare(vTx<=25);
    b = mldivide(x_regress, y_regress);
    x_plot = 5:.01:31;
    y_plot = b(1) + b(2) * x_plot;
    y_normalize = y_plot((22-5)*100)*25/22;
    y_plot = y_plot./y_normalize;

    figure(234); clf; hold on
    plot(vTx, p_peakRare./y_normalize, 'LineWidth', 2, 'color', color5{2})
    plot(x_plot, y_plot, 'Color', [.6, .6, .6])
    xlim([0, 33])
    ylim([0, 1.3])
    title('normalized acoustic pressure amplitude at spatial peak versus transmit voltage (small transducer)')
    xlabel('Voltage Tx (V)')
    ylabel('normalized pressure amplitude)')
    fn_format_and_save_figure(234, [pwd filesep figuresFoldername filesep 'Vtx_hydrophone_small_normalized_amp'], [900*.75 550*.75])
end




fprintf('\n');





