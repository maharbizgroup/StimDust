function [] = fn_format_and_save_figure(figHandle, filename, figSizeInput)
%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2018-06-26
% All rights reserved.
%========================================

    %  figSize = NaN means don't touch size
    %  figSize = 0 means use standard
    %  figSize = [width height] means use a specific value

    QUICKSAVE = 0;
    % figure settings
%     figSize = [1000 550];  [width height]
    figFontSize = 16;
    figTitleSize = 4;
    figAxislabelSize = 4;
    figTicklabelSize = 4;
    
    pauseTime = .0001;

    figure(figHandle)
    allAxesInFigure = findall(figHandle,'type','axes');

    %======== set size and position
    if isnan(figSizeInput(1))
        incommingPos = get(figHandle, 'pos');
        figSize = incommingPos(3:4);
    elseif figSizeInput(1) == 0
        figSize = [900 550.*length(allAxesInFigure)]; % [width height]
    else
        figSize = figSizeInput;
    end
    monitorPos = get(groot, 'MonitorPositions');
    monitorPixels = (monitorPos(:,4) - monitorPos(:,2)).*(monitorPos(:,3) - monitorPos(:,1));
    [~,monitorToUse] = max(monitorPixels);  % use largest monitor
    figPos = [monitorPos(monitorToUse, 1) + 40, monitorPos(monitorToUse, 4) - figSize(2) - 80, figSize(1), figSize(2)];
    set(figHandle, 'pos', figPos); 

    
    for n = 1:length(allAxesInFigure)
        currAxis = allAxesInFigure(n);

        % Adjust axes properties
        set(currAxis, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.025 .025], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XColor', [.35 .35 .35], 'YColor', [.35 .35 .35],  ...
            'LineWidth', 2)
        set(currAxis, 'Color', 'None');
        set(currAxis, 'FontSize', figFontSize);

        hTitle=get(currAxis,'title');
        set(hTitle, 'color', [0.1 0.1 0.1]);

        hLegend = findobj(currAxis, 'Type', 'Legend');
        set(hLegend, 'color', [0.35 0.35 0.35]);

    end

%     myFigPos = [40, 40, 40+figSize(1), 40+figSize(2)]
%     set(figHandle, 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); 
%     box off; 
    set(figHandle, 'Color', 'None');
    set(figHandle, 'Renderer', 'painters');
    figure(figHandle); drawnow; pause(pauseTime);

    if (~QUICKSAVE)
        % save
        figure(figHandle); drawnow; pause(pauseTime);
        saveas(figHandle, [filename '.eps'], 'epsc');
        figure(figHandle); drawnow; pause(pauseTime);
        saveas(figHandle, [filename '.pdf'], 'pdf');
        figure(figHandle); drawnow; pause(pauseTime);
        saveas(figHandle, [filename '.svg'], 'svg');
        figure(figHandle); drawnow; pause(pauseTime);
    %     saveas(figHandle, [filename '.tif'], 'tiffn');
        figure(figHandle); drawnow; pause(pauseTime);

    warning('off', 'MATLAB:print:FigureTooLargeForPage')
    %     set(figHandle,'Units','Inches');
    %     posPdf = get(figHandle,'Position');
    %     set(figHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[posPdf(3), posPdf(4)])
    %     figure(figHandle); drawnow; pause(pauseTime);
        print(figHandle,[filename '.pdf'],'-dpdf','-r0')
    %     figure(figHandle); drawnow; pause(pauseTime);
    end

    saveas(figHandle, [filename '.png'], 'png');
    figure(figHandle); drawnow; pause(pauseTime);

    % set figure background to white again
    set(figHandle, 'Color', 'white');
    for n = 1:length(allAxesInFigure)
        set(allAxesInFigure(n), 'Color', 'white');
    end
    
    if(~QUICKSAVE)
        figure(figHandle); drawnow; pause(pauseTime);
        saveas(figHandle, [filename '.jpg'], 'jpeg');
        figure(figHandle); drawnow; pause(pauseTime);
        savefig(figHandle, [filename '.fig'], 'compact')

    end    
end





%{

set(0,'defaultAxesFontName', '<fontname>')
set(0,'defaultTextFontName', '<fontname>')




xlim([t1*dt*1e3 t2*dt*1e3]);
ylim([-.2 3]);
set(gca,'fontsize',16);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.025 .025], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 1)

set(figHandle, 'Color', 'None');





semilogx(Istim,(Pstim)*100,'d-','linewidth',3,'markerfacecolor','w', ...
    'markersize',8,'color',[.5 .5 .5],'markeredgecolor','k');
% xlabel('I_{Stim,AVG} (\muA)');
% ylabel('Chip Efficiency %');
set(gca,'xtick',[.1 1 10 100]);
set(gca,'fontsize',20);
% set(gca,'fontname','georgia');

xlim([.06 200]);
%}
