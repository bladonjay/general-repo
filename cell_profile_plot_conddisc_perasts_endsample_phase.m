function cell_profile_plot_conddisc_perasts_endsample_phase(session,units,u,secbefore, secafter, nbins)
% I EDITED OUT ALL THE RED FOR INCORRECT RESPONSES, JUST UNCOMMENT LINE 297

% Added to take timestamps out of cell array.
% JHB edit, looking to s ee what epoch this samples,
% input looks like cell_profile_plot_conddisc_perasts_endsample(Session,Session.units,1,2,2,10)
% for ii = 1:numel(units); units(ii).ts = units(ii).ts{1}; end

% why the fuck does he have to grid out the plots so much
fig_cols = 4;
fig_rows = 2;
subplotwidth = .6/fig_cols;%same for hist and rast
subplotheight = .6/fig_rows;
rast_height = subplotheight*.6;
hist_height = subplotheight*.4;
leftmargin = .1;
topmargin = .15;
maxrate=[];

for sp = 1:8%8 subplots
    
    % raster 1
    if(sp==1)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 1;
        gridy  = 1;
        
        for t = 1:length(session.trial)
            if isfield(session.trial(t).X, 'ts')
                if session.trial(t).X.pos == 1
                    for  x = 1:length(session.trial(t).X.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).X.ts(x);
                        end_interval(ev_num) =session.trial(t).X.sample_end(x);
                        if isfield(session.trial(t).Y, 'ts')
                            y_before_x = session.trial(t).X.ts(x)-session.trial(t).Y.ts;
                            short_ys = session.trial(t).Y.ts(find(y_before_x < secbefore & y_before_x > 0));
                            if(short_ys)
                                short_pre_sample(ev_num) = (min(short_ys));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                    
                end
            end
        end % for t = 1:length(session.trial)
        
    elseif(sp==2)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 2;
        gridy  = 1;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).Y, 'ts')
                if session.trial(t).Y.pos == 1
                    for  y = 1:length(session.trial(t).Y.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).Y.ts(y);
                        end_interval(ev_num) =session.trial(t).Y.sample_end(y);
                        if isfield(session.trial(t).X, 'ts')
                            x_before_y = session.trial(t).Y.ts(y)-session.trial(t).X.ts;
                            short_xs = session.trial(t).X.ts(find(x_before_y < secbefore & x_before_y > 0));
                            if(short_xs)
                                short_pre_sample(ev_num) = (min(short_xs));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                    
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==3)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 3;
        gridy  = 1;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).X, 'ts')
                if session.trial(t).X.pos == 2
                    for  x = 1:length(session.trial(t).X.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).X.ts(x);
                        end_interval(ev_num) =session.trial(t).X.sample_end(x);
                        if isfield(session.trial(t).Y, 'ts')
                            y_before_x = session.trial(t).X.ts(x)-session.trial(t).Y.ts;
                            short_ys = session.trial(t).Y.ts(find(y_before_x < secbefore & y_before_x > 0));
                            if(short_ys)
                                short_pre_sample(ev_num) = (min(short_ys));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==4)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 4;
        gridy  = 1;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).Y, 'ts')
                if session.trial(t).Y.pos == 2
                    for  y = 1:length(session.trial(t).Y.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).Y.ts(y);
                        end_interval(ev_num) =session.trial(t).Y.sample_end(y);
                        if isfield(session.trial(t).X, 'ts')
                            x_before_y = session.trial(t).Y.ts(y)-session.trial(t).X.ts;
                            short_xs = session.trial(t).X.ts(find(x_before_y < secbefore & x_before_y > 0));
                            if(short_xs)
                                short_pre_sample(ev_num) = (min(short_xs));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                    
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==5)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 1;
        gridy  = 2;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).X, 'ts')
                if session.trial(t).X.pos == 3
                    for  x = 1:length(session.trial(t).X.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).X.ts(x);
                        end_interval(ev_num) =session.trial(t).X.sample_end(x);
                        if isfield(session.trial(t).Y, 'ts')
                            y_before_x = session.trial(t).X.ts(x)-session.trial(t).Y.ts;
                            short_ys = session.trial(t).Y.ts(find(y_before_x < secbefore & y_before_x > 0));
                            if(short_ys)
                                short_pre_sample(ev_num) = (min(short_ys));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==6)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 2;
        gridy  = 2;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).Y, 'ts')
                if session.trial(t).Y.pos == 3
                    for  y = 1:length(session.trial(t).Y.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).Y.ts(y);
                        end_interval(ev_num) =session.trial(t).Y.sample_end(y);
                        if isfield(session.trial(t).X, 'ts')
                            x_before_y = session.trial(t).Y.ts(y)-session.trial(t).X.ts;
                            short_xs = session.trial(t).X.ts(find(x_before_y < secbefore & x_before_y > 0));
                            if(short_xs)
                                short_pre_sample(ev_num) = (min(short_xs));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==7)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 3;
        gridy  = 2;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).X, 'ts')
                if session.trial(t).X.pos == 4
                    for  x = 1:length(session.trial(t).X.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).X.ts(x);
                        end_interval(ev_num) =session.trial(t).X.sample_end(x);
                        if isfield(session.trial(t).Y, 'ts')
                            y_before_x = session.trial(t).X.ts(x)-session.trial(t).Y.ts;
                            short_ys = session.trial(t).Y.ts(find(y_before_x < secbefore & y_before_x > 0));
                            if(short_ys)
                                short_pre_sample(ev_num) =(min(short_ys));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                end
            end
        end % for t = 1:length(session.trial)
    elseif(sp==8)
        events = [];
        end_interval = [];
        short_pre_sample = [];
        correct = [];
        ev_num = 0;
        gridx  = 4;
        gridy  = 2;
        for t = 1:length(session.trial)
            if isfield(session.trial(t).Y, 'ts')
                if session.trial(t).Y.pos == 4
                    for  y = 1:length(session.trial(t).Y.ts)
                        ev_num = ev_num +1;
                        correct(ev_num)  = session.trial(t).correct;
                        events(ev_num) = session.trial(t).Y.ts(y);
                        end_interval(ev_num) =session.trial(t).Y.sample_end(y);
                        if isfield(session.trial(t).X, 'ts')
                            x_before_y = session.trial(t).Y.ts(y)-session.trial(t).X.ts;
                            short_xs = session.trial(t).X.ts(find(x_before_y < secbefore & x_before_y > 0));
                            if(short_xs)
                                short_pre_sample(ev_num) = (min(short_xs));
                            else
                                short_pre_sample(ev_num) = 0;
                            end
                        else
                            short_pre_sample(ev_num) = 0;
                        end
                    end
                    
                end
            end
        end % for t = 1:length(session.trial)
    end
    
    
    %then determine the actual coords
    if gridx ==3
        leftside = .5;
    elseif gridx == 4
        leftside = .5 + subplotwidth + .044;
    elseif gridx==2
        leftside = .005+leftmargin + ((gridx-1)*subplotwidth) + (gridx-1)*(.087/fig_cols);
    else
        leftside = leftmargin + ((gridx-1)*subplotwidth) + (gridx-1)*(.087/fig_cols);
    end
    hist_bottom = 1 - topmargin - gridy*subplotheight - (gridy-1)*(.2/fig_rows);
    rast_bottom = hist_bottom + hist_height;
    
    
    %%%%%%%%%%%%%%%%%%%
    %%plot the raster%%
    %%%%%%%%%%%%%%%%%%%
    % added in late to force out coloring
    correct=ones(1,length(correct));
    
    % heres where i add phase
    rast_ax(u,sp) = subplot('Position', [leftside rast_bottom subplotwidth rast_height]);
    if isfield(units(u),'phases')
        phase=units(u).phases;
        pe_raster_endsample_phase(units(u).ts,events,secbefore,secafter,end_interval,short_pre_sample,phase);
    else
        pe_raster_endsample(units(u).ts,events,secbefore,secafter,end_interval,short_pre_sample,correct);
    end
        
    set(gca,'FontName', 'Arial','Fontsize',12);
    set(gca,'XTickLabel',[]);%don't print numbers on the raster x-axis
    set(gca,'YTick', 0:5:length(events),'YDir','reverse')
    
    %      set(gca,'YTickLabel',1:20);%don't print numbers on the raster y-axis
    if (gridx == 1 | gridx == 3)
        ylabel('Samples','FontName', 'Arial', 'FontSize', 10);
        %         else
        %             set(gca,'YTickLabel',[]);%don't print Hz unless it's the left column
    end
    %set(gca,'Xcolor',[1 1 1]);%make the x tics white (not visible)
    %put the unit number on the graph
    unitnumtx = sprintf('%d',u);
    %     text(-secbefore+(secbefore+secafter)/(nbins*2),length(events)*.9,unitnumtx,'Color', [0 1 0]);%in green
    
    %%%%%%%%%%%%%%%%%%
    %%plot histogram%%
    %%%%%%%%%%%%%%%%%%
    
    hist_ax(u,sp) = subplot('Position', [leftside (hist_bottom-.002)  subplotwidth hist_height]);
    [~,~,tempmaxrate]=pe_th(units(u).ts,events,secbefore,secafter,nbins);
    maxrate =[maxrate tempmaxrate];
    %                     if (gridy == fig_rows)
    %                            xlabel('Time (seconds)', 'FontSize', 12);
    %                     else
    %                         set(gca,'XTickLabel',[]);%don't print second unless it's the bottom row
    %                     end
    set(gca,'FontName', 'Arial','Fontsize',12);
    %set(gca,'Xtick',-(secbefore):.5:secafter);
    if (gridx == 1 | gridx == 3)
        ylabel('Hz','FontName', 'Arial', 'FontSize', 10);
    else
        set(gca,'YTickLabel',[]);%don't print Hz unless it's the left column
    end
    
    %make all the histograms for one unit have the same Y scale
    %yrange(sp,1:2) = get(hist_ax(u,sp), 'ylim');
    
    
    % annotate the outside plots
    if (sp==1)
        p = text(-.37,.7,'','FontSize', 20,'FontName', 'Arial', 'Units', 'normalized'); %coded out Position 1 label in ''
        set(p,'rotation',90);
        %         myline = annotation('line',[.15 .15],[0.509 .93]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
        %     elseif (sp==2)
        %         myline = annotation('line',[.3725 .3725],[0.509 .93]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
        %     elseif (sp==3)
        %         myline = annotation('line',[.63 .63],[0.509 .93]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
        %     elseif (sp==4)
        %         myline = annotation('line',[.852 .852],[0.509 .93]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
    elseif (sp==5)
        p = text(-.37,.7,'','FontSize', 20,'FontName', 'Arial', 'Units', 'normalized');  %coded out Position 2 label in ''
        set(p,'rotation',90);
        %         myline = annotation('line',[.15 .15],[0.05 .47]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
        %     elseif (sp==6)
        %         myline = annotation('line',[.3725 .3725],[0.05 .47]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
        %     elseif (sp==7)
        %         myline = annotation('line',[.63 .63],[0.05 .47]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
    elseif (sp==8)
        linkaxes(hist_ax(u,:),'y');%this way, changing one will change them all
        set(hist_ax(u,1),'FontName', 'Arial','ylim',[0 ceil(max(maxrate))+1]); % Temp edit to manually adjust Y axis scale CK 5/25/13 original line below
        %         set(hist_ax(u,1),'FontName', 'Arial','ylim',[0 max(yrange(:,2))]);
        %         myline = annotation('line',[.852 .852],[0.05 .47]);
        %         set(myline, 'Color', [0 0 0], 'LineWidth', 2);
    end
  
end%end sp1:8
text(-3.1,-.4,'Time (s)','FontName', 'Arial','FontSize', 15, 'Units', 'normalized');
text(-.29,-.4,'Time (s)','FontName', 'Arial','FontSize', 15, 'Units', 'normalized');
text(.3,11.5,num2str(u),'FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(-3.2,6.7,'Context A','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(-.45,6.7,'Context B','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(-3.7,6.3,'Item X','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(-2.5,6.3,'Item Y','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(-.95,6.3,'Item X','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
text(.27,6.3,'Item Y','FontName', 'Arial','FontSize', 20, 'Units', 'normalized');
if isfield(session,'name')
text(0,-1,[session.name ' unit ' num2str(u)],'FontName', 'Arial','FontSize', 12, 'Units', 'normalized');
end
%
set(gcf,'OuterPosition',[80, 400, 800, 600]);




end %end function
