function leaky_integrate_fire
%
% GUI for CN210/NE204 - Introduction to Computational Models of Brain and Behavior
% Leaky Intrgate and Fire Model
%
% Trouble printing?
%  - use 'Save As...' from the 'File' menu to save the figure as a picture (jpg) file
%  or
%  - use 'Copy Figure' from the 'Edit' menu to paste the figure into a text editor
%

set(0,'DefaultAxesFontSize',10);

%
mycolors.Iin = [0.6 0 0.8];
mycolors.Ires = [0.8 0 0];
mycolors.Icap = [0 0 0.8];
mycolors.gray = [0.6 0.6 0.6];

   %  Create and then hide the GUI as it is being constructed.
   f = figure('Visible','off','Position',[10,10,675,450], 'Color', [.7 .7 .7], 'Units', 'pixel');
   set(f, 'Name', 'CN210/NE204: Leaky Integrate and Fire Lab','NumberTitle','off', 'Resize', 'off');
   set(f, 'DefaultUIPanelFontName', 'Helvetica', 'DefaultUIPanelFontSize', 10);
   set(f, 'DefaultUIControlFontName', 'Helvetica', 'DefaultUIControlFontSize', 10);
   set(f, 'DefaultTextFontName', 'Helvetica', 'DefaultTextFontSize', 10)

   %  Circuit Diagram
   pcircuit = uipanel('Title', 'Circuit diagram', 'Units', 'pixels', 'Position', [10 290 300 150]);
   axCircuit = axes('Parent', pcircuit, 'units', 'pixels', 'Position', [0 0 300 150]);
   axis(axCircuit, [0 300 0 150]);
   set(axCircuit, 'Visible', 'off');
   hold(axCircuit, 'on');
   drawCircuit();

   % Slider to set speed
   pslider = uipanel('Title', 'Simulation speed', 'Units', 'pixels', 'Position', [320 10 345 45]);
   myslider = uicontrol('Parent', pslider, 'Style', 'slider', 'Min', 0, 'Max', 1, 'Value', 0.5, 'SliderStep', [0.05 0.05], 'Units', 'normalized', 'Position', [0.15 0.25 0.7 0.5], 'Callback', {@mycallbacks});
   uicontrol('Parent', pslider, 'Style', 'text', 'String', 'slow', 'Units', 'normalized', 'Position', [0 0.25 0.15 0.5]);
   uicontrol('Parent', pslider, 'Style', 'text', 'String', 'fast', 'Units', 'normalized', 'Position', [0.85 0.25 0.15 0.5]);
   
   %  Control Panel
   pcontrols = uipanel('Title', 'Controls', 'Units', 'pixels', 'Position', [10 10 300 270]);
   
   % resistor  set
   pres = uipanel('Parent', pcontrols, 'Title', 'resistor', 'TitlePosition', 'lefttop', 'Units', 'pixels', 'Position', [4 190 144 60]);
   uicontrol('Parent', pres, 'Style', 'text', 'String', 'R=', 'Units', 'pixels', 'Position', [1 10 19 20]);
   uicontrol('Parent', pres, 'Style', 'text', 'String', 'MegaOhms', 'Units', 'pixels', 'Position', [80 10 60 20], 'HorizontalAlignment', 'left');
   editres = uicontrol('Parent', pres, 'Style', 'edit', 'String', '0', 'Units', 'pixels', 'Position', [40 4 40 40], 'Callback', {@mycallbacks});
   resplus = uicontrol('Parent', pres, 'Style', 'pushbutton', 'String', '+', 'Units', 'pixels', 'Position', [20, 24, 20, 15], 'Callback', {@mycallbacks});
   resminus = uicontrol('Parent', pres, 'Style', 'pushbutton', 'String', '-', 'Units', 'pixels', 'Position', [20, 9, 20, 15], 'Callback', {@mycallbacks});

   % capacitance set
   pcap = uipanel('Parent', pcontrols, 'Title', 'capacitor', 'TitlePosition', 'lefttop', 'Units', 'pixels', 'Position', [152 190 144 60]);
   uicontrol('Parent', pcap, 'Style', 'text', 'String', 'C=', 'Units', 'pixels', 'Position', [10 0 0 0]+[1 10 19 20]);
   uicontrol('Parent', pcap, 'Style', 'text', 'String', 'nanoF', 'Units', 'pixels', 'Position', [10 0 0 0]+[80 10 50 20], 'HorizontalAlignment', 'left');
   editcap = uicontrol('Parent', pcap, 'Style', 'edit', 'String', '0', 'Units', 'pixels', 'Position', [10 0 0 0]+[40 4 40 40], 'Callback', {@mycallbacks});
   capplus = uicontrol('Parent', pcap, 'Style', 'pushbutton', 'String', '+', 'Units', 'pixels', 'Position', [10 0 0 0]+[20 24 20 15], 'Callback', {@mycallbacks});
   capminus = uicontrol('Parent', pcap, 'Style', 'pushbutton', 'String', '-', 'Units', 'pixels', 'Position', [10 0 0 0]+[20 9 20 15], 'Callback', {@mycallbacks});
   
   % Input current options
   bgcurrent = uipanel('Parent', pcontrols, 'Title', 'input current I(t)', 'Units', 'pixels', 'Position', [4 4 292 176], 'ForegroundColor', mycolors.Iin);

   % for constant current
   bcurrent1 = uicontrol( 'Parent', bgcurrent, 'Style', 'radiobutton', 'String', ' ', 'Units', 'pixels', 'Posit', [10 110 20 45], 'Callback', {@mycallbacks});
   tcurrent1 = uicontrol( 'Parent', bgcurrent, 'Style', 'text', 'String', 'constant current:', 'Units', 'pixels', 'Posit', [35 120 110 20], 'HorizontalAlignment', 'left');
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'I=', 'Units', 'pixels', 'Position', [150 110 0 0]+[1 10 19 20]);
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'mA', 'Units', 'pixels', 'Position', [150 110 0 0]+[80 10 40 20],'HorizontalAlignment', 'left');
   editcur1 = uicontrol('Parent', bgcurrent, 'Style', 'edit', 'String', '0', 'Units', 'pixels', 'Position', [150 110 0 0]+[40 4 40 40], 'Callback', {@mycallbacks});
   cur1plus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '+', 'Units', 'pixels', 'Position', [150 110 0 0]+[20 24 20 15], 'Callback', {@mycallbacks});
   cur1minus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '-', 'Units', 'pixels', 'Position', [150 110 0 0]+[20 9 20 15], 'Callback', {@mycallbacks});

   % or pulsed current (only one can be activated at a time)
   bcurrent2 = uicontrol( 'Parent', bgcurrent, 'Style', 'radiobutton', 'String', ' ', 'Units', 'pixels', 'Posit', [10 60 20 45], 'Callback', {@mycallbacks});
   tcurrent2 = uicontrol( 'Parent', bgcurrent, 'Style', 'text', 'String', 'current pulse:', 'Units', 'pixels', 'Posit', [35 70 110 20], 'HorizontalAlignment', 'left');
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'pulse height', 'Units', 'pixels', 'Position', [150 0 0 0]+[0 45 100 20], 'HorizontalAlignment', 'left');
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'I=', 'Units', 'pixels', 'Position', [150 0 0 0]+[1 10 19 20]);
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'mA', 'Units', 'pixels', 'Position', [150 0 0 0]+[80 10 40 20],'HorizontalAlignment', 'left');
   editcur2 = uicontrol('Parent', bgcurrent, 'Style', 'edit', 'String', '0', 'Units', 'pixels', 'Position', [150 0 0 0]+[40 4 40 40], 'Callback', {@mycallbacks});
   cur2plus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '+', 'Units', 'pixels', 'Position', [150 0 0 0]+[20 24 20 15], 'Callback', {@mycallbacks});
   cur2minus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '-', 'Units', 'pixels', 'Position', [150 0 0 0]+[20 9 20 15], 'Callback', {@mycallbacks});
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'pulse width', 'Units', 'pixels', 'Position', [20 0 0 0]+[0 45 100 20], 'HorizontalAlignment', 'left');
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'T=', 'Units', 'pixels', 'Position', [20 0 0 0]+[1 10 19 20]);
   uicontrol('Parent', bgcurrent, 'Style', 'text', 'String', 'ms', 'Units', 'pixels', 'Position', [20 0 0 0]+[80 10 40 20],'HorizontalAlignment', 'left');
   edittau2 = uicontrol('Parent', bgcurrent, 'Style', 'edit', 'String', '0', 'Units', 'pixels', 'Position', [20 0 0 0]+[40 4 40 40], 'Callback', {@mycallbacks});
   tau2plus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '+', 'Units', 'pixels', 'Position', [20 0 0 0]+[20 24 20 15], 'Callback', {@mycallbacks});
   tau2minus = uicontrol('Parent', bgcurrent, 'Style', 'pushbutton', 'String', '-', 'Units', 'pixels', 'Position', [20 0 0 0]+[20 9 20 15], 'Callback', {@mycallbacks});
   
   % output panel, with graphs
   poutput = uipanel('Title', 'Dynamics', 'Units', 'pixels', 'Position', [320 65 345 375]);
   
   % axTextout = axes('Parent', poutput, 'units', 'normalized', 'Position', [0 0 1 1]);
   % axis(axTextout, [0 1 0 1]);
   % set(axTextout, 'Visible', 'off');
   % hold(axTextout, 'on');

%    text('Parent', axTextout, 'String', 'time [ms]', 'Position', [300 4], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'k');
%    text('Parent', axTextout, 'String', 'I', 'Position', [10 30+40], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Iin);
%    text('Parent', axTextout, 'String', 'input', 'Position', [10 30+40], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Iin);
%    text('Parent', axTextout, 'String', ' [mA]', 'Position', [2 30+15], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k');
%    text('Parent', axTextout, 'String', 'I', 'Position', [10 125+55], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Ires);
%    text('Parent', axTextout, 'String', 'res', 'Position', [10 125+55], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Ires);
%    text('Parent', axTextout, 'String', 'I', 'Position', [10 125+30], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Icap);
%    text('Parent', axTextout, 'String', 'cap', 'Position', [10 125+30], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Icap);
%    text('Parent', axTextout, 'String', ' [mA]', 'Position', [2 125+5], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k');
%    text('Parent', axTextout, 'String', 'V', 'Position', [10 220+80], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k');
%    text('Parent', axTextout, 'String', ' [mV]', 'Position', [2 220+55], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k');
%    set(axTextout, 'Visible', 'on');

   
   
   axV = axes('Parent',poutput,'units','pixels', 'Position',[50 220 285 135], 'Box', 'on');
    axis(axV, [-2 42 -93 -47]); %-13 33]);
    set(axV, 'XTick', [0,10,20,30,40])%, 'XTickLabel', '0||||40')
    set(axV, 'YTick', [-90,-80,-70,-60,-50])%, 'YTickLabel', '|-80|||-50')%[0,30], 'YTickLabel', '0|30')
    set(axV, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(axV, 'on');
    lineV = plot(axV, 0, 0, '-', 'Color', 'k');
    set( lineV, 'XData', []); 
    set( lineV, 'YData', []);
%    axV.YLabel.String='Voltage';
   axIout = axes('Parent',poutput,'units','pixels', 'Position',[50 125 285 75], 'Box', 'on');
    axis(axIout, [-2 42 -12 12]);
    set(axIout, 'XTick', [0,10,20,30,40])%, 'XTickLabel', '0||||40')
    %set(axIout, 'YTick', [-10,10])%, 'YTickLabel', '-10|10')
    set(axIout, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(axIout, 'on');
    % lineIcap and res are the plot lines for the results
    lineIcap = plot(axIout, 0,0,'-', 'Color', mycolors.Icap);
    lineIres = plot(axIout, 0,0,'-', 'Color', mycolors.Ires);
%    axIout.YLabel.String='Current flow';
   axIin = axes('Parent',poutput,'units','pixels', 'Position',[50 30 285 75], 'Box', 'on');
    axis(axIin, [-2 42 -12 12]);
    set(axIin, 'XTick', [0,10,20,30,40])%, 'XTickLabel', '0||||40')
    set(axIin, 'YTick', [-10,10])%, 'YTickLabel', '-10|10')
    set(axIin, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(axIin, 'on');
    lineIin = plot(axIin, 0,0,'-', 'Color', mycolors.Iin);
%    axIin.YLabel.String='Input Current';

% set the initial values
vars.plotting=0;
vars.resistance = 10;
vars.resmin = 5;
vars.resmax = 100;
vars.dres = 0.1;
vars.capacitance = 0.6;
vars.capmin = 0.1;
vars.capmax = 1.2;
vars.dcap = 0.05;
vars.T = 0;
vars.Tmax = 40;
vars.dT = 0.053;
vars.V = 0;
vars.Vthreshold = -50; %30;
vars.Vreset = -80;
vars.Vrest = -70;
vars.Iinput = 0;
vars.currentswitch = 1;
vars.current1 = 0;
vars.cur1min = -2;
vars.cur1max = 10;
vars.dcur1 = 0.1;
vars.current2 = 5;
vars.cur2min = -1;
vars.cur2max = 10;
vars.dcur2 = 0.1;
vars.tau2 = 3;
vars.tau2min = 0.1;
vars.tau2max = vars.Tmax;
vars.dtau2 = 0.1;
vars.cacheT = [];
vars.cacheV = [];
vars.cacheIres = [];
vars.cacheIcap = [];
vars.cacheNmax = 10; % the initial value is changed later
vars.axImaxN = 10;
vars.axIminN = -10;
vars.axImax = 10.5;
vars.axImin = -10.5;

 % initialize the buttons, etc
 set( editres, 'String', num2str(vars.resistance));
 set( editcap, 'String', num2str(vars.capacitance));
 set( bcurrent1, 'Value', 1);
 set( bcurrent2, 'Value', 0);
 set( editcur1, 'String', num2str(vars.current1));      
 set( editcur2, 'String', num2str(vars.current2));
 set( edittau2, 'String', num2str(vars.tau2));
 set( myslider, 'Value', 0); updateslider;
 % a few other things
 plot(axV, [-2 42], [0 0], '--', 'Color', mycolors.gray);
 plot(axV, [-2 42], [0 0]+vars.Vthreshold, '--', 'Color', mycolors.gray);
 plot(axV, [-2 42], [0 0]+vars.Vreset, '--', 'Color', mycolors.gray);
 plot(axIin, [-2 42], [0 0], '--', 'Color', mycolors.gray);
 plot(axIout, [-2 42], [0 0], '--', 'Color', mycolors.gray);
 
 
 
   % Make the GUI visible.
   movegui(f,'center');
   set(f,'Visible','on');
   resetgraphs; 

   function mycallbacks(source,eventdata) 
       switch source
           case myslider
               updateslider; return;
           case editres
               updateres;
           case resplus
               set(editres, 'String', num2str(vars.resistance + vars.dres));
               updateres;
           case resminus
               set(editres, 'String', num2str(vars.resistance - vars.dres));
               updateres;
           case editcap
               updatecap;
           case capplus
               set(editcap, 'String', num2str(vars.capacitance + vars.dcap));
               updatecap;
           case capminus
               set(editcap, 'String', num2str(vars.capacitance - vars.dcap));
               updatecap;
           case bcurrent1
               setcurrentfocus( 1);
           case bcurrent2
               setcurrentfocus( 2);
           case editcur1
               updatecur1;
           case cur1plus
               set(editcur1, 'String', num2str(vars.current1 + vars.dcur1));
               updatecur1;
           case cur1minus
               set(editcur1, 'String', num2str(vars.current1 - vars.dcur1));
               updatecur1;
           case editcur2
               updatecur2;
           case cur2plus
               set(editcur2, 'String', num2str(vars.current2 + vars.dcur2));
               updatecur2;
           case cur2minus
               set(editcur2, 'String', num2str(vars.current2 - vars.dcur2));
               updatecur2;
           case edittau2
               updatetau2;
           case tau2plus
               set(edittau2, 'String', num2str(vars.tau2 + vars.dtau2));
               updatetau2;
           case tau2minus
               set(edittau2, 'String', num2str(vars.tau2 - vars.dtau2));
               updatetau2;
           otherwise
               disp( 'error');
       end
       resetgraphs();          
   end

    function updateslider 
        vars.cacheNmax = floor( 1 + ((80-1)* get( myslider, 'Value')) ) ;
    end

    function updateres
        [temp status] = str2num( get(editres, 'String'));
        if status
            if temp < vars.resmin set(editres, 'String', vars.resmin); end;
            if temp > vars.resmax set(editres, 'String', vars.resmax); end;
        else
            set(editres, 'String', vars.resistance);
        end
        vars.resistance = str2num( get(editres, 'String'));
    end

    function updatecap
        [temp status] = str2num( get(editcap, 'String'));
        if status
            if temp < vars.capmin set(editcap, 'String', vars.capmin); end;
            if temp > vars.capmax set(editcap, 'String', vars.capmax); end;
        else
            set(editcap, 'String', vars.capacitance);
        end
        vars.capacitance = str2num( get(editcap, 'String'));
    end

    function updatecur1
        [temp status] = str2num( get(editcur1, 'String'));
        if status
            if temp < vars.cur1min set(editcur1, 'String', vars.cur1min); end;
            if temp > vars.cur1max set(editcur1, 'String', vars.cur1max); end;
        else
            set(editcur1, 'String', vars.current1);
        end
        vars.current1 = str2num( get(editcur1, 'String'));
        setcurrentfocus( 1);
    end

    function updatecur2
        [temp status] = str2num( get(editcur2, 'String'));
        if status
            if temp < vars.cur2min set(editcur2, 'String', vars.cur2min); end;
            if temp > vars.cur2max set(editcur2, 'String', vars.cur2max); end;
        else
            set(editcur2, 'String', vars.current2);
        end
        vars.current2 = str2num( get(editcur2, 'String'));
        setcurrentfocus( 2);
    end

    function updatetau2
        [temp status] = str2num( get(edittau2, 'String'));
        if status
            if temp < vars.tau2min set(edittau2, 'String', vars.tau2min); end;
            if temp > vars.tau2max set(edittau2, 'String', vars.tau2max); end;
        else
            set(edittau2, 'String', vars.tau2);
        end
        vars.tau2 = str2num( get(edittau2, 'String'));
        setcurrentfocus( 2);
    end

    % this just sets our currents depending on if its pulse or continuous
    function setcurrentfocus( ww)
        switch ww
            case 2
                set( bcurrent1, 'Value', 0);
                set( bcurrent2, 'Value', 1);
                vars.currentswitch = 2;
            otherwise
                set( bcurrent1, 'Value', 1);
                set( bcurrent2, 'Value', 0);
                vars.currentswitch = 1;
        end
    end

    function resetgraphs
        set( lineV, 'XData', []);
        set( lineV, 'YData', []);
        set( lineIres, 'XData', []);
        set( lineIres, 'YData', []);
        set( lineIcap, 'XData', []);
        set( lineIcap, 'YData', []);
        switch vars.currentswitch
            case 2
                if 5+vars.tau2>=vars.Tmax
                    set( lineIin, 'XData', [0, 5, 5, vars.Tmax]);
                    set( lineIin, 'YData', [0, 0, vars.current2, vars.current2]);
                else
                    set( lineIin, 'XData', [0, 5, 5, 5+vars.tau2, 5+vars.tau2, vars.Tmax]);
                    set( lineIin, 'YData', [0, 0, vars.current2, vars.current2, 0, 0]);
                end
            otherwise
                set( lineIin, 'XData', [0 vars.Tmax]);
                set( lineIin, 'YData', [vars.current1, vars.current1]);
        end
        vars.T = 0;
        vars.V = vars.Vrest;
        myiterate;
    end

    function myiterate
      if vars.plotting==1
          return;
      else vars.plotting=1;
      end
      while vars.T <= vars.Tmax
          while length(vars.cacheT)<vars.cacheNmax & vars.T<=vars.Tmax
              vars.T = vars.T+vars.dT;
              % set the input current
              switch vars.currentswitch
                  case 2 % current pulse
                      if (vars.T>5 && vars.T<5+vars.tau2) vars.Iinput = vars.current2; else vars.Iinput = 0; end;
                  otherwise % constant current
                      vars.Iinput = vars.current1;
              end
              % iterate
              %vars.V = vars.V + vars.dT * (vars.Iinput*vars.resistance - vars.V)/(vars.resistance*vars.capacitance);
              vars.V = vars.V + vars.dT * (vars.Iinput*vars.resistance - vars.V + vars.Vrest)/(vars.resistance*vars.capacitance);
              if vars.V >vars.Vthreshold vars.V=vars.Vreset; end;
              vars.cacheT = [vars.cacheT, vars.T];
              vars.cacheV = [vars.cacheV, vars.V];
              vars.cacheIres = [vars.cacheIres, (vars.V-vars.Vrest)/vars.resistance];  %vars.V/vars.resistance];
              vars.cacheIcap = [vars.cacheIcap, vars.Iinput-(vars.V-vars.Vrest)/vars.resistance];%vars.Iinput-vars.V/vars.resistance];
          end
          set( lineV, 'XData', [get(lineV, 'XData'), vars.cacheT]);
          set( lineV, 'YData', [get(lineV, 'YData'), vars.cacheV]);
          set( lineIres, 'XData', [get(lineIres, 'XData'), vars.cacheT]);
          set( lineIres, 'YData', [get(lineIres, 'YData'), vars.cacheIres]);
          set( lineIcap, 'XData', [get(lineIcap, 'XData'), vars.cacheT]);
          set( lineIcap, 'YData', [get(lineIcap, 'YData'), vars.cacheIcap]);
          
          adjustaxI;
          drawnow;
          vars.cacheT = [];
          vars.cacheV = [];
          vars.cacheIres = [];
          vars.cacheIcap = [];
        if vars.plotting==0, return; end
      end
      vars.plotting=0;
    end
    
    % 
    function adjustaxI
        yvals1 = get(lineIcap, 'YData');
        yvals2 = get(lineIres, 'YData');
        vars.axImaxN = ceil(min( 20, max([yvals1, yvals2, 1])));
        vars.axIminN = floor(max(-20, min([yvals1, yvals2,-1])));
        gp = 0.1*(vars.axImaxN-vars.axIminN);
        set(axIout, 'YLim', [min(-1, min([yvals1,yvals2])), max(1, max([yvals1,yvals2]))]);
        
        % set(axIout, 'YLim', [vars.axIminN-gp, vars.axImaxN+gp]);
        % set(axIout, 'YTick', [vars.axIminN, vars.axImaxN]);
        % set(axIout, 'YTickLabel', [num2str(vars.axIminN), '|', num2str(vars.axImaxN)]);
    end

    function drawCircuit

        % input arrow
        text('Parent', axCircuit, 'String', 'I', 'Position', [15 115], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Iin);
        text('Parent', axCircuit, 'String', 'input', 'Position', [15 115], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Iin);
        plot(axCircuit, [70-20 150-50], [105+10 105+10], '-k', 'LineWidth', 1, 'Color', mycolors.Iin)
        patch(150-50+[0 10 0], 105+10+[-4 0 4], mycolors.Iin, 'EdgeColor', mycolors.Iin, 'Parent', axCircuit)
        % capcitor arrow
        text('Parent', axCircuit, 'String', 'I', 'Position', [210+30 75], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Icap);
        text('Parent', axCircuit, 'String', 'cap', 'Position', [210+30 75], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Icap);
        plot(axCircuit, [210-10 210+10 210+10], [85+10 85+10 85-10], '-k', 'LineWidth', 1, 'Color', mycolors.Icap)
        patch(210+10+[-4 0 4], 85-10-[0 10 0], mycolors.Icap, 'EdgeColor', mycolors.Icap, 'Parent', axCircuit)
        % resistor arrow
        text('Parent', axCircuit, 'String', 'I', 'Position', [90-33 75], 'units', 'pixels', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', mycolors.Ires);
        text('Parent', axCircuit, 'String', 'res', 'Position', [90-33 75], 'units', 'pixels', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', mycolors.Ires);
        plot(axCircuit, [90+10 90-10 90-10], [85+10 85+10 85-10], '-k', 'LineWidth', 1, 'Color', mycolors.Ires)
        patch(90-10+[-4 0 4], 85-10-[0 10 0], mycolors.Ires, 'EdgeColor', mycolors.Ires, 'Parent', axCircuit)

        % input
        plot(axCircuit, [70 150], [105 105], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [150 150], [85 105], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [90 210], [85 85], '-k', 'LineWidth', 1.5)
        
        % res
        plot(axCircuit, [90 90], [60 85], '-k', 'LineWidth', 1.5)
        plot(axCircuit, 90+5*[0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0], linspace(30,60,15), '-k', 'LineWidth', 1.5)
        plot(axCircuit, [90 90], [22 30], '-k', 'LineWidth', 1.5)

        % cap
        plot(axCircuit, [210 210], [85 45+5], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [210-15 210+15], [45+5 45+5], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [210-15 210+15], [45-5 45-5], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [210 210], [22 45-5], '-k', 'LineWidth', 1.5)

        % ground
        plot(axCircuit, [150-60 150+60], [22 22], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [150 150], [12 22], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [150-1 150+1], [4 4], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [150-4 150+4], [8 8], '-k', 'LineWidth', 1.5)
        plot(axCircuit, [150-7 150+7], [12 12], '-k', 'LineWidth', 1.5)
        
        uicontrol('Parent', pcircuit, 'Style', 'text', 'String', 'R', 'Units', 'pixels', 'Position', [90+10 35 20 20], 'HorizontalAlignment', 'left', 'FontSize', 15);
        uicontrol('Parent', pcircuit, 'Style', 'text', 'String', 'C', 'Units', 'pixels', 'Position', [210-45 35 20 20], 'HorizontalAlignment', 'right', 'FontSize', 15);
        
    end


end