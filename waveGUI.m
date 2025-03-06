function waveGUI()
% waveGUI
% Creates a GUI that allows the user to:
%   - Select wave type (regular / irregular + specific model)
%   - Enter simulation parameters
%   - Save/Load config to .MAT files
%   - Plot wave time series in the GUI (preview)
%   - (Optional) Generate ASCII output file for SIMA

%% === Create the figure & main layout ===
% A simple figure with two main areas: controls (left) and axes (right).
figWidth  = 950;
figHeight = 500;
fig = figure('Name','Wave Generator GUI',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'Position',[100 100 figWidth figHeight]);

% Panel for user controls
% We'll place it on the left side
ctrlPanelWidth = 300;
ctrlPanel = uipanel('Parent',fig,...
    'Units','pixels',...
    'Title','User Inputs',...
    'FontWeight','bold',...
    'Position',[10 10 ctrlPanelWidth figHeight-20]);

% Axes for the preview plot (on the right)
ax = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[ctrlPanelWidth+40, 60, figWidth - ctrlPanelWidth - 70, figHeight - 120]);
title(ax, 'Wave Elevation Preview');
xlabel(ax, 'Time [s]');
ylabel(ax, '\eta(t) [m]');
grid(ax,'on');

%% === Common UI elements (shared by both wave types) ===
% We will store references to all UI components in a struct, and store it in the figure's UserData

% 1) waveType selection (Regular / Irregular)
uicontrol('Parent',ctrlPanel,'Style','text','String','Wave Type:',...
    'HorizontalAlignment','left','Position',[10 figHeight-60-10 80 20]);
waveTypePopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',{'Regular','Irregular'},...
    'Value',1,...
    'Position',[100 figHeight-60-10 150 25],...
    'Callback',@onWaveTypeChanged);

% 2) For irregular wave, we need an extra selection for the spectrum model
modelList = {'JONSWAP 2-param', 'JONSWAP 3-param','JONSWAP 6-param',...
    'Torsethaugen','Ochi-Hubble','PM (1-param)',...
    'PM (2-param, Tavg)','PM (2-param, Tzc)'};
uicontrol('Parent',ctrlPanel,'Style','text','String','Spectrum Model:',...
    'HorizontalAlignment','left','Position',[10 figHeight-90-10 100 20]);
spectrumPopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',modelList,...
    'Value',1,...
    'Position',[120 figHeight-90-10 130 25]);

% 3) Common parameters: time step, N, ramp time, random seed
commonYStart = figHeight - 130;
uicontrol('Parent',ctrlPanel,'Style','text','String','Time Step:',...
    'HorizontalAlignment','left','Position',[10 commonYStart 100 20]);
edtTimeStep = uicontrol('Parent',ctrlPanel,'Style','edit','String','0.5',...
    'Position',[120 commonYStart 50 25]);

uicontrol('Parent',ctrlPanel,'Style','text','String','Num Samples:',...
    'HorizontalAlignment','left','Position',[10 commonYStart-30 100 20]);
edtNumSamples = uicontrol('Parent',ctrlPanel,'Style','edit','String','2048',...
    'Position',[120 commonYStart-30 50 25]);

uicontrol('Parent',ctrlPanel,'Style','text','String','Ramp Duration:',...
    'HorizontalAlignment','left','Position',[10 commonYStart-60 100 20]);
edtRamp = uicontrol('Parent',ctrlPanel,'Style','edit','String','20',...
    'Position',[120 commonYStart-60 50 25]);

uicontrol('Parent',ctrlPanel,'Style','text','String','Random Seed:',...
    'HorizontalAlignment','left','Position',[10 commonYStart-90 100 20]);
edtRandomSeed = uicontrol('Parent',ctrlPanel,'Style','edit','String','42',...
    'Position',[120 commonYStart-90 50 25]);

%% === REGULAR wave parameter UI ===
% We'll place them further down in the panel. We'll show/hide them as needed.
regularYStart = commonYStart-130;
txtH = uicontrol('Parent',ctrlPanel,'Style','text','String','Wave Height, H:',...
    'HorizontalAlignment','left','Position',[10 regularYStart 120 20]);
edtH = uicontrol('Parent',ctrlPanel,'Style','edit','String','2.0',...
    'Position',[130 regularYStart 50 25]);

txtT = uicontrol('Parent',ctrlPanel,'Style','text','String','Wave Period, T:',...
    'HorizontalAlignment','left','Position',[10 regularYStart-30 120 20]);
edtT = uicontrol('Parent',ctrlPanel,'Style','edit','String','8.0',...
    'Position',[130 regularYStart-30 50 25]);

txtBeta = uicontrol('Parent',ctrlPanel,'Style','text','String','Phase, beta:',...
    'HorizontalAlignment','left','Position',[10 regularYStart-60 120 20]);
edtBeta = uicontrol('Parent',ctrlPanel,'Style','edit','String','0',...
    'Position',[130 regularYStart-60 50 25]);

%% === IRREGULAR wave parameter UI (some are model-specific) ===
% We'll create a group box or just group them together. Show/hide as needed.
irrYStart = regularYStart - 100;

% 2-Param/3-Param etc. mostly revolve around Hs, Tpeak, gamma, etc.
txtSIWAHE = uicontrol('Parent',ctrlPanel,'Style','text','String','Significant Hs:',...
    'HorizontalAlignment','left','Position',[10 irrYStart 110 20]);
edtSIWAHE = uicontrol('Parent',ctrlPanel,'Style','edit','String','3.0',...
    'Position',[130 irrYStart 50 25]);

txtTpeak = uicontrol('Parent',ctrlPanel,'Style','text','String','Peak Period, Tp:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-30 110 20]);
edtTpeak = uicontrol('Parent',ctrlPanel,'Style','edit','String','10.0',...
    'Position',[130 irrYStart-30 50 25]);

txtGamma = uicontrol('Parent',ctrlPanel,'Style','text','String','Gamma:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-60 110 20]);
edtGamma = uicontrol('Parent',ctrlPanel,'Style','edit','String','3.3',...
    'Position',[130 irrYStart-60 50 25]);

% 6-param JONSWAP:
txtOmegaP = uicontrol('Parent',ctrlPanel,'Style','text','String','Omega_p:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-90 110 20]);
edtOmegaP = uicontrol('Parent',ctrlPanel,'Style','edit','String','0.6283',...
    'Position',[130 irrYStart-90 50 25]); % ~ 2*pi/10

txtAlpha = uicontrol('Parent',ctrlPanel,'Style','text','String','Alpha:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-120 110 20]);
edtAlpha = uicontrol('Parent',ctrlPanel,'Style','edit','String','0.0081',...
    'Position',[130 irrYStart-120 50 25]);

txtBetaVal = uicontrol('Parent',ctrlPanel,'Style','text','String','Beta:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-150 110 20]);
edtBetaVal = uicontrol('Parent',ctrlPanel,'Style','edit','String','1.25',...
    'Position',[130 irrYStart-150 50 25]);

txtSigA = uicontrol('Parent',ctrlPanel,'Style','text','String','Sigma_a:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-180 110 20]);
edtSigA = uicontrol('Parent',ctrlPanel,'Style','edit','String','0.07',...
    'Position',[130 irrYStart-180 50 25]);

txtSigB = uicontrol('Parent',ctrlPanel,'Style','text','String','Sigma_b:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-210 110 20]);
edtSigB = uicontrol('Parent',ctrlPanel,'Style','edit','String','0.09',...
    'Position',[130 irrYStart-210 50 25]);

% Additional fields: e.g. Ochi-Hubble, Torsethaugen might just rely on Hs and Tpeak
% PM 2 param with Tavg or Tzc
txtTavg = uicontrol('Parent',ctrlPanel,'Style','text','String','Tavg:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-240 110 20]);
edtTavg = uicontrol('Parent',ctrlPanel,'Style','edit','String','6.0',...
    'Position',[130 irrYStart-240 50 25]);

txtTzc = uicontrol('Parent',ctrlPanel,'Style','text','String','Tzc:',...
    'HorizontalAlignment','left','Position',[10 irrYStart-270 110 20]);
edtTzc = uicontrol('Parent',ctrlPanel,'Style','edit','String','7.0',...
    'Position',[130 irrYStart-270 50 25]);

%% === Buttons (Plot, Save Config, Load Config, Generate ASCII) ===
btnYPos = 30;

% Plot button
btnPlot = uicontrol('Parent',ctrlPanel,'Style','pushbutton','String','Plot Wave',...
    'FontWeight','bold',...
    'Position',[10 btnYPos 80 40],...
    'Callback',@onPlotWave);

% Generate ASCII button
btnGen = uicontrol('Parent',ctrlPanel,'Style','pushbutton','String','Generate ASCII',...
    'Position',[100 btnYPos 100 40],...
    'Callback',@onGenerateASCII);

% Save config
btnSave = uicontrol('Parent',ctrlPanel,'Style','pushbutton','String','Save Config',...
    'Position',[210 btnYPos 80 20],...
    'Callback',@onSaveConfig);

% Load config
btnLoad = uicontrol('Parent',ctrlPanel,'Style','pushbutton','String','Load Config',...
    'Position',[210 btnYPos+20 80 20],...
    'Callback',@onLoadConfig);

%% === Store handles in the figure's UserData, then do initial update
handles.ax = ax;
handles.waveTypePopup    = waveTypePopup;
handles.spectrumPopup    = spectrumPopup;
handles.edtTimeStep     = edtTimeStep;
handles.edtNumSamples   = edtNumSamples;
handles.edtRamp         = edtRamp;
handles.edtRandomSeed   = edtRandomSeed;
handles.txtH            = txtH;
handles.edtH            = edtH;
handles.txtT            = txtT;
handles.edtT            = edtT;
handles.txtBeta         = txtBeta;
handles.edtBeta         = edtBeta;
handles.txtSIWAHE       = txtSIWAHE;
handles.edtSIWAHE       = edtSIWAHE;
handles.txtTpeak        = txtTpeak;
handles.edtTpeak        = edtTpeak;
handles.txtGamma        = txtGamma;
handles.edtGamma        = edtGamma;
handles.txtOmegaP       = txtOmegaP;
handles.edtOmegaP       = edtOmegaP;
handles.txtAlpha        = txtAlpha;
handles.edtAlpha        = edtAlpha;
handles.txtBetaVal      = txtBetaVal;
handles.edtBetaVal      = edtBetaVal;
handles.txtSigA         = txtSigA;
handles.edtSigA         = edtSigA;
handles.txtSigB         = txtSigB;
handles.edtSigB         = edtSigB;
handles.txtTavg         = txtTavg;
handles.edtTavg         = edtTavg;
handles.txtTzc          = txtTzc;
handles.edtTzc          = edtTzc;
handles.btnPlot         = btnPlot;
handles.btnGen          = btnGen;
handles.btnSave         = btnSave;
handles.btnLoad         = btnLoad;

set(fig, 'UserData', handles);

% Hide all irregular param fields initially if 'Regular' is default
updateUIVisibility();

% --- End of waveGUI (main function) ---

%% =============================================================
%  Nested Callback Functions & Utility
%% =============================================================

    function onWaveTypeChanged(~,~)
        % Callback when user changes waveType popup
        updateUIVisibility();
    end

    function updateUIVisibility()
        % Show/hide controls depending on "Regular" vs "Irregular" and the chosen model
        handles = get(fig,'UserData');

        val = get(handles.waveTypePopup,'Value');
        waveTypeList = get(handles.waveTypePopup,'String');
        waveTypeStr = waveTypeList{val};

        if strcmpi(waveTypeStr,'Regular')
            % Show regular wave controls
            set([handles.txtH, handles.edtH, ...
                handles.txtT, handles.edtT, ...
                handles.txtBeta, handles.edtBeta], 'Visible','on');
            % Hide spectrum popup and irregular fields
            set(handles.spectrumPopup, 'Enable','off');
            set([handles.txtSIWAHE, handles.edtSIWAHE, ...
                handles.txtTpeak,  handles.edtTpeak, ...
                handles.txtGamma,  handles.edtGamma, ...
                handles.txtOmegaP, handles.edtOmegaP, ...
                handles.txtAlpha,  handles.edtAlpha, ...
                handles.txtBetaVal,handles.edtBetaVal, ...
                handles.txtSigA,   handles.edtSigA, ...
                handles.txtSigB,   handles.edtSigB, ...
                handles.txtTavg,   handles.edtTavg, ...
                handles.txtTzc,    handles.edtTzc], 'Visible','off');
        else
            % Irregular
            set([handles.txtH, handles.edtH, ...
                handles.txtT, handles.edtT, ...
                handles.txtBeta, handles.edtBeta], 'Visible','off');
            set(handles.spectrumPopup, 'Enable','on');
            % Now show/hide sub-fields based on the chosen model
            modelVal = get(handles.spectrumPopup,'Value');
            modelStr = modelList{modelVal};

            % Start by turning them all off:
            allIrrCtrls = [handles.txtSIWAHE, handles.edtSIWAHE,...
                handles.txtTpeak,  handles.edtTpeak,...
                handles.txtGamma,  handles.edtGamma,...
                handles.txtOmegaP, handles.edtOmegaP,...
                handles.txtAlpha,  handles.edtAlpha,...
                handles.txtBetaVal,handles.edtBetaVal,...
                handles.txtSigA,   handles.edtSigA,...
                handles.txtSigB,   handles.edtSigB,...
                handles.txtTavg,   handles.edtTavg,...
                handles.txtTzc,    handles.edtTzc];
            set(allIrrCtrls, 'Visible','off');

            % Turn on relevant fields:
            switch modelStr
                case 'JONSWAP 2-param'
                    set([handles.txtSIWAHE, handles.edtSIWAHE,...
                        handles.txtTpeak,  handles.edtTpeak], 'Visible','on');
                case 'JONSWAP 3-param'
                    set([handles.txtSIWAHE, handles.edtSIWAHE,...
                        handles.txtTpeak,  handles.edtTpeak,...
                        handles.txtGamma,  handles.edtGamma], 'Visible','on');
                case 'JONSWAP 6-param'
                    set([handles.txtOmegaP, handles.edtOmegaP,...
                        handles.txtAlpha,  handles.edtAlpha,...
                        handles.txtBetaVal,handles.edtBetaVal,...
                        handles.txtGamma,  handles.edtGamma,...
                        handles.txtSigA,   handles.edtSigA,...
                        handles.txtSigB,   handles.edtSigB], 'Visible','on');
                case 'Torsethaugen'
                    % Typically needs Hs, Tp
                    set([handles.txtSIWAHE, handles.edtSIWAHE,...
                        handles.txtTpeak,  handles.edtTpeak], 'Visible','on');
                case 'Ochi-Hubble'
                    % Typically needs Hs, possibly other sub-params
                    set([handles.txtSIWAHE, handles.edtSIWAHE], 'Visible','on');
                case 'PM (1-param)'
                    % Typically just Hs
                    set([handles.txtSIWAHE, handles.edtSIWAHE], 'Visible','on');
                case 'PM (2-param, Tavg)'
                    set([handles.txtSIWAHE, handles.edtSIWAHE,...
                        handles.txtTavg,   handles.edtTavg], 'Visible','on');
                case 'PM (2-param, Tzc)'
                    set([handles.txtSIWAHE, handles.edtSIWAHE,...
                        handles.txtTzc,    handles.edtTzc], 'Visible','on');
            end
        end
    end

    function onPlotWave(~,~)
        % Callback: read params, generate wave, plot on the axes
        [t, eta] = generateWaveFromGUI();
        % Plot
        axes(handles.ax); %#ok<LAXES>
        cla(handles.ax);
        plot(handles.ax, t, eta,'b-','LineWidth',1.2);
        grid(handles.ax,'on');
        xlabel(handles.ax,'Time [s]');
        ylabel(handles.ax,'\eta(t) [m]');
        title(handles.ax,'Preview of Wave Elevation');
    end

    function onGenerateASCII(~,~)
        % Generate wave, then prompt user to save ASCII in SIMA format
        [t, eta] = generateWaveFromGUI();
        prompt = {'Enter output ASCII file name:'};
        dlgTitle = 'Generate SIMA ASCII';
        def = {'waveElevation_sima.dat'};
        answer = inputdlg(prompt,dlgTitle,[1 60],def);
        if isempty(answer), return; end
        outFileName = answer{1};

        dt = t(2) - t(1);
        writeSIMAAscii(outFileName, eta, dt);
        msgbox(['Wrote wave time series to ', outFileName],'Success','modal');
    end

    function onSaveConfig(~,~)
        % Save all current GUI parameters to a .MAT file
        [fileName, pathName] = uiputfile('*.mat','Save Configuration');
        if fileName == 0, return; end
        cfg = getAllParams();
        save(fullfile(pathName, fileName), 'cfg');
    end

    function onLoadConfig(~,~)
        % Load .MAT file and populate the GUI
        [fileName, pathName] = uigetfile('*.mat','Load Configuration');
        if fileName == 0, return; end
        S = load(fullfile(pathName, fileName),'cfg');
        if isfield(S,'cfg')
            setParams(S.cfg);
            updateUIVisibility();
        end
    end

    function cfg = getAllParams()
        % Read all GUI fields into a struct
        cfg.waveTypeVal = get(handles.waveTypePopup,'Value');
        cfg.spectrumVal = get(handles.spectrumPopup,'Value');
        cfg.tStep       = str2double(get(handles.edtTimeStep,'String'));
        cfg.nSamples    = str2double(get(handles.edtNumSamples,'String'));
        cfg.tRamp       = str2double(get(handles.edtRamp,'String'));
        cfg.rSeed       = str2double(get(handles.edtRandomSeed,'String'));

        cfg.H       = str2double(get(handles.edtH,'String'));
        cfg.T       = str2double(get(handles.edtT,'String'));
        cfg.beta    = str2double(get(handles.edtBeta,'String'));

        cfg.Hs      = str2double(get(handles.edtSIWAHE,'String'));
        cfg.Tpeak   = str2double(get(handles.edtTpeak,'String'));
        cfg.gamma   = str2double(get(handles.edtGamma,'String'));
        cfg.omegaP  = str2double(get(handles.edtOmegaP,'String'));
        cfg.alpha   = str2double(get(handles.edtAlpha,'String'));
        cfg.betaVal = str2double(get(handles.edtBetaVal,'String'));
        cfg.sigA    = str2double(get(handles.edtSigA,'String'));
        cfg.sigB    = str2double(get(handles.edtSigB,'String'));
        cfg.Tavg    = str2double(get(handles.edtTavg,'String'));
        cfg.Tzc     = str2double(get(handles.edtTzc,'String'));
    end

    function setParams(cfg)
        % Populate GUI fields from cfg struct
        set(handles.waveTypePopup,'Value', cfg.waveTypeVal);
        set(handles.spectrumPopup,'Value', cfg.spectrumVal);
        set(handles.edtTimeStep,'String', num2str(cfg.tStep));
        set(handles.edtNumSamples,'String', num2str(cfg.nSamples));
        set(handles.edtRamp,'String', num2str(cfg.tRamp));
        set(handles.edtRandomSeed,'String', num2str(cfg.rSeed));

        set(handles.edtH,'String', num2str(cfg.H));
        set(handles.edtT,'String', num2str(cfg.T));
        set(handles.edtBeta,'String', num2str(cfg.beta));

        set(handles.edtSIWAHE,'String', num2str(cfg.Hs));
        set(handles.edtTpeak,'String', num2str(cfg.Tpeak));
        set(handles.edtGamma,'String', num2str(cfg.gamma));
        set(handles.edtOmegaP,'String', num2str(cfg.omegaP));
        set(handles.edtAlpha,'String', num2str(cfg.alpha));
        set(handles.edtBetaVal,'String', num2str(cfg.betaVal));
        set(handles.edtSigA,'String', num2str(cfg.sigA));
        set(handles.edtSigB,'String', num2str(cfg.sigB));
        set(handles.edtTavg,'String', num2str(cfg.Tavg));
        set(handles.edtTzc,'String', num2str(cfg.Tzc));
    end

    function [t, eta] = generateWaveFromGUI()
        % Gather all GUI params, produce the time vector and wave series
        cfg = getAllParams();
        dt = cfg.tStep;
        N = cfg.nSamples;
        t = (0:N-1).'*dt;

        switch cfg.waveTypeVal
            case 1  % Regular
                % from the popup 'String', 'Value'=1 -> 'Regular'
                omega = 2*pi/cfg.T;
                rampFactor = rampUp(t, cfg.tRamp);
                eta = 0.5*cfg.H * cos(omega*t + cfg.beta) .* rampFactor;
            otherwise
                % Irregular
                rng(cfg.rSeed);
                % Which model?
                modelStr = modelList{cfg.spectrumVal};
                specFun = @(w) zeros(size(w)); % default
                switch modelStr
                    case 'JONSWAP 2-param'
                        specFun = @(w) jonswap2(w, cfg.Hs, cfg.Tpeak);
                    case 'JONSWAP 3-param'
                        specFun = @(w) jonswap3(w, cfg.Hs, cfg.Tpeak, cfg.gamma);
                    case 'JONSWAP 6-param'
                        specFun = @(w) jonswap6(w, cfg.omegaP, cfg.alpha, cfg.betaVal, cfg.gamma, cfg.sigA, cfg.sigB);
                    case 'Torsethaugen'
                        specFun = @(w) torsethaugenSpectrum(w, cfg.Hs, cfg.Tpeak);
                    case 'Ochi-Hubble'
                        specFun = @(w) ochiHubbleSpectrum(w, cfg.Hs);
                    case 'PM (1-param)'
                        specFun = @(w) pmSpectrum1(w, cfg.Hs);
                    case 'PM (2-param, Tavg)'
                        specFun = @(w) pmSpectrum2_AVWAPE(w, cfg.Hs, cfg.Tavg);
                    case 'PM (2-param, Tzc)'
                        specFun = @(w) pmSpectrum2_ZCWAPE(w, cfg.Hs, cfg.Tzc);
                end
                eta = generateIrregularWave(t, cfg.tRamp, specFun);
        end
    end
end

%% ========================================================================
%  Wave Utility Functions
%% ========================================================================

function r = rampUp(t, t_ramp)
% Simple linear ramp from 0 to 1 over [0, t_ramp].
r = ones(size(t));
idx = (t < t_ramp);
r(idx) = t(idx)/t_ramp;
end

function eta = generateIrregularWave(t, t_ramp, specFun)
N = length(t);
dt = t(2)-t(1);

df = 1/(N*dt);
f  = (0:(N-1)).' * df;
w  = 2*pi*f;

% Evaluate spectral density
S = specFun(w);
S(w<=0) = 0;

dw = 2*pi*df;
A  = sqrt(2*S*dw);
phi = 2*pi*rand(size(w));

Npos = floor(N/2);
eta_temp = zeros(size(t));
for i = 2:Npos
    eta_temp = eta_temp + A(i)*cos(w(i)*t + phi(i));
end

rampFactor = rampUp(t, t_ramp);
eta = rampFactor.*eta_temp;
end

%% --------------------- Example Spectral Models -------------------------%%
function S = jonswap2(w, Hs, Tp)
g = 9.81;
wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
gamma = 3.3;
sigma = 0.07.*(w<=wp) + 0.09.*(w>wp);

S = alpha_PM ./ (w.^5).*exp(-5/4*(wp./w).^4);
exponent = exp(-((w - wp).^2)./(2*sigma.^2*wp.^2));
S = S.*(gamma.^exponent);

% Scale to match Hs
S = scaleSpectrumToHs(S, w, Hs);
end

function S = jonswap3(w, Hs, Tp, gammaVal)
g = 9.81;
wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
sigma = 0.07.*(w<=wp) + 0.09.*(w>wp);

S = alpha_PM./(w.^5).*exp(-5/4*(wp./w).^4);
exponent = exp(-((w-wp).^2)./(2*sigma.^2.*wp.^2));
S = S.*(gammaVal.^exponent);

S = scaleSpectrumToHs(S, w, Hs);
end

function S = jonswap6(w, w_p, alpha, betaVal, gammaVal, sigA, sigB)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_  = w(idx);
sigma = sigA.*(w_<=w_p) + sigB.*(w_>w_p);

S_PM = alpha*g^2./(w_.^5).*exp(-betaVal*(w_p./w_).^4);
peakExp = exp(-((w_-w_p).^2)./(2*sigma.^2*w_p.^2));
G = gammaVal.^peakExp;
S_(idx) = S_PM.*G;

S = S_;
% We do NOT automatically scale to Hs here, because the user
% might specify alpha, gamma, etc. to get a known Hs.
% If you do want to scale to an Hs, call scaleSpectrumToHs.
end

function S = torsethaugenSpectrum(w, Hs, Tp)
% Placeholder: Return a JONSWAP-like shape
S = jonswap2(w, Hs, Tp);
% Real Torsethaugen is double-peaked. Replace with your actual formula.
end

function S = ochiHubbleSpectrum(w, Hs)
% Placeholder. Real Ochi-Hubble is typically a sum of two gamma distributions.
% Here, we do a JONSWAP(2) with a guessed peak.
TpeakFake = 8.0;
S = jonswap2(w, Hs, TpeakFake);
end

function S = pmSpectrum1(w, Hs)
g = 9.81;
alpha = 0.0081;
% Just guess a peak freq:
w_p = 0.44*g/((Hs/4)^1.0);
S = alpha*g^2./(w.^5).*exp(-1.25*(w_p./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S, w, Hs);
end

function S = pmSpectrum2_AVWAPE(w, Hs, Tavg)
g = 9.81;
w_ave = 2*pi/Tavg;
alpha = 0.0081;
S = alpha*g^2./(w.^5).*exp(-1.25*(w_ave./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S, w, Hs);
end

function S = pmSpectrum2_ZCWAPE(w, Hs, Tzc)
g = 9.81;
w_zc = 2*pi/Tzc;
alpha = 0.0081;
S = alpha*g^2./(w.^5).*exp(-1.25*(w_zc./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S, w, Hs);
end

function Sscaled = scaleSpectrumToHs(S, w, Hs)
% Numeric integration of S(omega) with respect to omega => m0
% Then Hs = 4 * sqrt(m0).
dw = mean(diff(w));
m0 = trapz(w, S);
Hs_model = 4*sqrt(m0);
scaleFactor = (Hs / Hs_model)^2;
Sscaled = S * scaleFactor;
end

%% ========================================================================
%  Subfunction: Write ASCII file for SIMA
%% ========================================================================
function writeSIMAAscii(outFileName, eta, t_step)
fid = fopen(outFileName, 'w');
if fid < 0
    error('Could not open file: %s', outFileName);
end
N = length(eta);

fprintf(fid, '%d\n', N);
fprintf(fid, '%.4f\n', t_step);
fprintf(fid, 'Generated by waveGUI\n');
fprintf(fid, 'Wave elevation time series\n');

valsPerLine = 5;
for i = 1:N
    fprintf(fid, '%12.6f', eta(i));
    if mod(i, valsPerLine) == 0
        fprintf(fid, '\n');
    else
        fprintf(fid, ' ');
    end
end
if mod(N, valsPerLine) ~= 0
    fprintf(fid, '\n');
end
fclose(fid);
end
