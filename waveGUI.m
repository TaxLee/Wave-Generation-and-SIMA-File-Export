function waveGUI()
% waveGUI - Wave Generator GUI
%
% Description:
%   A MATLAB GUI to generate regular and irregular wave elevation signals,
%   plot the wave, save wave data, and generate SIMA-compatible ASCII files.
%   The tool features a ramp function to gradually introduce wave amplitudes,
%   ensuring a smooth start to the simulated waves.
%
% Features:
%   - Displays all parameters at once.
%   - Disables irrelevant fields based on the selected wave/spectrum model.
%   - Improved alignment and intuitive UI layout.
%   - Default window size of 1024 x 666 pixels (adjust as preferred).
%   - Load/Save configuration options at the top.
%   - Supports generating Regular waves and various Irregular waves, including:
%       * JONSWAP (2-param, 3-param, 6-param)
%       * Torsethaugen (testing)
%       * Ochi-Hubble (testing)
%       * PM variants (1-param, 2-param with Tavg or Tzc)
%
% Changes in this version (v1.1):
%   - Adds a "please wait" popup during wave generation, and a success popup
%     that auto‐closes after 5 seconds.
%   - Fixes JONSWAP 2-param, JONSWAP 3-param, Torsethaugen (testing), and Ochi-Hubble (testing)
%     so they skip w = 0 in the spectrum calculation (avoiding division by zero).
%
% **Author:** Shuijin Li
% **Email:** lishuijin@nbu.edu.cn
% **Date:** 2025.03.06

%% --------------------------------------------------------------------
%                     Figure & Panel Setup
%% --------------------------------------------------------------------
figWidth  = 1024;
figHeight = 666;
fig = figure('Name','Wave Generator GUI',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','figure',... % Includes Zoom/Pan/Home
    'Resize','off',...
    'Position',[100 100 figWidth figHeight]);

panelWidth  = 380;
panelHeight = figHeight - 20;
ctrlPanel = uipanel('Parent',fig,...
    'Units','pixels',...
    'Title','User Inputs',...
    'FontWeight','bold',...
    'Position',[10 10 panelWidth panelHeight]);

% Axes for plotting on the right side
ax = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[panelWidth+60, 80, figWidth - (panelWidth+80), figHeight - 120]);
title(ax, 'Wave Elevation Preview');
xlabel(ax, 'Time [s]');
ylabel(ax, '\eta(t) [m]');
grid(ax,'on');

%% --------------------------------------------------------------------
%                  Common Layout Parameters
%% --------------------------------------------------------------------
yTop = panelHeight - 40;
rowHeight = 20;
labelWidth = 200;
editWidth  = 60;
xLabel = 10;
xEdit  = xLabel + labelWidth + 10;

%% --------------------------------------------------------------------
%              1) Load/Save Config Buttons (Top of Panel)
%% --------------------------------------------------------------------
btnLoadCfg = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Load Config','Position',[xLabel yTop 90 rowHeight],...
    'Callback',@onLoadConfig);

btnSaveCfg = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Save Config','Position',[xLabel+100 yTop 90 rowHeight],...
    'Callback',@onSaveConfig);

yTop = yTop - rowHeight - 10;

%% --------------------------------------------------------------------
%              2) Wave Type & Spectrum Model
%% --------------------------------------------------------------------
uicontrol('Parent',ctrlPanel,'Style','text','String','Wave Type:',...
    'HorizontalAlignment','left',...
    'Position',[xLabel yTop labelWidth rowHeight]);
waveTypePopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',{'Regular','Irregular'},...
    'Value',1,...
    'Position',[xEdit yTop editWidth+40 rowHeight],...
    'Callback',@onWaveTypeChanged);

yTop = yTop - rowHeight - 5;

% modelList = {'JONSWAP 2-param','JONSWAP 3-param','JONSWAP 6-param',...
%     'Torsethaugen','Ochi-Hubble','PM (1-param)',...
%     'PM (2-param, Tavg)','PM (2-param, Tzc)'};
modelList = {'JONSWAP 2-param','JONSWAP 3-param','JONSWAP 6-param',...
    'PM (1-param)','PM (2-param, Tavg)','PM (2-param, Tzc)'};
uicontrol('Parent',ctrlPanel,'Style','text','String','Spectrum Model:',...
    'HorizontalAlignment','left',...
    'Position',[xLabel yTop labelWidth rowHeight]);
spectrumPopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',modelList,...
    'Value',1,...
    'Position',[xEdit yTop editWidth+70 rowHeight],...
    'Callback',@onSpectrumModelChanged);

yTop = yTop - rowHeight - 15;

%% --------------------------------------------------------------------
%         3) Common Simulation Parameters
%% --------------------------------------------------------------------
labels = {'Time Step (s):', 'Number of Samples:', 'Ramp Duration (sec):', 'Random Seed (for irregular waves):'};
defaultVals = {'0.1','40000','0','1'};
handlesCommon = [];

for ii = 1:numel(labels)
    uicontrol('Parent',ctrlPanel,'Style','text','String',labels{ii},...
        'HorizontalAlignment','left',...
        'Position',[xLabel yTop labelWidth rowHeight]);
    handlesCommon(ii) = uicontrol('Parent',ctrlPanel,'Style','edit',...
        'String',defaultVals{ii},...
        'Position',[xEdit yTop editWidth rowHeight]);
    yTop = yTop - rowHeight - 5;
end

edtTimeStep   = handlesCommon(1);
edtNumSamples = handlesCommon(2);
edtRamp       = handlesCommon(3);
edtRandomSeed = handlesCommon(4);

yTop = yTop - 10;

%% --------------------------------------------------------------------
%         4) Regular Wave Parameters (H, T, beta)
%% --------------------------------------------------------------------
regLabels = {'Wave Height, H:', 'Wave Period, T:', 'Phase, beta:'};
regDefault = {'9.96','7.99','0'};
handlesRegular = [];

for ii = 1:numel(regLabels)
    uicontrol('Parent',ctrlPanel,'Style','text','String',regLabels{ii},...
        'HorizontalAlignment','left',...
        'Position',[xLabel yTop labelWidth rowHeight]);
    handlesRegular(ii) = uicontrol('Parent',ctrlPanel,'Style','edit',...
        'String',regDefault{ii},...
        'Position',[xEdit yTop editWidth rowHeight]);
    yTop = yTop - rowHeight - 5;
end

edtH    = handlesRegular(1);
edtT    = handlesRegular(2);
edtBeta = handlesRegular(3);

yTop = yTop - 15;

%% --------------------------------------------------------------------
%         5) Irregular Wave Parameters (Plain text labels)
%% --------------------------------------------------------------------
irrLabels = { ...
    'Significant wave height, H_s:', ...
    'Peak period, T_p:', ...
    'Peakedness parameter, γ:', ...
    'Peak frequency, ω_p:', ...
    'Spectrum parameter, α:', ...
    'Form parameter, β:', ...
    'Spectrum parameter, σ_a:', ...
    'Spectrum parameter, σ_b:', ...
    'Average wave period, T_{avg}:', ...
    'Zero-crossing wave period, T_z:'};

irrDefaults = {'3.0','10.0','3.3','0.6283','0.0081','1.25','0.07','0.09','6.0','7.0'};
handlesIrr = [];

for ii = 1:numel(irrLabels)
    uicontrol('Parent', ctrlPanel, ...
        'Style', 'text', ...
        'String', irrLabels{ii}, ...
        'HorizontalAlignment', 'left', ...
        'Position', [xLabel, yTop, labelWidth, rowHeight]);

    handlesIrr(ii) = uicontrol('Parent', ctrlPanel, ...
        'Style', 'edit', ...
        'String', irrDefaults{ii}, ...
        'Position', [xEdit, yTop, editWidth, rowHeight]);

    yTop = yTop - (rowHeight + 5);
end

edtSIWAHE  = handlesIrr(1);
edtTpeak   = handlesIrr(2);
edtGamma   = handlesIrr(3);
edtOmegaP  = handlesIrr(4);
edtAlpha   = handlesIrr(5);
edtBetaVal = handlesIrr(6);
edtSigA    = handlesIrr(7);
edtSigB    = handlesIrr(8);
edtTavg    = handlesIrr(9);
edtTzc     = handlesIrr(10);

%% --------------------------------------------------------------------
%    6) Bottom Buttons: Plot, Generate ASCII, Save Wave Data
%% --------------------------------------------------------------------
btnWidth  = 100;
btnHeight = 30;
spacing   = 10;
yBottom = 40;  % vertical position for bottom row buttons

% Plot Wave
btnPlot = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'Position',[xLabel yBottom btnWidth+20 btnHeight],...
    'String','Plot Wave','FontWeight','bold',...
    'Callback',@onPlotWave);

% Save Wave Data
btnSaveData = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Save Wave Data',...
    'Position',[xLabel+btnWidth+spacing+20 yBottom btnWidth btnHeight],...
    'Callback',@onSaveWaveData);

% Generate ASCII
btnGen = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Generate ASCII',...
    'Position',[xLabel+(btnWidth+spacing+20)*1.85 yBottom btnWidth+20 btnHeight],...
    'Callback',@onGenerateASCII);

%% --------------------------------------------------------------------
%             Store Handles in the Figure's UserData
%% --------------------------------------------------------------------
handles.ax = ax;

handles.waveTypePopup  = waveTypePopup;
handles.spectrumPopup  = spectrumPopup;
handles.modelList      = modelList;

handles.edtTimeStep   = edtTimeStep;
handles.edtNumSamples = edtNumSamples;
handles.edtRamp       = edtRamp;
handles.edtRandomSeed = edtRandomSeed;

handles.edtH      = edtH;
handles.edtT      = edtT;
handles.edtBeta   = edtBeta;

handles.edtSIWAHE  = edtSIWAHE;
handles.edtTpeak   = edtTpeak;
handles.edtGamma   = edtGamma;
handles.edtOmegaP  = edtOmegaP;
handles.edtAlpha   = edtAlpha;
handles.edtBetaVal = edtBetaVal;
handles.edtSigA    = edtSigA;
handles.edtSigB    = edtSigB;
handles.edtTavg    = edtTavg;
handles.edtTzc     = edtTzc;

handles.btnLoadCfg = btnLoadCfg;
handles.btnSaveCfg = btnSaveCfg;
handles.btnSaveData= btnSaveData;
handles.btnPlot    = btnPlot;
handles.btnGen     = btnGen;

set(fig,'UserData',handles);

% Initial UI refresh
updateUIEnable();

%% ====================================================================
%  Callbacks
%% ====================================================================

    function onLoadConfig(~,~)
        [fileName, pathName] = uigetfile('*.mat','Load Configuration');
        if fileName == 0, return; end
        S = load(fullfile(pathName,fileName),'cfg');
        if isfield(S,'cfg')
            setParams(S.cfg);
            updateUIEnable();
        end
    end

    function onSaveConfig(~,~)
        [fileName, pathName] = uiputfile('*.mat','Save Configuration');
        if fileName == 0, return; end
        cfg = getAllParams();
        save(fullfile(pathName,fileName),'cfg');
    end

    function onWaveTypeChanged(~,~)
        updateUIEnable();
    end

    function onSpectrumModelChanged(~,~)
        updateUIEnable();
    end

    function onPlotWave(~,~)
        % Show "Please wait" box
        hWait = msgbox('Generating wave. Please wait...','Please wait','help');
        drawnow;  % Ensure the message appears

        % Generate wave
        [t, eta] = generateWaveFromGUI();

        % Close "Please wait" box
        if isvalid(hWait), close(hWait); end

        % Show "success" box that auto-closes in 5 seconds
        hSuccess = msgbox('Wave generation complete!','Success');
        pause(5);
        if isvalid(hSuccess), close(hSuccess); end

        % Plot the result
        axes(handles.ax); %#ok<LAXES>
        cla(handles.ax);
        plot(handles.ax, t, eta, 'b-','LineWidth',1.2);
        grid(handles.ax,'on');
        xlabel(handles.ax,'Time [s]');
        ylabel(handles.ax,'\eta(t) [m]');
        title(handles.ax,'Preview of Wave Elevation');
    end

    function onGenerateASCII(~,~)
        % Show "Please wait" box
        hWait = msgbox('Generating wave. Please wait...','Please wait','help');
        drawnow;

        [t, eta] = generateWaveFromGUI();

        if isvalid(hWait), close(hWait); end

        % Show success (auto-close)
        hSuccess = msgbox('Wave generation complete!','Success');
        pause(5);
        if isvalid(hSuccess), close(hSuccess); end

        % Now ask user for file name
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

    function onSaveWaveData(~,~)
        % Show "Please wait" box
        hWait = msgbox('Generating wave. Please wait...','Please wait','help');
        drawnow;

        [t, eta] = generateWaveFromGUI();

        if isvalid(hWait), close(hWait); end

        % Show success (auto-close)
        hSuccess = msgbox('Wave generation complete!','Success');
        pause(5);
        if isvalid(hSuccess), close(hSuccess); end

        % Then prompt for file
        cfg = getAllParams();
        [fileName, pathName] = uiputfile('*.mat','Save Wave Data');
        if fileName == 0, return; end

        waveData.time   = t;
        waveData.eta    = eta;
        waveData.config = cfg;
        save(fullfile(pathName, fileName),'waveData');
        msgbox(['Saved wave data + config to ', fullfile(pathName, fileName)],...
            'Wave Data Saved','modal');
    end

%% ====================================================================
%  Enabling/Disabling Controls
%% ====================================================================
    function updateUIEnable()
        handles = get(fig,'UserData');
        waveTypeVal = get(handles.waveTypePopup,'Value');
        isRegular   = (waveTypeVal == 1); 

        % Regular wave fields
        set([handles.edtH, handles.edtT, handles.edtBeta], 'Enable', bool2OnOff(isRegular));

        % Irregular wave fields
        set(handles.spectrumPopup, 'Enable', bool2OnOff(~isRegular));

        if ~isRegular
            modelVal = get(handles.spectrumPopup,'Value');
            modelStr = handles.modelList{modelVal};
            relevantFields = {};
            switch modelStr
                case 'JONSWAP 2-param'
                    relevantFields = {'edtSIWAHE','edtTpeak'};
                case 'JONSWAP 3-param'
                    relevantFields = {'edtSIWAHE','edtTpeak','edtGamma'};
                case 'JONSWAP 6-param'
                    relevantFields = {'edtOmegaP','edtAlpha','edtBetaVal','edtGamma','edtSigA','edtSigB'};
                % case 'Torsethaugen'
                %     relevantFields = {'edtSIWAHE','edtTpeak'};
                % case 'Ochi-Hubble'
                %     relevantFields = {'edtSIWAHE'};
                case 'PM (1-param)'
                    relevantFields = {'edtSIWAHE'};
                case 'PM (2-param, Tavg)'
                    relevantFields = {'edtSIWAHE','edtTavg'};
                case 'PM (2-param, Tzc)'
                    relevantFields = {'edtSIWAHE','edtTzc'};
            end

            allIrrHandles = {'edtSIWAHE','edtTpeak','edtGamma','edtOmegaP','edtAlpha',...
                'edtBetaVal','edtSigA','edtSigB','edtTavg','edtTzc'};

            for ii=1:numel(allIrrHandles)
                set(handles.(allIrrHandles{ii}), 'Enable','off');
            end
            for ii=1:numel(relevantFields)
                set(handles.(relevantFields{ii}), 'Enable','on');
            end
        else
            allIrrHandles = {'edtSIWAHE','edtTpeak','edtGamma','edtOmegaP','edtAlpha',...
                'edtBetaVal','edtSigA','edtSigB','edtTavg','edtTzc'};
            for ii=1:numel(allIrrHandles)
                set(handles.(allIrrHandles{ii}), 'Enable','off');
            end
        end
    end

    function str = bool2OnOff(tf)
        if tf
            str = 'on';
        else
            str = 'off';
        end
    end

%% ====================================================================
%  Get/Set Parameter Struct
%% ====================================================================
    function cfg = getAllParams()
        handles = get(fig,'UserData');
        cfg.waveTypeVal = get(handles.waveTypePopup,'Value');
        cfg.spectrumVal = get(handles.spectrumPopup,'Value');

        cfg.tStep   = str2double(get(handles.edtTimeStep,'String'));
        cfg.nSamples= str2double(get(handles.edtNumSamples,'String'));
        cfg.tRamp   = str2double(get(handles.edtRamp,'String'));
        cfg.rSeed   = str2double(get(handles.edtRandomSeed,'String'));

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
        handles = get(fig,'UserData');

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
        cfg = getAllParams();
        dt  = cfg.tStep;
        N   = cfg.nSamples;
        t   = (0:N-1).'*dt;

        switch cfg.waveTypeVal
            case 1 % Regular wave
                omega = 2*pi/cfg.T;
                rampFactor = rampUp(t, cfg.tRamp);
                eta = 0.5*cfg.H * cos(omega*t + cfg.beta) .* rampFactor;
            otherwise
                rng(cfg.rSeed);
                modelStr = handles.modelList{cfg.spectrumVal};
                switch modelStr
                    case 'JONSWAP 2-param'
                        specFun = @(w) jonswap2_skipZero(w, cfg.Hs, cfg.Tpeak);
                    case 'JONSWAP 3-param'
                        specFun = @(w) jonswap3_skipZero(w, cfg.Hs, cfg.Tpeak, cfg.gamma);
                    case 'JONSWAP 6-param'
                        specFun = @(w) jonswap6(w, cfg.omegaP, cfg.alpha, cfg.betaVal, cfg.gamma, cfg.sigA, cfg.sigB);
                    % case 'Torsethaugen'
                    %     specFun = @(w) torsethaugenSpectrum_skipZero(w, cfg.Hs, cfg.Tpeak);
                    % case 'Ochi-Hubble'
                    %     specFun = @(w) ochiHubbleSpectrum_skipZero(w, cfg.Hs);
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
r = ones(size(t));
idx = (t < t_ramp);
r(idx) = t(idx)/t_ramp;
end

function eta = generateIrregularWave(t, t_ramp, specFun)
N = length(t);
dt = t(2)-t(1);
df = 1/(N*dt);
f  = (0:(N-1)).'*df;
w  = 2*pi*f;

S = specFun(w);
S(w<=0) = 0;  % Hard‐zero at w <= 0, just in case

dw = 2*pi*df;
A  = sqrt(2*S.*dw);
phi = 2*pi*rand(size(w));

Npos = floor(N/2);
eta_temp = zeros(size(t));
for i = 2:Npos
    eta_temp = eta_temp + A(i)*cos(w(i)*t + phi(i));
end

rampFactor = rampUp(t, t_ramp);
eta = rampFactor .* eta_temp;
end

%% -- JONSWAP 2/3 with zero-frequency skip --
function S = jonswap2_skipZero(w, Hs, Tp)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
gammaVal = 3.3;
sigma = 0.07.*(w_<=wp) + 0.09.*(w_>wp);

S_temp = alpha_PM ./ (w_.^5) .* exp(-5/4*(wp./w_).^4);
exponent = exp(- ((w_ - wp).^2)./(2*sigma.^2*wp.^2));
S_temp = S_temp .* (gammaVal.^exponent);

S(idx) = scaleSpectrumToHs(S_temp, w_, Hs);
end

function S = jonswap3_skipZero(w, Hs, Tp, gammaVal)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
sigma = 0.07.*(w_<=wp) + 0.09.*(w_>wp);

S_temp = alpha_PM./(w_.^5).*exp(-5/4*(wp./w_).^4);
exponent = exp(-((w_-wp).^2)./(2*sigma.^2.*wp.^2));
S_temp = S_temp.*(gammaVal.^exponent);

S(idx) = scaleSpectrumToHs(S_temp, w_, Hs);
end

function S = torsethaugenSpectrum(w, Hs, Tp)
% torsethaugenSpectrum
% A simplified double-peaked Torsethaugen-type spectrum.
%
% Inputs:
%   w  : angular frequencies (vector), w >= 0
%   Hs : total significant wave height
%   Tp : (user-chosen) peak period of primary system
%
% Outputs:
%   S  : spectral density [m^2·s/rad] for each w
%
% References:
%  - Torsethaugen (2004), "Simplified Double Peak Spectrum" approach
%  - Paper No. 2004-JSC-193

% ---- 1) Decide how to split Hs between two peaks
% For demonstration, let's do 80% energy in primary peak, 20% in secondary
Hs1 = 0.8 * Hs; 
Hs2 = sqrt(Hs^2 - Hs1^2);  % so that Hs^2 = Hs1^2 + Hs2^2

% ---- 2) Assign periods
Tp1 = Tp;        % primary system
deltaOffset = 3; % e.g. 3-second offset for secondary peak
Tp2 = Tp + deltaOffset;

% ---- 3) Peakedness factors
gamma1 = 3.3;   % typical for primary
gamma2 = 1.0;   % secondary often set to 1 (no peak enhancement)

% ---- 4) Build each partial spectrum with a JONSWAP-like formula
S1 = jonswapPartial(w, Hs1, Tp1, gamma1);
S2 = jonswapPartial(w, Hs2, Tp2, gamma2);

% ---- 5) Sum to get final double-peaked spectrum
S = S1 + S2;

end

function Sj = jonswapPartial(w, Hs_j, Tp_j, gamma_j)
% jonswapPartial: a JONSWAP-like partial spectrum, scaled to (Hs_j, Tp_j).
g = 9.81;
Sj = zeros(size(w));
idx = (w > 0);
w_ = w(idx);

wp_j = 2*pi / Tp_j;
alpha_PM = 0.0081*g^2;  
sigma = 0.07.*(w_<=wp_j) + 0.09.*(w_>wp_j);

% Pierson-Moskowitz core
S_base = alpha_PM ./ (w_.^5) .* exp(-5/4*(wp_j./w_).^4);
% JONSWAP peak enhancement
exponent = exp(-((w_-wp_j).^2)./(2*sigma.^2.*wp_j.^2));
S_g = S_base .* (gamma_j.^exponent);

% Scale to match partial Hs_j
S_g = scaleSpectrumToHs(S_g, w_, Hs_j);

% Insert back into original vector
Sj(idx) = S_g;
end

%% ------------------- JONSWAP 6-param, PM variants (unchanged) -----------
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
end

function S = pmSpectrum1(w, Hs)
g = 9.81;
alpha = 0.0081;
% guess a peak freq:
w_p = 0.44*g/((Hs/4)^1.0);
S = alpha*g^2./(w.^5).*exp(-1.25*(w_p./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S,w,Hs);
end

function S = pmSpectrum2_AVWAPE(w, Hs, Tavg)
g = 9.81;
w_ave = 2*pi/Tavg;
alpha = 0.0081;
S = alpha*g^2./(w.^5).*exp(-1.25*(w_ave./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S,w,Hs);
end

function S = pmSpectrum2_ZCWAPE(w, Hs, Tzc)
g = 9.81;
w_zc = 2*pi/Tzc;
alpha = 0.0081;
S = alpha*g^2./(w.^5).*exp(-1.25*(w_zc./w).^4);
S(w<=0) = 0;
S = scaleSpectrumToHs(S,w,Hs);
end

%% ------------------------------------------------------------------------
%  scaleSpectrumToHs for sub-spectra that skip zero freq
%% ------------------------------------------------------------------------
function Sscaled = scaleSpectrumToHs(S_partial, w_partial, Hs)
% S_partial is already zero where w<=0. w_partial = w(idx>0).
% We do a numeric integral of S_partial over w_partial:
dw = mean(diff(w_partial));
m0 = trapz(w_partial, S_partial);
Hs_model = 4*sqrt(m0);
if Hs_model == 0
    Sscaled = S_partial; % edge case: no scaling
else
    scaleFactor = (Hs / Hs_model)^2;
    Sscaled = S_partial * scaleFactor;
end
end

%% ========================================================================
%  Write ASCII file for SIMA
%% ========================================================================
function writeSIMAAscii(outFileName, eta, t_step)
fid = fopen(outFileName, 'w');
if fid < 0
    error('Could not open file: %s', outFileName);
end

N = length(eta);

fprintf(fid, '%d\n', N);            % Number of samples
fprintf(fid, '%.4f\n', t_step);    % Time step
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

