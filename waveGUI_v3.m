function waveGUI_v3()
% waveGUI_v3 - Enhanced Wave Generator GUI with Torsethaugen & Ochi-Hubble
%
% Description:
%   A MATLAB GUI for generating both regular and irregular wave elevation signals,
%   plotting the time-domain wave and its theoretical power spectrum, saving wave data,
%   and generating SIMA-compatible ASCII files. The tool incorporates a ramp function
%   to gradually introduce wave amplitudes for a smooth simulation start.
%
% Features:
%   - Supports generating Regular waves and various Irregular waves.
%   - Multiple spectral models available:
%       * JONSWAP (2-param, 3-param, 6-param)
%       * Torsethaugen (enhanced implementation)
%       * Ochi-Hubble (enhanced implementation)
%       * PM variants (1-param, 2-param with Tavg or Tzc)
%   - Displays all parameters in one window with an improved, intuitive layout.
%   - Dynamically disables irrelevant input fields based on the selected wave/spectrum model.
%   - Default window size set to 1024 x 666 pixels.
%   - Includes load/save configuration options.
%   - Provides feedback using a "please wait" popup during wave generation along with an
%     auto‐closing success message.
%
% Version Log:
%   v3.0 - 2025.03.07 : Integrated updated Torsethaugen and Ochi-Hubble implementations.
%   v2.1 - 2025.03.06 : Updated simulation parameter controls and layout adjustments.
%   v2.0 - 2025.03.06 : Introduced both regular and irregular wave generation methods with improved logic.
%   v1.1 - 2025.03.06 : Added popups for wave generation and fixed skip-zero issues in spectrum calculations.
%   v1.0 - 2025.03.06 : Initial version.
%
% Author: Shuijin Li
% Email: lishuijin@nbu.edu.cn
%

%% --------------------------------------------------------------------
%                     Figure & Panel Setup
%% --------------------------------------------------------------------
figWidth  = 1024;
figHeight = 666;
fig = figure('Name','Wave Generator GUI (Torsethaugen & Ochi-Hubble)',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','figure',...
    'Resize','off',...
    'Position',[100 100 figWidth figHeight]);

panelWidth  = 380;
panelHeight = figHeight - 20;
ctrlPanel = uipanel('Parent',fig,...
    'Units','pixels',...
    'Title','User Inputs',...
    'FontWeight','bold',...
    'Position',[10 10 panelWidth panelHeight]);

% Split the plotting area into top (spectrum) and bottom (time-domain)
plotPanelHeight = figHeight - 120;
axWidth  = figWidth - (panelWidth + 70);
axHeight = floor(plotPanelHeight / 2) - 40;

% Top axis: spectrum
axSpec = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[panelWidth+60, 80 + axHeight + 80, axWidth, axHeight]);
title(axSpec, 'Wave Spectrum');
xlabel(axSpec, '\omega [rad/s]');
ylabel(axSpec, 'S(\omega) [m^2·s]');
grid(axSpec,'on');

% Bottom axis: time-domain wave
axTime = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[panelWidth+60, 80, axWidth, axHeight]);
title(axTime, 'Wave Elevation Preview');
xlabel(axTime, 'Time [s]');
ylabel(axTime, '\eta(t) [m]');
grid(axTime,'on');

% A status label at the bottom
statusLabel = uicontrol('Parent',fig,...
    'Style','text',...
    'String','Status: Ready',...
    'HorizontalAlignment','left',...
    'Position',[panelWidth+60, 10, axWidth, 20]);

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

% Full set of possible spectrum models (including Torsethaugen & Ochi-Hubble)
modelList = {'JONSWAP 2-param','JONSWAP 3-param','JONSWAP 6-param',...
    'Torsethaugen (testing)','Ochi-Hubble (testing)',...
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
%         5) Irregular Wave Parameters
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

btnPlot = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'Position',[xLabel yBottom btnWidth+20 btnHeight],...
    'String','Plot Wave','FontWeight','bold',...
    'Callback',@onPlotWave);

btnSaveData = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Save Wave Data',...
    'Position',[xLabel+btnWidth+spacing+20 yBottom btnWidth btnHeight],...
    'Callback',@onSaveWaveData);

btnGen = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Generate ASCII',...
    'Position',[xLabel+(btnWidth+spacing+20)*1.85 yBottom btnWidth+20 btnHeight],...
    'Callback',@onGenerateASCII);

%% --------------------------------------------------------------------
%             Store Handles in the Figure's UserData
%% --------------------------------------------------------------------
handles.axSpec = axSpec;
handles.axTime = axTime;
handles.statusLabel = statusLabel;

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

% For caching the last wave data & spectrum to avoid repeated generation
handles.lastWave = [];
handles.lastSpec = [];

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
            clearWaveCache();
            setStatus('Configuration loaded. Old wave data cleared.');
        end
    end

    function onSaveConfig(~,~)
        [fileName, pathName] = uiputfile('*.mat','Save Configuration');
        if fileName == 0, return; end
        cfg = getAllParams();
        save(fullfile(pathName,fileName),'cfg');
        setStatus('Configuration saved successfully.');
    end

    function onWaveTypeChanged(~,~)
        updateUIEnable();
        clearWaveCache();
        setStatus('Wave type changed. Old wave data cleared.');
    end

    function onSpectrumModelChanged(~,~)
        updateUIEnable();
        clearWaveCache();
        setStatus('Spectrum model changed. Old wave data cleared.');
    end

    function onPlotWave(~,~)
        [t, eta] = getOrGenerateWaveData();
        [wPlot, SPlot] = getOrComputeSpectrum();

        % Plot time-domain wave
        axes(handles.axTime); %#ok<LAXES>
        cla(handles.axTime);
        plot(handles.axTime, t, eta, 'b-','LineWidth',1.2);
        grid(handles.axTime,'on');
        xlabel(handles.axTime,'Time [s]');
        ylabel(handles.axTime,'\eta(t) [m]');
        title(handles.axTime,'Preview of Wave Elevation');

        % Plot spectrum
        axes(handles.axSpec); %#ok<LAXES>
        cla(handles.axSpec);
        plot(handles.axSpec, wPlot, SPlot, 'r-','LineWidth',1.2);
        grid(handles.axSpec,'on');
        xlabel(handles.axSpec,'\omega [rad/s]');
        ylabel(handles.axSpec,'S(\omega) [m^2·s]');
        title(handles.axSpec,'Wave Spectrum');

        setStatus('Plot complete.');
    end

    function onSaveWaveData(~,~)
        [t, eta] = getOrGenerateWaveData();
        cfg = getAllParams();

        [fileName, pathName] = uiputfile('*.mat','Save Wave Data');
        if fileName == 0, return; end

        waveData.time   = t;
        waveData.eta    = eta;
        waveData.config = cfg;
        save(fullfile(pathName, fileName),'waveData');
        setStatus(['Saved wave data + config to ', fullfile(pathName, fileName)]);
    end

    function onGenerateASCII(~,~)
        [t, eta] = getOrGenerateWaveData();
        getOrComputeSpectrum();  % ensures the spectrum is also up to date (if desired)

        [fileName, pathName] = uiputfile('waveElevation_sima.dat','Generate SIMA ASCII');
        if fileName == 0, return; end
        outFileName = fullfile(pathName,fileName);

        dt = t(2) - t(1);
        writeSIMAAscii(outFileName, eta, dt);
        setStatus(sprintf('Wrote wave time series to %s', outFileName));
    end

%% ====================================================================
%  Helper: Clear cached wave data
%% ====================================================================
    function clearWaveCache()
        handles = get(fig,'UserData');
        handles.lastWave = [];
        handles.lastSpec = [];
        set(fig,'UserData',handles);
    end

%% ====================================================================
%  Helper: Set status message
%% ====================================================================
    function setStatus(msg)
        handles = get(fig,'UserData');
        set(handles.statusLabel, 'String', ['Status: ', msg]);
        drawnow;
    end

%% ====================================================================
%  Get or Generate Wave Data (Avoid Re-computation if Possible)
%% ====================================================================
    function [t, eta] = getOrGenerateWaveData()
        handlesLocal = get(fig,'UserData');
        cfg = getAllParams();

        if ~isempty(handlesLocal.lastWave) && cfgStructsEqual(handlesLocal.lastWave.cfg, cfg)
            t   = handlesLocal.lastWave.t;
            eta = handlesLocal.lastWave.eta;
            setStatus('Reusing cached wave data.');
        else
            setStatus('Generating wave data...');
            [t, eta] = generateWaveFromGUI(cfg);
            handlesLocal.lastWave.t   = t;
            handlesLocal.lastWave.eta = eta;
            handlesLocal.lastWave.cfg = cfg;
            set(fig,'UserData',handlesLocal);
            setStatus('Wave generation complete.');
        end
    end

%% ====================================================================
%  Get or Compute Theoretical Spectrum (Avoid Re-computation)
%% ====================================================================
    function [wPlot, SPlot] = getOrComputeSpectrum()
        handlesLocal = get(fig,'UserData');
        cfg = getAllParams();

        if ~isempty(handlesLocal.lastSpec) && cfgStructsEqual(handlesLocal.lastSpec.cfg, cfg)
            wPlot = handlesLocal.lastSpec.w;
            SPlot = handlesLocal.lastSpec.S;
            setStatus('Reusing cached spectrum.');
        else
            setStatus('Computing wave spectrum...');
            [wPlot, SPlot] = computeSpectrumFromCfg(cfg);
            handlesLocal.lastSpec.w = wPlot;
            handlesLocal.lastSpec.S = SPlot;
            handlesLocal.lastSpec.cfg = cfg;
            set(fig,'UserData',handlesLocal);
            setStatus('Spectrum computation complete.');
        end
    end

%% ====================================================================
%  Check if two configs are identical
%% ====================================================================
    function tf = cfgStructsEqual(cfg1, cfg2)
        fA = sort(fieldnames(cfg1));
        fB = sort(fieldnames(cfg2));
        if ~isequal(fA,fB)
            tf = false; return;
        end
        for fn = 1:numel(fA)
            if ~isequal(cfg1.(fA{fn}), cfg2.(fA{fn}))
                tf = false; return;
            end
        end
        tf = true;
    end

%% ====================================================================
%  Enabling/Disabling Controls
%% ====================================================================
    function updateUIEnable()
        handlesLocal = get(fig,'UserData');
        waveTypeVal = get(handlesLocal.waveTypePopup,'Value');
        isRegular   = (waveTypeVal == 1);

        % Regular wave fields
        set([handlesLocal.edtH, handlesLocal.edtT, handlesLocal.edtBeta], 'Enable', bool2OnOff(isRegular));

        % Irregular wave fields
        set(handlesLocal.spectrumPopup, 'Enable', bool2OnOff(~isRegular));

        if ~isRegular
            modelVal = get(handlesLocal.spectrumPopup,'Value');
            modelStr = handlesLocal.modelList{modelVal};
            relevantFields = {};
            switch modelStr
                case 'JONSWAP 2-param'
                    relevantFields = {'edtSIWAHE','edtTpeak'};
                case 'JONSWAP 3-param'
                    relevantFields = {'edtSIWAHE','edtTpeak','edtGamma'};
                case 'JONSWAP 6-param'
                    relevantFields = {'edtOmegaP','edtAlpha','edtBetaVal','edtGamma','edtSigA','edtSigB'};
                case 'Torsethaugen'
                    % We will assume user only needs: Hs (SIWAHE) & Tpeak
                    % (Similar to the original placeholders).
                    relevantFields = {'edtSIWAHE','edtTpeak'};
                case 'Ochi-Hubble'
                    % Minimal user input: just Hs (SIWAHE).
                    % You could add more fields if you want 6 param inputs.
                    relevantFields = {'edtSIWAHE'};
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
                set(handlesLocal.(allIrrHandles{ii}), 'Enable','off');
            end
            for ii=1:numel(relevantFields)
                set(handlesLocal.(relevantFields{ii}), 'Enable','on');
            end
        else
            % If wave type is Regular
            allIrrHandles = {'edtSIWAHE','edtTpeak','edtGamma','edtOmegaP','edtAlpha',...
                'edtBetaVal','edtSigA','edtSigB','edtTavg','edtTzc'};
            for ii=1:numel(allIrrHandles)
                set(handlesLocal.(allIrrHandles{ii}), 'Enable','off');
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
        handlesLocal = get(fig,'UserData');
        cfg.waveTypeVal = get(handlesLocal.waveTypePopup,'Value');
        cfg.spectrumVal = get(handlesLocal.spectrumPopup,'Value');

        cfg.tStep   = str2double(get(handlesLocal.edtTimeStep,'String'));
        cfg.nSamples= str2double(get(handlesLocal.edtNumSamples,'String'));
        cfg.tRamp   = str2double(get(handlesLocal.edtRamp,'String'));
        cfg.rSeed   = str2double(get(handlesLocal.edtRandomSeed,'String'));

        cfg.H       = str2double(get(handlesLocal.edtH,'String'));
        cfg.T       = str2double(get(handlesLocal.edtT,'String'));
        cfg.beta    = str2double(get(handlesLocal.edtBeta,'String'));

        cfg.Hs      = str2double(get(handlesLocal.edtSIWAHE,'String'));
        cfg.Tpeak   = str2double(get(handlesLocal.edtTpeak,'String'));
        cfg.gamma   = str2double(get(handlesLocal.edtGamma,'String'));
        cfg.omegaP  = str2double(get(handlesLocal.edtOmegaP,'String'));
        cfg.alpha   = str2double(get(handlesLocal.edtAlpha,'String'));
        cfg.betaVal = str2double(get(handlesLocal.edtBetaVal,'String'));
        cfg.sigA    = str2double(get(handlesLocal.edtSigA,'String'));
        cfg.sigB    = str2double(get(handlesLocal.edtSigB,'String'));
        cfg.Tavg    = str2double(get(handlesLocal.edtTavg,'String'));
        cfg.Tzc     = str2double(get(handlesLocal.edtTzc,'String'));
    end

    function setParams(cfg)
        handlesLocal = get(fig,'UserData');

        set(handlesLocal.waveTypePopup,'Value', cfg.waveTypeVal);
        set(handlesLocal.spectrumPopup,'Value', cfg.spectrumVal);

        set(handlesLocal.edtTimeStep,'String', num2str(cfg.tStep));
        set(handlesLocal.edtNumSamples,'String', num2str(cfg.nSamples));
        set(handlesLocal.edtRamp,'String', num2str(cfg.tRamp));
        set(handlesLocal.edtRandomSeed,'String', num2str(cfg.rSeed));

        set(handlesLocal.edtH,'String', num2str(cfg.H));
        set(handlesLocal.edtT,'String', num2str(cfg.T));
        set(handlesLocal.edtBeta,'String', num2str(cfg.beta));

        set(handlesLocal.edtSIWAHE,'String', num2str(cfg.Hs));
        set(handlesLocal.edtTpeak,'String', num2str(cfg.Tpeak));
        set(handlesLocal.edtGamma,'String', num2str(cfg.gamma));
        set(handlesLocal.edtOmegaP,'String', num2str(cfg.omegaP));
        set(handlesLocal.edtAlpha,'String', num2str(cfg.alpha));
        set(handlesLocal.edtBetaVal,'String', num2str(cfg.betaVal));
        set(handlesLocal.edtSigA,'String', num2str(cfg.sigA));
        set(handlesLocal.edtSigB,'String', num2str(cfg.sigB));
        set(handlesLocal.edtTavg,'String', num2str(cfg.Tavg));
        set(handlesLocal.edtTzc,'String', num2str(cfg.Tzc));
    end

%% ====================================================================
%  Generate Wave According to GUI Parameters
%% ====================================================================
    function [t, eta] = generateWaveFromGUI(cfg)
        dt = cfg.tStep;
        N  = cfg.nSamples;
        t  = (0:N-1).'*dt;

        if cfg.waveTypeVal == 1
            % Regular wave
            omega = 2*pi/cfg.T;
            rampFactor = rampUp(t, cfg.tRamp);
            eta = 0.5*cfg.H * sin(omega*t + cfg.beta) .* rampFactor;
        else
            % Irregular wave
            rng(cfg.rSeed);
            modelStr = handles.modelList{cfg.spectrumVal};
            switch modelStr
                case 'JONSWAP 2-param'
                    specFun = @(w) jonswap2_skipZero(w, cfg.Hs, cfg.Tpeak);
                case 'JONSWAP 3-param'
                    specFun = @(w) jonswap3_skipZero(w, cfg.Hs, cfg.Tpeak, cfg.gamma);
                case 'JONSWAP 6-param'
                    specFun = @(w) jonswap6(w, cfg.omegaP, cfg.alpha, cfg.betaVal, cfg.gamma, cfg.sigA, cfg.sigB);
                case 'Torsethaugen'
                    specFun = @(w) torsethaugenSpectrum_skipZero(w, cfg.Hs, cfg.Tpeak);
                case 'Ochi-Hubble'
                    specFun = @(w) ochiHubbleSpectrum_skipZero(w, cfg.Hs);
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

%% ====================================================================
%  Compute Theoretical Spectrum for Current Config
%% ====================================================================
    function [wPlot, SPlot] = computeSpectrumFromCfg(cfg)
        if cfg.waveTypeVal == 1
            % Regular wave -> approximate a delta spike near w0 = 2*pi/T
            w0 = 2*pi/cfg.T;
            wMax = max(2.5*w0, 2.5);
            dw = 0.01;
            wPlot = (0:dw:wMax).';
            SPlot = zeros(size(wPlot));
            idxPeak = round(w0/dw)+1;
            if idxPeak>0 && idxPeak<=numel(SPlot)
                % approximate delta function integral of H^2/8 (peak area)
                SPlot(idxPeak) = (cfg.H^2/8)/dw;
            end
        else
            % Irregular wave -> typical range from 0..some max
            modelStr = handles.modelList{cfg.spectrumVal};
            wMax = 5*(2*pi/cfg.Tpeak); % a default guess for max freq
            if strcmp(modelStr,'JONSWAP 6-param') && (cfg.omegaP>0)
                wMax = max(wMax, 4*cfg.omegaP);
            elseif strcmp(modelStr,'PM (2-param, Tavg)') && (cfg.Tavg>0)
                wMax = max(wMax, 5*(2*pi/cfg.Tavg));
            elseif strcmp(modelStr,'PM (2-param, Tzc)') && (cfg.Tzc>0)
                wMax = max(wMax, 5*(2*pi/cfg.Tzc));
            end
            if wMax<2.5, wMax=2.5; end
            dw = 0.01;
            wPlot = (0:dw:wMax).';

            switch modelStr
                case 'JONSWAP 2-param'
                    sfun = @(w) jonswap2_skipZero(w, cfg.Hs, cfg.Tpeak);
                case 'JONSWAP 3-param'
                    sfun = @(w) jonswap3_skipZero(w, cfg.Hs, cfg.Tpeak, cfg.gamma);
                case 'JONSWAP 6-param'
                    sfun = @(w) jonswap6(w, cfg.omegaP, cfg.alpha, cfg.betaVal, cfg.gamma, cfg.sigA, cfg.sigB);
                case 'Torsethaugen'
                    sfun = @(w) torsethaugenSpectrum_skipZero(w, cfg.Hs, cfg.Tpeak);
                case 'Ochi-Hubble'
                    sfun = @(w) ochiHubbleSpectrum_skipZero(w, cfg.Hs);
                case 'PM (1-param)'
                    sfun = @(w) pmSpectrum1(w, cfg.Hs);
                case 'PM (2-param, Tavg)'
                    sfun = @(w) pmSpectrum2_AVWAPE(w, cfg.Hs, cfg.Tavg);
                case 'PM (2-param, Tzc)'
                    sfun = @(w) pmSpectrum2_ZCWAPE(w, cfg.Hs, cfg.Tzc);
            end
            SPlot = sfun(wPlot);
            SPlot(wPlot<=0) = 0;
        end
    end

end % waveGUI_v3 function

%% ========================================================================
%  Utility sub-functions
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
S(w<=0) = 0;  % Hard-zero at w<=0

dw = 2*pi*df;
A  = sqrt(2*S.*dw);
phi = 2*pi*rand(size(w));

Npos = floor(N/2);
eta_temp = zeros(size(t));
for i = 2:Npos
    eta_temp = eta_temp + A(i)*sin(w(i)*t + phi(i));
end

rampFactor = rampUp(t, t_ramp);
eta = rampFactor .* eta_temp;
end

%% ------------------- JONSWAP 2/3 skipZero & JONSWAP 6-param -----------
function S = jonswap2_skipZero(w, Hs, Tp)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
gammaVal = 3.3;
sigma = 0.07.*(w_<=wp) + 0.09.*(w_>wp);

S_tmp = alpha_PM ./ (w_.^5).*exp(-5/4*(wp./w_).^4);
expnt = exp(-((w_-wp).^2)./(2*sigma.^2.*wp.^2));
S_tmp = S_tmp .* (gammaVal.^expnt);

S(idx) = scaleSpectrumToHs(S_tmp, w_, Hs);
end

function S = jonswap3_skipZero(w, Hs, Tp, gammaVal)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

wp = 2*pi/Tp;
alpha_PM = 0.0081*g^2;
sigma = 0.07.*(w_<=wp) + 0.09.*(w_>wp);

S_tmp = alpha_PM./(w_.^5).*exp(-5/4*(wp./w_).^4);
expnt = exp(-((w_-wp).^2)./(2*sigma.^2.*wp.^2));
S_tmp = S_tmp.*(gammaVal.^expnt);

S(idx) = scaleSpectrumToHs(S_tmp, w_, Hs);
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

S(idx) = S_PM.*G;
end

%% ------------------- Torsethaugen Implementation ----------------------
function S = torsethaugenSpectrum_skipZero(w, Hs, Tp)
% Convert w->f, compute Torsethaugen in [m^2/Hz], then convert to [m^2*s/rad]
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);
f_ = w_/(2*pi);

S_hz = torsethaugenSpectrum(f_, Hs, Tp);  % in m^2/Hz
% Convert to m^2 s/rad: multiply by 2*pi
S_rad = S_hz .* (2*pi);

S(idx) = S_rad;
end

function S = torsethaugenSpectrum(f, Hs, Tp)
% A simplified Torsethaugen approach, focusing on Hs, Tp inputs
% Double-peaked combination of a wind-sea and a swell system.

f = f(:);
% If needed, define how we partition wind vs. swell.
% For demonstration, we do a straightforward approach:
[S1, S2] = torsetComponents(f, Hs, Tp);
S = S1 + S2;
end

function [S1, S2] = torsetComponents(f, Hs, Tp)
% Demo partition: wind-sea dominated if Tp < Tpf, else swell-dominated
% We'll define Tpf ~ 4.9*(Hs^(1/3)) (simple approach),
% or you can refine with actual Torsethaugen references.

Tpf = 4.9 * (Hs^(1/3));
if Tp <= Tpf
    % wind-sea dominated
    % 80% of variance in primary, 20% in secondary
    Hswind = 0.8 * Hs;
    Hsswell= sqrt(Hs^2 - Hswind^2);
    Tp_wind = Tp;
    Tp_swell= Tpf + 2.0;
else
    % swell-dominated
    Hsswell = 0.7 * Hs;
    Hswind  = sqrt(Hs^2 - Hsswell^2);
    Tp_swell= Tp;
    Tp_wind = Tpf;
end

S1 = jonswapPartial(f, Hswind, Tp_wind, 3.3);  % treat wind-sea
S2 = jonswapPartial(f, Hsswell, Tp_swell, 1.0);% treat swell (gamma=1 => broader)
end

function Sj = jonswapPartial(f, Hs_j, Tp_j, gamma_j)
% Simple partial JONSWAP used inside Torsethaugen
g = 9.81;
f = f(:);

% alpha_PM ~ 0.0081 ...
% We'll do a direct scaling to match Hs_j
% JONSWAP shape
wp_j = 2*pi/Tp_j;
fp_j = 1/Tp_j;
sigma = 0.07.*(2*pi*f <= wp_j) + 0.09.*(2*pi*f > wp_j);

S_base = 0.0081*g^2 ./ ((2*pi*f).^5) .* exp(-5/4 * (wp_j ./ (2*pi*f)).^4);
peakEnh = gamma_j.^ (exp(- ( (2*pi*f - wp_j).^2 ) ./ (2*sigma.^2 .* wp_j.^2 )));
Sj_ = S_base .* peakEnh;

% scale to match Hs_j
Sj = scaleSpectrumToHs(Sj_, 2*pi*f, Hs_j);
end

%% ------------------- Ochi-Hubble Implementation -----------------------
function S = ochiHubbleSpectrum_skipZero(w, Hs)
% For demonstration, we'll define a simple 2-peak Ochi-Hubble with
% user-specified total Hs.
% Real usage might require specifying all 6 parameters.
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);
f_ = w_/(2*pi);

% Example: define 2 sets of (Hs_j, Tpj, lambda_j)
% Let's do a small wind sea plus a bigger swell, just as an example:
Hs1 = 0.6 * Hs;  Tp1 = 14;  lam1 = 2;
Hs2 = sqrt(Hs^2 - Hs1^2);
if Hs2 < 0.1, Hs2=0; end
Tp2 = 8;        lam2 = 3.3;

S_hz = ochiHubbleSpectrum(f_, Hs1, Tp1, lam1, Hs2, Tp2, lam2);
S_rad = S_hz.*(2*pi);
S(idx) = S_rad;
end

function S = ochiHubbleSpectrum(f, Hs1, Tp1, lam1, Hs2, Tp2, lam2)
% Ochi-Hubble two-term spectrum in [m^2/Hz].
f = f(:);
fp1 = 1/Tp1;
fp2 = 1/Tp2;

Gamma1 = gamma(lam1);
Gamma2 = gamma(lam2);

A1 = ((4*lam1 + 1)/4)^(lam1) * (fp1)^(4*lam1) * Hs1^2/(4 * Gamma1);
A2 = ((4*lam2 + 1)/4)^(lam2) * (fp2)^(4*lam2) * Hs2^2/(4 * Gamma2);

S1 = A1.* f.^(-(4*lam1 + 1)) .* exp(-(4*lam1 + 1)/4.* (fp1./f).^4);
S2 = A2.* f.^(-(4*lam2 + 1)) .* exp(-(4*lam2 + 1)/4.* (fp2./f).^4);

S = S1 + S2;
% Zero out negative/zero freq if any
S(f<=0)=0;
end

%% ------------------- PM variants, etc. -------------------------------
function S = pmSpectrum1(w, Hs)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

% guess a peak freq from Hs
A = Hs/4;
wp = 0.44*g/(A^1.0);
alpha = 0.0081;

S_tmp = alpha*g^2./(w_.^5).*exp(-1.25*(wp./w_).^4);
S(idx) = scaleSpectrumToHs(S_tmp, w_, Hs);
end

function S = pmSpectrum2_AVWAPE(w, Hs, Tavg)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

w_ave = 2*pi/Tavg;
alpha = 0.0081;
S_tmp = alpha*g^2./(w_.^5).*exp(-1.25*(w_ave./w_).^4);
S(idx) = scaleSpectrumToHs(S_tmp, w_, Hs);
end

function S = pmSpectrum2_ZCWAPE(w, Hs, Tzc)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

w_zc = 2*pi/Tzc;
alpha = 0.0081;
S_tmp = alpha*g^2./(w_.^5).*exp(-1.25*(w_zc./w_).^4);
S(idx) = scaleSpectrumToHs(S_tmp, w_, Hs);
end

%% ---------------------------------------------------------------------
%  scaleSpectrumToHs
%% ---------------------------------------------------------------------
function Sscaled = scaleSpectrumToHs(S_partial, w_partial, Hs)
dw = mean(diff(w_partial));
m0 = trapz(w_partial, S_partial);
Hs_model = 4*sqrt(m0);
if Hs_model == 0
    Sscaled = S_partial;
else
    scaleFactor = (Hs / Hs_model)^2;
    Sscaled = S_partial * scaleFactor;
end
end

%% ---------------------------------------------------------------------
%  Write ASCII for SIMA
%% ---------------------------------------------------------------------
function writeSIMAAscii(outFileName, eta, t_step)
fid = fopen(outFileName, 'w');
if fid < 0
    error('Could not open file: %s', outFileName);
end

N = length(eta);
fprintf(fid, '%d\n', N);
fprintf(fid, '%.4f\n', t_step);
fprintf(fid, 'Generated by waveGUI_v3\n');
fprintf(fid, 'Author: Shuijin LI; Email: lishuijin@nbu.edu.cn; Github: https://github.com/TaxLee\n');

valsPerLine = 1;
for i = 1:N
    fprintf(fid, '%12.6f', eta(i));
    if mod(i, valsPerLine)==0
        fprintf(fid,'\n');
    else
        fprintf(fid,'');
    end
end
if mod(N,valsPerLine)~=0
    fprintf(fid,'\n');
end
fclose(fid);
end
