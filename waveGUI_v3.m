function waveGUI_v3()
% waveGUI_v3_refactored - Enhanced Wave Generator GUI (Refactored)
%
% Description:
%   A MATLAB GUI for generating both regular and irregular wave elevation signals,
%   plotting the time-domain wave and its theoretical power spectrum, saving wave data,
%   and generating SIMA-compatible ASCII files. The tool incorporates a ramp function
%   to gradually introduce wave amplitudes for a smooth simulation start.
%
% What's New (v4.0 - 2025.03.14, Refactored):
%   - Code refactoring for improved maintainability.
%   - Core functionality unchanged, but reorganized into modular subfunctions.
%
% What's New (v3.2 - 2025.03.14):
%   - Displays measured maximum wave height from the generated time series,
%     plus a Rayleigh-based statistical estimate of maximum wave height (for irregular waves)
%     and the zero-crossing period Tz in the GUI.
%
% Features:
%   - Supports generating Regular waves and various Irregular waves.
%   - Multiple spectral models available:
%       * JONSWAP (2-param, 3-param, 6-param)
%       * Torsethaugen (testing)
%       * Ochi-Hubble (testing)
%       * PM variants (1-param, 2-param with Tavg or Tzc)
%   - Displays all parameters in one window with an improved, intuitive layout.
%   - Dynamically disables irrelevant input fields based on the selected wave/spectrum model.
%   - Default window size set to 1024 x 666 pixels.
%   - Includes load/save configuration options.
%   - Provides feedback using a "please wait" popup during wave generation along with an
%     auto‐closing success message.
%   - Updated ASCII file generator: The output file begins with the sample count and time step,
%     followed by a header containing tool info, a parameter description line, and then the wave data.
%
% Output File Format:
%   - Line 1: Number of data points (integer).
%   - Line 2: Time step (floating point with 4 decimal places).
%   - Line 3: Header comment line detailing generator info and author contact.
%   - Line 4: Parameter settings string.
%   - Subsequent lines: Wave elevation values (each printed in a fixed 12.6f format).
%
% Version Log:
%   v4.0 - 2025.03.14 : **Refactored** for easier maintenance.
%   v3.3 - 2025.03.14 : Reallocated results box, raised the lower plot, fixed Tz formula.
%   v3.2 - 2025.03.14 : Added real-time display of measured Hmax, estimated Hmax, and Tz.
%   v3.1 - 2025.03.10 : Updated ASCII file generator with expanded header and output file format.
%   v3.0 - 2025.03.07 : Integrated updated Torsethaugen and Ochi-Hubble implementations.
%   v2.1 - 2025.03.06 : Updated simulation parameter controls and layout adjustments.
%   v2.0 - 2025.03.06 : Introduced both regular and irregular wave generation methods.
%   v1.1 - 2025.03.06 : Added popups for wave generation and fixed skip-zero issues in spectrum calculations.
%   v1.0 - 2025.03.06 : Initial version.
%
% Author:  Shuijin Li
% Email:   lishuijin@nbu.edu.cn
%

% -------------------------------------------------------------------------
% REF NOTE: Build the GUI and store all handles in a struct.
% -------------------------------------------------------------------------
[fig, handles] = buildGUI();

% -------------------------------------------------------------------------
% REF NOTE: Set up callbacks for all UI controls. 
%           By passing 'handles' and 'fig' into nested subfunctions, 
%           we can keep code more organized.
% -------------------------------------------------------------------------
initializeCallbacks(handles, fig);

% -------------------------------------------------------------------------
% REF NOTE: Store the handles struct into the figure's UserData so that
%           all subfunctions can retrieve it easily.
% -------------------------------------------------------------------------
set(fig,'UserData',handles);

% -------------------------------------------------------------------------
% REF NOTE: One-time UI refresh to make sure fields are enabled/disabled
%           correctly at startup.
% -------------------------------------------------------------------------
updateUIEnable(fig);

end % waveGUI_v3_refactored

%% =========================================================================
%  SUBFUNCTION: BUILD THE GUI
%% =========================================================================
function [fig, handles] = buildGUI()
% buildGUI: Creates the figure and all UI controls, returning a structure
%           of handles and references to important GUI components.

figWidth  = 1024;
figHeight = 666;

fig = figure('Name','Wave Generator GUI (Refactored)',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','figure',...
    'Resize','off',...
    'Position',[100 100 figWidth figHeight]);

panelWidth  = 380;
panelHeight = figHeight - 20;
ctrlPanel   = uipanel('Parent',fig,...
    'Units','pixels',...
    'Title','User Inputs',...
    'FontWeight','bold',...
    'Position',[10 10 panelWidth panelHeight]);

% We'll leave some extra vertical space at the bottom for the wave metrics
plotPanelHeight = figHeight - 120;
axWidth  = figWidth - (panelWidth + 70);

% Adjusted: Raise the bottom plot from 80 -> 110
axHeight = floor(plotPanelHeight / 2) - 40;

% Top axis: spectrum
handles.axSpec = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[panelWidth+60, 80 + axHeight + 80, axWidth, axHeight]);
title(handles.axSpec, 'Wave Spectrum');
xlabel(handles.axSpec, '\omega [rad/s]');
ylabel(handles.axSpec, 'S(\omega) [m^2·s]');
grid(handles.axSpec,'on');

% Bottom axis: time-domain wave
handles.axTime = axes('Parent',fig,...
    'Units','pixels',...
    'Position',[panelWidth+60, 110, axWidth, axHeight]);
title(handles.axTime, 'Preview of Wave Elevation');
xlabel(handles.axTime, 'Time [s]');
ylabel(handles.axTime, '\eta(t) [m]');
grid(handles.axTime,'on');

% A status label near the bottom
handles.statusLabel = uicontrol('Parent',fig,...
    'Style','text',...
    'String','Status: Ready',...
    'HorizontalAlignment','left',...
    'Position',[panelWidth+60, 10, 220, 20]);

% Layout parameters
yTop = panelHeight - 40;
rowHeight = 20;
labelWidth = 200;
editWidth  = 60;
xLabel = 10;
xEdit  = xLabel + labelWidth + 10;

%% 1) Load/Save Config Buttons
handles.btnLoadCfg = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Load Config','Position',[xLabel yTop 90 rowHeight]);

handles.btnSaveCfg = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Save Config','Position',[xLabel+100 yTop 90 rowHeight]);

yTop = yTop - rowHeight - 10;

%% 2) Wave Type & Spectrum Model
uicontrol('Parent',ctrlPanel,'Style','text','String','Wave Type:',...
    'HorizontalAlignment','left',...
    'Position',[xLabel yTop labelWidth rowHeight]);
handles.waveTypePopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',{'Regular','Irregular'},...
    'Value',1,...
    'Position',[xEdit yTop editWidth+40 rowHeight]);

yTop = yTop - rowHeight - 5;

modelList = {'JONSWAP 2-param','JONSWAP 3-param','JONSWAP 6-param',...
    'Torsethaugen (testing)','Ochi-Hubble (testing)',...
    'PM (1-param)','PM (2-param, Tavg)','PM (2-param, Tzc)'};

uicontrol('Parent',ctrlPanel,'Style','text','String','Spectrum Model:',...
    'HorizontalAlignment','left',...
    'Position',[xLabel yTop labelWidth rowHeight]);
handles.spectrumPopup = uicontrol('Parent',ctrlPanel,'Style','popupmenu',...
    'String',modelList,...
    'Value',1,...
    'Position',[xEdit yTop editWidth+70 rowHeight]);

yTop = yTop - rowHeight - 15;

%% 3) Common Simulation Parameters
labels = {'Time Step (s):', 'Number of Samples:', 'Ramp Duration (sec):', 'Random Seed (Irregular):'};
defaultVals = {'0.1','40000','0','1'};
for ii = 1:numel(labels)
    uicontrol('Parent',ctrlPanel,'Style','text','String',labels{ii},...
        'HorizontalAlignment','left',...
        'Position',[xLabel yTop labelWidth rowHeight]);
    handles.edtCommon(ii) = uicontrol('Parent',ctrlPanel,'Style','edit',...
        'String',defaultVals{ii},...
        'Position',[xEdit yTop editWidth rowHeight]);
    yTop = yTop - rowHeight - 5;
end
% Unpack them more meaningfully
handles.edtTimeStep   = handles.edtCommon(1);
handles.edtNumSamples = handles.edtCommon(2);
handles.edtRamp       = handles.edtCommon(3);
handles.edtRandomSeed = handles.edtCommon(4);

yTop = yTop - 10;

%% 4) Regular Wave Parameters
regLabels = {'Wave Height, H:', 'Wave Period, T:', 'Phase, beta:'};
regDefault = {'9.96','7.99','0'};
for ii = 1:numel(regLabels)
    uicontrol('Parent',ctrlPanel,'Style','text','String',regLabels{ii},...
        'HorizontalAlignment','left',...
        'Position',[xLabel yTop labelWidth rowHeight]);
    handles.edtRegular(ii) = uicontrol('Parent',ctrlPanel,'Style','edit',...
        'String',regDefault{ii},...
        'Position',[xEdit yTop editWidth rowHeight]);
    yTop = yTop - rowHeight - 5;
end
handles.edtH    = handles.edtRegular(1);
handles.edtT    = handles.edtRegular(2);
handles.edtBeta = handles.edtRegular(3);

yTop = yTop - 15;

%% 5) Irregular Wave Parameters
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

irrDefaults = {'5.43','7.99','3.3','0.6283','0.0081','1.25','0.07','0.09','6.0','7.0'};
for ii = 1:numel(irrLabels)
    uicontrol('Parent', ctrlPanel, ...
        'Style', 'text', ...
        'String', irrLabels{ii}, ...
        'HorizontalAlignment', 'left', ...
        'Position', [xLabel, yTop, labelWidth, rowHeight]);

    handles.edtIrregular(ii) = uicontrol('Parent', ctrlPanel, ...
        'Style', 'edit', ...
        'String', irrDefaults{ii}, ...
        'Position', [xEdit, yTop, editWidth, rowHeight]);

    yTop = yTop - (rowHeight + 5);
end
% Unpack them more meaningfully
handles.edtSIWAHE  = handles.edtIrregular(1);
handles.edtTpeak   = handles.edtIrregular(2);
handles.edtGamma   = handles.edtIrregular(3);
handles.edtOmegaP  = handles.edtIrregular(4);
handles.edtAlpha   = handles.edtIrregular(5);
handles.edtBetaVal = handles.edtIrregular(6);
handles.edtSigA    = handles.edtIrregular(7);
handles.edtSigB    = handles.edtIrregular(8);
handles.edtTavg    = handles.edtIrregular(9);
handles.edtTzc     = handles.edtIrregular(10);

%% 6) Bottom Buttons: Plot, Generate ASCII, Save Wave Data
btnWidth  = 100;
btnHeight = 30;
spacing   = 10;
yBottom   = 40;

handles.btnPlot = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'Position',[xLabel yBottom btnWidth+20 btnHeight],...
    'String','Plot Wave','FontWeight','bold');

handles.btnSaveData = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Save Wave Data',...
    'Position',[xLabel+btnWidth+spacing+20 yBottom btnWidth btnHeight]);

handles.btnGen = uicontrol('Parent',ctrlPanel,'Style','pushbutton',...
    'String','Generate ASCII',...
    'Position',[xLabel+(btnWidth+spacing+20)*1.85 yBottom btnWidth+20 btnHeight]);

%% 7) Wave Metrics placed near the Status Label
xStatus = panelWidth+60;
yStatus = 50;
xOffset = 0; 
metricLabelW = 120;
metricEditW  = 60;
gap          = 15;

uicontrol('Parent',fig,'Style','text',...
    'String','Measured Max (m):',...
    'HorizontalAlignment','left',...
    'Position',[xStatus+xOffset, yStatus, metricLabelW, 20]);
handles.edtMeasuredMaxH = uicontrol('Parent',fig,'Style','edit',...
    'String','-', 'Enable','inactive',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Position',[xStatus+xOffset+metricLabelW, yStatus, metricEditW, 20]);

uicontrol('Parent',fig,'Style','text',...
    'String','Estimated Max (m):',...
    'HorizontalAlignment','left',...
    'Position',[xStatus+xOffset+metricLabelW+metricEditW+gap, yStatus, metricLabelW, 20]);
handles.edtEstMaxH = uicontrol('Parent',fig,'Style','edit',...
    'String','-', 'Enable','inactive',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Position',[xStatus+xOffset+2*metricLabelW+metricEditW+gap, yStatus, metricEditW, 20]);

uicontrol('Parent',fig,'Style','text',...
    'String','T_z (s):',...
    'HorizontalAlignment','left',...
    'Position',[xStatus+xOffset+2*metricLabelW+2*metricEditW+2*gap, yStatus, 50, 20]);
handles.edtTzDisplay = uicontrol('Parent',fig,'Style','edit',...
    'String','-', 'Enable','inactive',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Position',[xStatus+xOffset+2*metricLabelW+2*metricEditW+2*gap+50, yStatus, metricEditW, 20]);

% Additional references:
handles.modelList = {'JONSWAP 2-param','JONSWAP 3-param','JONSWAP 6-param',...
    'Torsethaugen (testing)','Ochi-Hubble (testing)',...
    'PM (1-param)','PM (2-param, Tavg)','PM (2-param, Tzc)'};

handles.lastWave = [];
handles.lastSpec = [];

end % buildGUI


%% =========================================================================
%  SUBFUNCTION: SET/INITIALIZE CALLBACKS
%% =========================================================================
function initializeCallbacks(handles, fig)
% Attach callbacks to the UI controls, forwarding "fig" so each callback
% can retrieve the updated handles via get(fig,'UserData').

set(handles.btnLoadCfg,    'Callback', @(src,evt) onLoadConfig(fig));
set(handles.btnSaveCfg,    'Callback', @(src,evt) onSaveConfig(fig));
set(handles.waveTypePopup, 'Callback', @(src,evt) onWaveTypeChanged(fig));
set(handles.spectrumPopup, 'Callback', @(src,evt) onSpectrumModelChanged(fig));
set(handles.btnPlot,       'Callback', @(src,evt) onPlotWave(fig));
set(handles.btnSaveData,   'Callback', @(src,evt) onSaveWaveData(fig));
set(handles.btnGen,        'Callback', @(src,evt) onGenerateASCII(fig));

end % initializeCallbacks

%% =========================================================================
%  CALLBACKS
%% =========================================================================
function onLoadConfig(fig)
[fileName, pathName] = uigetfile('*.mat','Load Configuration');
if fileName == 0, return; end
S = load(fullfile(pathName,fileName),'cfg');
if isfield(S,'cfg')
    setParams(fig, S.cfg);
    updateUIEnable(fig);
    clearWaveCache(fig);
    setStatus(fig,'Configuration loaded. Old wave data cleared.');
end
end

function onSaveConfig(fig)
[fileName, pathName] = uiputfile('*.mat','Save Configuration');
if fileName == 0, return; end
cfg = getAllParams(fig);
save(fullfile(pathName,fileName),'cfg');
setStatus(fig,'Configuration saved successfully.');
end

function onWaveTypeChanged(fig)
updateUIEnable(fig);
clearWaveCache(fig);
setStatus(fig,'Wave type changed. Old wave data cleared.');
end

function onSpectrumModelChanged(fig)
updateUIEnable(fig);
clearWaveCache(fig);
setStatus(fig,'Spectrum model changed. Old wave data cleared.');
end

function onPlotWave(fig)
[t, eta] = getOrGenerateWaveData(fig);
[wPlot, SPlot] = getOrComputeSpectrum(fig);

% Plot time-domain wave
handles = get(fig,'UserData');
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

% Update wave metrics
updateWaveMetrics(fig, t, eta, wPlot, SPlot);

setStatus(fig,'Plot complete.');
end

function onSaveWaveData(fig)
[t, eta] = getOrGenerateWaveData(fig);
cfg = getAllParams(fig);

[fileName, pathName] = uiputfile('*.mat','Save Wave Data');
if fileName == 0, return; end

waveData.time   = t;
waveData.eta    = eta;
waveData.config = cfg;
save(fullfile(pathName, fileName),'waveData');
setStatus(fig,['Saved wave data + config to ', fullfile(pathName, fileName)]);
end

function onGenerateASCII(fig)
[t, eta] = getOrGenerateWaveData(fig);
[~, SPlot] = getOrComputeSpectrum(fig); %#ok<ASGLU>

[fileName, pathName] = uiputfile('waveElevation_sima.dat','Generate SIMA ASCII');
if fileName == 0, return; end
outFileName = fullfile(pathName,fileName);

dt = t(2) - t(1);
cfg = getAllParams(fig);
writeSIMAAscii(outFileName, eta, dt, cfg);
setStatus(fig,sprintf('Wrote wave time series to %s', outFileName));
end

%% =========================================================================
%  WAVE DATA / SPECTRUM CACHING & UTILITIES
%% =========================================================================
function clearWaveCache(fig)
handlesC = get(fig,'UserData');
handlesC.lastWave = [];
handlesC.lastSpec = [];
set(fig,'UserData',handlesC);
end

function [t, eta] = getOrGenerateWaveData(fig)
handlesLocal = get(fig,'UserData');
cfg = getAllParams(fig);

if ~isempty(handlesLocal.lastWave) && cfgStructsEqual(handlesLocal.lastWave.cfg, cfg)
    t   = handlesLocal.lastWave.t;
    eta = handlesLocal.lastWave.eta;
    setStatus(fig,'Reusing cached wave data.');
else
    setStatus(fig,'Generating wave data...');
    [t, eta] = generateWaveFromGUI(cfg, handlesLocal.modelList);
    handlesLocal.lastWave.t   = t;
    handlesLocal.lastWave.eta = eta;
    handlesLocal.lastWave.cfg = cfg;
    set(fig,'UserData',handlesLocal);
    setStatus(fig,'Wave generation complete.');
end
end

function [wPlot, SPlot] = getOrComputeSpectrum(fig)
handlesLocal = get(fig,'UserData');
cfg = getAllParams(fig);

if ~isempty(handlesLocal.lastSpec) && cfgStructsEqual(handlesLocal.lastSpec.cfg, cfg)
    wPlot = handlesLocal.lastSpec.w;
    SPlot = handlesLocal.lastSpec.S;
    setStatus(fig,'Reusing cached spectrum.');
else
    setStatus(fig,'Computing wave spectrum...');
    [wPlot, SPlot] = computeSpectrumFromCfg(cfg, handlesLocal.modelList);
    handlesLocal.lastSpec.w   = wPlot;
    handlesLocal.lastSpec.S   = SPlot;
    handlesLocal.lastSpec.cfg = cfg;
    set(fig,'UserData',handlesLocal);
    setStatus(fig,'Spectrum computation complete.');
end
end

function tf = cfgStructsEqual(cfg1, cfg2)
% Checks if the two config structs have the same fields/values
fA = sort(fieldnames(cfg1));
fB = sort(fieldnames(cfg2));
if ~isequal(fA,fB)
    tf = false;
    return;
end
for fn = 1:numel(fA)
    if ~isequal(cfg1.(fA{fn}), cfg2.(fA{fn}))
        tf = false;
        return;
    end
end
tf = true;
end

function setStatus(fig, msg)
handlesC = get(fig,'UserData');
set(handlesC.statusLabel, 'String', ['Status: ', msg]);
drawnow;
end

%% =========================================================================
%  ENABLE/DISABLE CONTROLS
%% =========================================================================
function updateUIEnable(fig)
handlesLocal = get(fig,'UserData');
waveTypeVal  = get(handlesLocal.waveTypePopup,'Value');
isRegular    = (waveTypeVal == 1);

% Regular wave fields
set([handlesLocal.edtH, handlesLocal.edtT, handlesLocal.edtBeta], ...
    'Enable', bool2OnOff(isRegular));

% Irregular wave fields
set(handlesLocal.spectrumPopup, 'Enable', bool2OnOff(~isRegular));

if ~isRegular
    modelVal = get(handlesLocal.spectrumPopup,'Value');
    modelStr = handlesLocal.modelList{modelVal};
    
    switch modelStr
        case 'JONSWAP 2-param'
            relevantFields = {'edtSIWAHE','edtTpeak'};
        case 'JONSWAP 3-param'
            relevantFields = {'edtSIWAHE','edtTpeak','edtGamma'};
        case 'JONSWAP 6-param'
            relevantFields = {'edtOmegaP','edtAlpha','edtBetaVal','edtGamma','edtSigA','edtSigB'};
        case 'Torsethaugen (testing)'
            relevantFields = {'edtSIWAHE','edtTpeak'};
        case 'Ochi-Hubble (testing)'
            relevantFields = {'edtSIWAHE'};
        case 'PM (1-param)'
            relevantFields = {'edtSIWAHE'};
        case 'PM (2-param, Tavg)'
            relevantFields = {'edtSIWAHE','edtTavg'};
        case 'PM (2-param, Tzc)'
            relevantFields = {'edtSIWAHE','edtTzc'};
        otherwise
            relevantFields = {};
    end
    
    allIrrHandles = {'edtSIWAHE','edtTpeak','edtGamma','edtOmegaP','edtAlpha',...
        'edtBetaVal','edtSigA','edtSigB','edtTavg','edtTzc'};
    % Disable all first
    for ii=1:numel(allIrrHandles)
        set(handlesLocal.(allIrrHandles{ii}), 'Enable','off');
    end
    % Then enable relevant ones
    for ii=1:numel(relevantFields)
        set(handlesLocal.(relevantFields{ii}), 'Enable','on');
    end
else
    % Regular wave => disable all irregular fields
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

%% =========================================================================
%  GET/SET PARAMS
%% =========================================================================
function cfg = getAllParams(fig)
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

function setParams(fig, cfg)
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

set(fig,'UserData',handlesLocal);  % Save changes
end

%% =========================================================================
%  WAVE GENERATION & SPECTRUM
%% =========================================================================
function [t, eta] = generateWaveFromGUI(cfg, modelList)
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
    modelStr = modelList{cfg.spectrumVal};
    specFun = selectSpectrumFunction(modelStr, cfg);
    eta = generateIrregularWave(t, cfg.tRamp, specFun);
end
end

function [wPlot, SPlot] = computeSpectrumFromCfg(cfg, modelList)
if cfg.waveTypeVal == 1
    % Regular wave -> approximate delta spike near w0
    w0 = 2*pi/cfg.T;
    wMax = max(2.5*w0, 2.5);
    dw = 0.01;
    wPlot = (0:dw:wMax).';
    SPlot = zeros(size(wPlot));
    idxPeak = round(w0/dw)+1;
    if idxPeak>0 && idxPeak<=numel(SPlot)
        SPlot(idxPeak) = (cfg.H^2/8)/dw; % approximate delta
    end
else
    modelStr = modelList{cfg.spectrumVal};
    wMax = estimateWmax(cfg, modelStr);
    dw   = 0.01;
    wPlot = (0:dw:wMax).';
    
    specFun = selectSpectrumFunction(modelStr, cfg);
    SPlot   = specFun(wPlot);
    SPlot(wPlot<=0) = 0;
end
end

function wMax = estimateWmax(cfg, modelStr)
% A basic guess for wMax based on the input parameters
if cfg.Tpeak > 0
    wMax = 5*(2*pi/cfg.Tpeak);
else
    wMax = 5*(2*pi/8); % fallback
end

if strcmp(modelStr,'JONSWAP 6-param') && (cfg.omegaP>0)
    wMax = max(wMax, 4*cfg.omegaP);
elseif strcmp(modelStr,'PM (2-param, Tavg)') && (cfg.Tavg>0)
    wMax = max(wMax, 5*(2*pi/cfg.Tavg));
elseif strcmp(modelStr,'PM (2-param, Tzc)') && (cfg.Tzc>0)
    wMax = max(wMax, 5*(2*pi/cfg.Tzc));
end
if wMax<2.5, wMax=2.5; end
end

function specFun = selectSpectrumFunction(modelStr, cfg)
switch modelStr
    case 'JONSWAP 2-param'
        specFun = @(w) jonswap2_skipZero(w, cfg.Hs, cfg.Tpeak);
    case 'JONSWAP 3-param'
        specFun = @(w) jonswap3_skipZero(w, cfg.Hs, cfg.Tpeak, cfg.gamma);
    case 'JONSWAP 6-param'
        specFun = @(w) jonswap6(w, cfg.omegaP, cfg.alpha, cfg.betaVal, cfg.gamma, cfg.sigA, cfg.sigB);
    case 'Torsethaugen (testing)'
        specFun = @(w) torsethaugenSpectrum_skipZero(w, cfg.Hs, cfg.Tpeak);
    case 'Ochi-Hubble (testing)'
        specFun = @(w) ochiHubbleSpectrum_skipZero(w, cfg.Hs);
    case 'PM (1-param)'
        specFun = @(w) pmSpectrum1(w, cfg.Hs);
    case 'PM (2-param, Tavg)'
        specFun = @(w) pmSpectrum2_AVWAPE(w, cfg.Hs, cfg.Tavg);
    case 'PM (2-param, Tzc)'
        specFun = @(w) pmSpectrum2_ZCWAPE(w, cfg.Hs, cfg.Tzc);
    otherwise
        % fallback: zero
        specFun = @(w) zeros(size(w));
end
end

%% =========================================================================
%  WAVE METRICS
%% =========================================================================
function updateWaveMetrics(fig, t, eta, w, S)
hMaxMeasured = max(eta) - min(eta);

handlesLocal = get(fig,'UserData');
cfg = handlesLocal.lastWave.cfg; 

if cfg.waveTypeVal == 1
    % Regular wave => Tz, estimated Hmax not especially relevant
    TzStr   = 'N/A';
    HmaxEst = 'N/A';
else
    % Irregular wave => Tz = 2*pi * sqrt(m0 / m2)
    m0 = trapz(w, S);
    m2 = trapz(w, w.^2 .* S);
    if m2 > 0
        Tz_val = 2*pi * sqrt(m0 / m2);
    else
        Tz_val = 0;
    end
    
    Tdur = t(end) - t(1);
    Nwaves = max(Tdur / max(Tz_val,1e-8), 1);
    HmaxEst_val = cfg.Hs * sqrt(0.5 * log(Nwaves));
    
    TzStr   = sprintf('%.3f', Tz_val);
    HmaxEst = sprintf('%.3f', HmaxEst_val);
end

set(handlesLocal.edtMeasuredMaxH,'String', sprintf('%.3f', hMaxMeasured));
set(handlesLocal.edtTzDisplay,   'String', TzStr);
set(handlesLocal.edtEstMaxH,     'String', HmaxEst);
end

%% =========================================================================
%  UTILITY SUBFUNCTIONS
%% =========================================================================
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
S(w<=0) = 0;

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

%% ----------------- SPECTRUM MODELS ----------------------------------
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

function S = torsethaugenSpectrum_skipZero(w, Hs, Tp)
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);
f_ = w_/(2*pi);

S_hz = torsethaugenSpectrum(f_, Hs, Tp);
S_rad= S_hz .* (2*pi);
S(idx) = S_rad;
end

function S = torsethaugenSpectrum(f, Hs, Tp)
f = f(:);
[S1, S2] = torsetComponents(f, Hs, Tp);
S = S1 + S2;
end

function [S1, S2] = torsetComponents(f, Hs, Tp)
Hs1 = 0.6 * Hs;
Tp1 = Tp;  
gamma1 = 3.3;

Hs2 = sqrt(Hs^2 - Hs1^2);
Tp2 = Tp + 4;
gamma2 = 1.0;

S1 = jonswapPartial(f, Hs1, Tp1, gamma1);
S2 = jonswapPartial(f, Hs2, Tp2, gamma2);
end

function Sj = jonswapPartial(f, Hs_j, Tp_j, gamma_j)
g = 9.81;
f = f(:);
wp_j = 2*pi/Tp_j;
sigma = 0.07.*(2*pi*f <= wp_j) + 0.09.*(2*pi*f > wp_j);

S_base = 0.0081*g^2 ./ ((2*pi*f).^5) .* exp(-5/4 * (wp_j ./ (2*pi*f)).^4);
peakEnh= gamma_j.^(exp(- ((2*pi*f - wp_j).^2) ./ (2*sigma.^2 .* wp_j.^2)));
Sj_ = S_base .* peakEnh;

Sj = scaleSpectrumToHs(Sj_, 2*pi*f, Hs_j);
end

function S = ochiHubbleSpectrum_skipZero(w, Hs)
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);
f_ = w_/(2*pi);

Hs1 = 0.6*Hs;  Tp1 = 12; lam1 = 2;
Hs2 = sqrt(Hs^2 - Hs1^2);  Tp2 = 8;  lam2 = 3.3;

S_hz = ochiHubbleSpectrum(f_, Hs1, Tp1, lam1, Hs2, Tp2, lam2);
S_rad= S_hz.*(2*pi);
S(idx) = S_rad;
end

function S = ochiHubbleSpectrum(f, Hs1, Tp1, lam1, Hs2, Tp2, lam2)
f = f(:);
S1 = singleOchiHubble(f, Hs1, Tp1, lam1);
S2 = singleOchiHubble(f, Hs2, Tp2, lam2);
S  = S1 + S2;
S(f<=0) = 0;
end

function S1 = singleOchiHubble(f, Hs, Tp, lam)
fp  = 1/Tp;
Gamma = gamma(lam);
A  = ((4*lam + 1)/4)^lam * fp^(4*lam) * (Hs^2)/(4*Gamma);
S1 = A.* f.^(-(4*lam + 1)) .* exp(-((4*lam + 1)/4).*(fp./f).^4);
S1(f<=0)=0;
end

function S = pmSpectrum1(w, Hs)
g = 9.81;
S = zeros(size(w));
idx = (w>0);
w_ = w(idx);

A  = Hs/4; 
wp = 0.44*g/A;  % Simplified approach
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

function Sscaled = scaleSpectrumToHs(S_partial, w_partial, Hs)
m0 = trapz(w_partial, S_partial);
Hs_model = 4*sqrt(m0);
if Hs_model == 0
    Sscaled = S_partial;
else
    scaleFactor = (Hs / Hs_model)^2;
    Sscaled = S_partial * scaleFactor;
end
end

%% =========================================================================
%  WRITE ASCII
%% =========================================================================
function writeSIMAAscii(outFileName, eta, dt, cfg)
fid = fopen(outFileName, 'w');
if fid < 0
    error('Could not open file: %s', outFileName);
end

N = length(eta);
fprintf(fid, '%d\n', N);
fprintf(fid, '%.4f\n', dt);

fprintf(fid, 'Generated by waveGUI_v3_refactored (v4.0). Author: Shuijin Li (Email: lishuijin@nbu.edu.cn)\n');

paramLine = buildParameterLine(cfg);
fprintf(fid, '%s\n', paramLine);

valsPerLine = 1;
for i = 1:N
    fprintf(fid, '%12.6f', eta(i));
    if mod(i, valsPerLine)==0
        fprintf(fid,'\n');
    end
end
if mod(N,valsPerLine)~=0
    fprintf(fid,'\n');
end
fclose(fid);
end

function paramLine = buildParameterLine(cfg)
if cfg.waveTypeVal == 1
    paramPairs = {
        'Time Step (s)',                cfg.tStep
        'Number of Samples',            cfg.nSamples
        'Ramp Duration (sec)',          cfg.tRamp
        'Wave Height, H',               cfg.H
        'Wave Period, T',               cfg.T
        'Phase, beta',                  cfg.beta
        };
else
    switch cfg.spectrumVal
        case 1 % 'JONSWAP 2-param'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                'T_p',                 cfg.Tpeak
                };
        case 2 % 'JONSWAP 3-param'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                'T_p',                 cfg.Tpeak
                'γ',                   cfg.gamma
                };
        case 3 % 'JONSWAP 6-param'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'ω_p',                 cfg.omegaP
                'α',                   cfg.alpha
                'β',                   cfg.betaVal
                'γ',                   cfg.gamma
                'σ_a',                 cfg.sigA
                'σ_b',                 cfg.sigB
                };
        case 4 % 'Torsethaugen (testing)'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                'T_p',                 cfg.Tpeak
                };
        case 5 % 'Ochi-Hubble (testing)'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                };
        case 6 % 'PM (1-param)'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                };
        case 7 % 'PM (2-param, Tavg)'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                'T_{avg}',             cfg.Tavg
                };
        case 8 % 'PM (2-param, Tzc)'
            paramPairs = {
                'Time Step (s)',       cfg.tStep
                'Number of Samples',   cfg.nSamples
                'Ramp Duration',       cfg.tRamp
                'Random Seed',         cfg.rSeed
                'H_s',                 cfg.Hs
                'T_z',                 cfg.Tzc
                };
        otherwise
            paramPairs = {
                'Time Step (s)',     cfg.tStep
                'Number of Samples', cfg.nSamples
                'Ramp Duration',     cfg.tRamp
                'Unknown ModelVal',  cfg.spectrumVal
                };
    end
end

paramStrings = cell(1, size(paramPairs,1));
for k = 1:size(paramPairs,1)
    name_k = paramPairs{k,1};
    val_k  = paramPairs{k,2};
    if abs(val_k - round(val_k)) < 1e-12
        paramStrings{k} = sprintf('%s=%d', name_k, round(val_k));
    else
        paramStrings{k} = sprintf('%s=%.4f', name_k, val_k);
    end
end

paramLine = ['Used parameters: ' strjoin(paramStrings, ', ')];
end
