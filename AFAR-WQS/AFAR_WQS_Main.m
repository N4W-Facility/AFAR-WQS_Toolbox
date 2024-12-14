% -------------------------------------------------------------------------
% Matlab - R2023b 
% -------------------------------------------------------------------------
%                            INFORMATION
% -------------------------------------------------------------------------
% Author        : Jonathan Nogales Pimentel
%                 Carlos Andrés Rogéliz Prada
% Email         : jonathannogales02@gmail.com
% 
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or option) any 
% later version. This program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% ee the GNU General Public License for more details. You should have 
% received a copy of the GNU General Public License along with this program.  
% If not, see http://www.gnu.org/licenses/.
% -------------------------------------------------------------------------
%                            DESCRIPCIÓN
% -------------------------------------------------------------------------
% This code estimates the concentrations, loads and assimilation factors for
% the basin for the following water quality determinants:
%   T       : Water Temperature
%   DO      : Dissolved oxygen
%   OM      : Organic Matter
%   TC      : Total coliforms
%   SS      : Suspended solids
%   NO      : Organic nitrogen
%   NH4     : Ammoniacal nitrogen
%   NO3     : Nitrates
%   PO      : Organic phosphorus
%   PI      : Inorganic Phosphorus
%   Hg0     : Elemental mercury (metallic)
%   Hg2     : Divalent Mercury
%   MeHg    : Methyl Mercury

%% Preliminary
warning off
clearvars
close all 
clc

%% Path de proyecto 
ProjectPath = pwd;

%% Read Reach 
Tmp         = readmatrix(fullfile(ProjectPath,'Example_Project','INPUTS','Main_AFAR-WQS.xlsx'),'Sheet','INPUTS');

ShpPath     = fullfile('Example_Project','INPUTS','Network','Network.shp');
ReachData   = shaperead(ShpPath);
CodeNet     = [ReachData.ReachID]';
[~,Posi]    = ismember(CodeNet,Tmp(:,1));

for ii = 1:length(Posi)   
    % Starting node of the reach
    ReachData(ii).('FromNode')      = Tmp(Posi(ii),2);
    % End node of the reach
    ReachData(ii).('ToNode')        = Tmp(Posi(ii),3);
    % River type [mountain - 1, plane 0]
    ReachData(ii).('RiverType')     = Tmp(Posi(ii),4);
    % River mouth code
    ReachData(ii).('RiverMouth')    = logical(Tmp(Posi(ii),5));
    % Elevation [m.a.s.l.]
    ReachData(ii).('Z')             = Tmp(Posi(ii),6);
    % Area [m^2]
    ReachData(ii).('A')             = Tmp(Posi(ii),7);
    % Reach length [m]
    ReachData(ii).('L')             = Tmp(Posi(ii),8);
    % Discharge [m3/s]
    ReachData(ii).('Qr')            = Tmp(Posi(ii),9);    
    % Wigth [m]
    ReachData(ii).('W')             = Tmp(Posi(ii),10);
    % Depth [m]
    ReachData(ii).('H')             = Tmp(Posi(ii),11);
    % Flow velocity [m/s]
    ReachData(ii).('v')             = Tmp(Posi(ii),12);
    % Slope [m/m]
    ReachData(ii).('S')             = Tmp(Posi(ii),13);
    % Temperature
    ReachData(ii).('T')             = Tmp(Posi(ii),14);
    % Wastewater Discharge [m3/s]
    ReachData(ii).('Qwwd')          = Tmp(Posi(ii),15);
end
[~,id]      = sort([ReachData.ReachID]');
ReachData   = ReachData(id);

%% Instanciar Objeto
AFAR_Obj = AFAR_WQS(ReachData);

% ---------------------------------------------------------------------
% Modelo - Temperature (T)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_T('Load_T',Tmp(:,16)*NaN)

% ---------------------------------------------------------------------
% Modelo - Suspended Solids (SS)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_SS('Load_SS',Tmp(:,17))

% ---------------------------------------------------------------------
% Modelo - Cotal Coliform (X)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_X('Load_X',Tmp(:,18));

% ---------------------------------------------------------------------
% Modelo - Organic Nitrogen (NO)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_NO('Load_NO',Tmp(:,19));

% ---------------------------------------------------------------------
% Modelo - Ammoniacal Nitrogen (NH4)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_NH4('Load_NH4',Tmp(:,20))

% ---------------------------------------------------------------------
% Modelo - Nitrate (NO3)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_NO3('Load_NO3',Tmp(:,21))

% ---------------------------------------------------------------------
% Modelo - Organic Phosphorus (PO)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_PO('Load_PO',Tmp(:,22));

% ---------------------------------------------------------------------
% Modelo - Onorganic Phosphorus (PI)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_PI('Load_PI',Tmp(:,23));

% ---------------------------------------------------------------------
% Modelo - Organic Matter (OM)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_OM('Load_OM',Tmp(:,24));

% ---------------------------------------------------------------------
% Modelo - Dissolved Oxygen Deficit (DO)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_DOD;
    
% ---------------------------------------------------------------------
% Modelo - Elemental Mercury (Hg0)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_Hg2('Load_Hg2',Tmp(:,25));

% ---------------------------------------------------------------------
% Modelo - Divalent Mercury (Hg2)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_Hg0('Load_Hg0',Tmp(:,26));

% ---------------------------------------------------------------------
% Modelo - Methyl Mercury (MeHg)
% ---------------------------------------------------------------------
AFAR_Obj.WQS_MeHg('Load_MeHg',Tmp(:,27));

% ---------------------------------------------------------------------
% Export Data
% ---------------------------------------------------------------------
% Create folder
mkdir(fullfile(ProjectPath,'Example_Project','OUTPUTS'))

for ii = 1:length(Posi) 
    % River discharge [m3/s]
    ReachData(ii).('Q')         = AFAR_Obj.Q(Posi(ii));
    % Temperature [°C]
    ReachData(ii).('T')         = AFAR_Obj.T(Posi(ii));
    % Suspended Solids Concentration [mg/l]
    ReachData(ii).('C_SS')      = AFAR_Obj.C_SS(Posi(ii));
    % Total Coliform Concentration [CFU/l]
    ReachData(ii).('C_X')       = AFAR_Obj.C_X(Posi(ii));
    % Turbidity [NTU]
    ReachData(ii).('Turb')      = AFAR_Obj.Turb(Posi(ii));
    % Inorganic Phosphorus Concentration [mg/l]
    ReachData(ii).('C_PO')      = AFAR_Obj.C_PO(Posi(ii));
    % Organic Phosphorus Concentration [mg/l]
    ReachData(ii).('C_PI')      = AFAR_Obj.C_PI(Posi(ii));
    % Total Phosphorus Concentration [mg/l]
    ReachData(ii).('C_PT')      = AFAR_Obj.C_PO(Posi(ii)) + AFAR_Obj.C_PI(Posi(ii));
    % Organic Nitrogen Concentration [mg/l]
    ReachData(ii).('C_NO')      = AFAR_Obj.C_NO(Posi(ii));
    % Ammonia Nitrogen Concentration [mg/l]
    ReachData(ii).('C_NH4')     = AFAR_Obj.C_NH4(Posi(ii));
    % Nitrate Concentration [mg/l]
    ReachData(ii).('C_NO3')     = AFAR_Obj.C_NO3(Posi(ii));
    % Total Nitrogen Concentration [mg/l]
    ReachData(ii).('C_NT')      = AFAR_Obj.C_NO(Posi(ii)) + AFAR_Obj.C_NH4(Posi(ii)) + AFAR_Obj.C_NO3(Posi(ii));
    % Organic Matter Concentration [mg/l]
    ReachData(ii).('C_OM')      = AFAR_Obj.C_OM(Posi(ii));
    % Saturation Dissolved Oxygen Concentration [mg/l]
    ReachData(ii).('C_DO')      = AFAR_Obj.C_DO(Posi(ii));
    % Saturation dissolved oxygen concentration [mg/l]
    ReachData(ii).('C_OS')      = AFAR_Obj.C_OS(Posi(ii));
    % Elemental Mercury Concentration [ng/l]
    ReachData(ii).('C_Hg0')     = AFAR_Obj.C_Hg0(Posi(ii));
    % Divalent Mercury Concentration [ng/l]
    ReachData(ii).('C_Hg2')     = AFAR_Obj.C_Hg2(Posi(ii));
    % Methyl Mercury Concentration [ng/l]
    ReachData(ii).('C_MeHg')    = AFAR_Obj.C_MeHg(Posi(ii));
    % Total Mercury Concentration [ng/l]
    ReachData(ii).('C_HgT')     = (AFAR_Obj.C_Hg0(Posi(ii)) + AFAR_Obj.C_Hg2(Posi(ii)) + AFAR_Obj.C_MeHg(Posi(ii)));
end
shapewrite(ReachData,fullfile(ProjectPath,'Example_Project','OUTPUTS','Network_AFAR-WQS.shp'));