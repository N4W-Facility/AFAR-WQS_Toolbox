classdef AFAR_WQS < ClassNetwork & matlab.mixin.Copyable
    % ---------------------------------------------------------------------
    % Matlab - R2023b 
    % ---------------------------------------------------------------------
    %                            INFORMATION
    % ----------------------------------------------------------------------
    % Autor         : Jonathan Nogales Pimentel
    %                 Carlos Andrés Rogéliz Prada
    % Email         : jonathannogales02@gmail.com
    % Company       : The Nature Conservancy - TNC
    % Date          : October, 2024
    % 
    % ----------------------------------------------------------------------
    % This program is free software: you can redistribute it and/or modify it 
    % under the terms of the GNU General Public License as published by the 
    % Free Software Foundation, either version 3 of the License, or option) any 
    % later version. This program is distributed in the hope that it will be 
    % useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    % ee the GNU General Public License for more details. You should have 
    % received a copy of the GNU General Public License along with this program.  
    % If not, see http://www.gnu.org/licenses/.

    % ---------------------------------------------------------------------
    % Temporal variables
    % ---------------------------------------------------------------------
    properties
        % Status to accumulate area [1 true or 0 false]
        StatusAccArea(1,1) logical = false
        % Boolean variable identifying whether a basin is a headwater basin 
        % or not [true or false]
        BC_Status(:,1) logical
        % Accumulated river discharge [m3/s]. This variable is calculated 
        % internally and is used to estimate pollutant loads.
        Qaccum(:,1) double
    end
    
    % ---------------------------------------------------------------------
    % River and associated watershed characteristics  
    % ---------------------------------------------------------------------
    properties
        % Bounding Box
        BoundingBox
        % X
        X
        % AND
        Y
        % X_FN
        X_FN
        % Y_FN
        Y_FN
        % X_TN
        X_TN
        % Y_TN
        Y_TN
        % Area [m^2]
        A(:,1) double         
        % River length [m]
        L(:,1) double 
        % Wastewater discharge [m3/s]
        Qwwd(:,1) double 
        % River discharge [m3/s]
        Qr(:,1) double
        % River discharge + wastewater discharge [m3/s]
        Q(:,1) double  
        % Average river elevation [m.s.n.m]
        Z(:,1) double  
        % Slope River [m/m]
        S(:,1) double                  
        % Boolean variable that identifies whether a river is a mountain  
        % river or not [true or false]. 
        % Note 1: if false is specified, the tool will assume that the river 
        % is a plain river.  
        % Note 2: To define whether a river is a plain or a mountain river, 
        % a first criterion may be to assume that the former is limited by 
        % capacity (slope ≤0.025 m/m) and the latter is limited by supply 
        % (slope >0.025 m /m). These slope thresholds are defined by 
        % Flores et al. (2006) 
        % Note 3: A second criterion can be to use the slope threshold 
        % defined by Wohl and Merritt (2005) to define whether a river is 
        % mountainous (slope > 0.002 m/m) or plain (slope < 0.002 m/m).
        RiverType(:,1) logical
    end
    
    % -------------------------------------------------- ------------------- 
    % Hydraulic and transport characteristics of the river.  
    % Note: These data can be provided by the user, however, if they are  
    % not provided, they can be calculated by indirect methods.  
    % -------------------------------------------------- -------------------
    properties
        % River width [m]
        W(:,1) double
        % Depth River [m]
        H(:,1) double
        % River flow velocity [m/s]
        v(:,1) double
        % Water temperature [°C]
        T(:,1) double
        % Arrival or advection time of the river [day]
        tao(:,1) double
        % River travel time [day]
        tb(:,1) double
        % Residence time in dead zone of the river [day]
        Tr(:,1) double
        % Dispersive fraction [dimensionless]
        DF(:,1) double
        % Solute velocity in the river [m/s]
        Vs(:,1) double
        % Effective delay coefficient [dimensionless]
        Beta          
        % Turbidity in nephelometric turbidity units [NTU]
        Turb
    end
    
    % --------------------------------------------------------------------- 
    % Model Execution Status
    % ---------------------- ----------------------------------------------
    properties         
        % Total Coliform Concentration [CFU/l]
        StatusWQS_X logical = false;
        % Suspended Solids Concentration [mg/l]
        StatusWQS_SS logical = false;
        % Inorganic Phosphorus Concentration [mg/l]
        StatusWQS_PO logical = false;
        % Organic Phosphorus Concentration [mg/l]
        StatusWQS_PI logical = false;
        % Organic Nitrogen Concentration [mg/l]
        StatusWQS_NO logical = false; 
        % Ammonia Nitrogen Concentration [mg/l]
        StatusWQS_NH4 logical = false;
        % Nitrate Concentration [mg/l]
        StatusWQS_NO3 logical = false;
        % Organic Matter Concentration [mg/l]
        StatusWQS_OM logical = false;
        % Dissolved Oxygen Deficit [mg/l]
        StatusWQS_DOD logical = false;
        % Saturation Dissolved Oxygen Concentration [mg/l]
        StatusWQS_DO logical = false;
        % Saturation dissolved oxygen concentration [mg/l]
        StatusWQS_OS logical = false;
        % Temperature [°C]
        StatusWQS_T logical = false;
        % Elemental Mercury Concentration [mg/l]
        StatusWQS_Hg0 logical = false;
        % Divalent Mercury Concentration [mg/l]
        StatusWQS_Hg2 logical = false;
        % Methyl Mercury Concentration [mg/l]
        StatusWQS_MeHg logical = false;
    end

    % --------------------------------------------------------------------- 
    % Variables that store the concentration of determinants in river 
    % ---------------------- ----------------------------------------------
    properties         
        % Total Coliform Concentration [CFU/l]
        C_X(:,1) double
        % Suspended Solids Concentration [mg/l]
        C_SS(:,1) double
        % Turbidity [NTU]
        C_Turb(:,1) double
        % Inorganic Phosphorus Concentration [mg/l]
        C_PO(:,1) double
        % Organic Phosphorus Concentration [mg/l]
        C_PI(:,1) double
        % Organic Nitrogen Concentration [mg/l]
        C_NO(:,1) double    
        % Ammonia Nitrogen Concentration [mg/l]
        C_NH4(:,1) double
        % Nitrate Concentration [mg/l]
        C_NO3(:,1) double
        % Organic Matter Concentration [mg/l]
        C_OM(:,1) double
        % Dissolved Oxygen Deficit [mg/l]
        C_DOD(:,1) double
        % Saturation Dissolved Oxygen Concentration [mg/l]
        C_DO(:,1) double
        % Saturation dissolved oxygen concentration [mg/l]
        C_OS(:,1) double
        % Temperature [°C]
        C_T(:,1) double
        % Elemental Mercury Concentration [mg/l]
        C_Hg0(:,1) double
        % Divalent Mercury Concentration [mg/l]
        C_Hg2(:,1) double
        % Methyl Mercury Concentration [mg/l]
        C_MeHg(:,1) double
    end
    
    % --------------------------------------------------------------------- 
    % Variables that store the load of determinants in the river 
    % ---------------------------------------------------------------------
    properties         
        % Downstream load of Total Coliform[CFU]
        W_X(:,1) double
        % Downstream load of Total Suspended Solids [mg]
        W_SS(:,1) double
        % Downstream load of Organic Phosphorus [mg]
        W_PO(:,1) double
        % Downstream load of Inorganic Phosphorus [mg]
        W_PI(:,1) double
        % Downstream load of Organic Nitrogen [mg]
        W_NO(:,1) double    
        % Downstream load of Ammoniacal Nitrogen [mg]
        W_NH4(:,1) double
        % Downstream load of Nitrate [mg]
        W_NO3(:,1) double
        % Downstream load of Organic Matter Loading [mg]
        W_OM(:,1) double
        % Downstream load of Dissolved Oxygen Deficit [mg]
        W_DOD(:,1) double
        % Downstream load of Dissolved Oxygen Saturation Oxygen [mg]
        W_DO(:,1) double
        % Downstream load of Dissolved Oxygen Dissolved Oxygen [mg]
        W_OS(:,1) double
        % Downstream load of Temperature [°C]
        W_T(:,1) double
        % Downstream load of Elemental Mercury [mg]
        W_Hg0(:,1) double
        % Downstream load of Divalent Mercury [mg]
        W_Hg2(:,1) double
        % Downstream load of Load of Methyl Mercury [mg]
        W_MeHg(:,1) double
    end
    
    % --------------------------------------------------------------------- 
    % Variables that store the concentration of determinants in a 
    % Waste Water Disposal 
    % ------------------- --------------------------------------------------
    properties         
        % Upstream load of Total Coliform[CFU]
        Load_X(:,1) double
        % Upstream load of Total Suspended Solids [mg]
        Load_SS(:,1) double
        % Upstream load of Organic Phosphorus [mg]
        Load_PO(:,1) double
        % Upstream load of Inorganic Phosphorus [mg]
        Load_PI(:,1) double
        % Upstream load of Organic Nitrogen [mg]
        Load_NO(:,1) double    
        % Upstream load of Ammoniacal Nitrogen [mg]
        Load_NH4(:,1) double
        % Upstream load of Nitrate [mg]
        Load_NO3(:,1) double
        % Upstream load of Organic Matter load [mg]
        Load_OM(:,1) double
        % Upstream load of Dissolved Oxygen Deficit [mg]
        Load_DOD(:,1) double
        % Upstream load of Dissolved Oxygen [mg]
        Load_DO(:,1) double
        % Upstream load of Saturation Dissolved Oxygen Saturation [mg]
        Load_OS(:,1) double
        % Upstream load of Temperature [°C]
        Load_T(:,1) double
        % Upstream load of Elemental Mercury [mg]
        Load_Hg0(:,1) double
        % Upstream load of Divalent Mercury [mg]
        Load_Hg2(:,1) double
        % Upstream load of Load of Methyl Mercury [mg]
        Load_MeHg(:,1) double
        % if the variable is “true” the determinant dump has not entered 
        % the reach. If it is “False”, the determinant dump has already 
        % entered.
        Load_Status(:,1) logical
    end
    
    % ---------------------------------------------------------------------
    % Load of determinants on boundary conditions 
    % ---------------------------------------------------------------------
    properties    
        % Total Coliform boundary condition [UFC]
        BC_X(:,1) double
        % Suspended Solids boundary condition [mg]
        BC_SS(:,1) double
        % Organic Phosphorus boundary condition [mg]
        BC_PO(:,1) double
        % Inorganic Phosphorus boundary condition [mg]
        BC_PI(:,1) double
        % Organic Nitrogen boundary condition [mg]
        BC_NO(:,1) double    
        % Ammoniacal Nitrogen boundary condition [mg]
        BC_NH4(:,1) double
        % Nitrate boundary condition [mg]
        BC_NO3(:,1) double
        % Organic Matter boundary condition [mg]
        BC_OM(:,1) double
        % Dissolved Oxygen Deficit boundary condition [mg]
        BC_DOD(:,1) double
        % Saturation Dissolved Oxygen boundary condition [mg]
        BC_DO(:,1) double
        % Dissolved Oxygen boundary condition[mg]
        BC_OS(:,1) double
        % Temperature boundary condition [°C]
        BC_T(:,1) double
        % Elemental Mercury boundary condition [mg]
        BC_Hg0(:,1) double
        % Divalent Mercury boundary condition [mg]
        BC_Hg2(:,1) double
        % Methyl Mercury boundary condition [mg]
        BC_MeHg(:,1) double
    end
    
    % --------------------------------------------------------------------- 
    % Determinant assimilation factors 
    % ---------------------------------------------------------------------
    properties         
        % Total Coliform Assimilation Factor (Liter)
        AF_X(:,1) double
        % Suspended Solids Assimilation Factor (Liter)
        AF_SS(:,1) double
        % Organic Phosphorus Assimilation Factor (Liter)
        AF_PO(:,1) double
        % Inorganic Phosphorus Assimilation Factor (Liter)
        AF_PI(:,1) double
        % Organic Nitrogen Assimilation Factor (Liter)
        AF_NO(:,1) double  
        % Ammoniacal Nitrogen Assimilation Factor (Liter)
        AF_NH4(:,1) double        
        % Nitrates Assimilation Factor (Liter)
        AF_NO3(:,1) double
        % Organic Matter Assimilation Factor (Liter)
        AF_OM(:,1) double
        % Dissolve Oxygen Deficit Assimilation Factor (Liter)
        AF_DOD(:,1) double
        % Temperature Assimilation Factor (Liter)
        AF_T(:,1) double
        % Elemental Mercury Assimilation Factor (Liter)
        AF_Hg0(:,1) double
        % Divalent Mercury Assimilation Factor (Liter)
        AF_Hg2(:,1) double
        % Methyl Mercury Assimilation Factor (Liter)
        AF_MeHg(:,1) double
        % Discharge conversion factor - [m^3/s] - [lts/day]
        ConFac_Q = 1000*3600*24;
        % Temporal Discharge
        QtmpT
    end    
    
    % --------------------------------------------------------------------- 
    % Initialization of Water Quality object 
    % ---------------------------------------------------------------------
    methods
        function obj = AFAR_WQS(ReachData,varargin)
                        
            % BoundingBox
            obj.BoundingBox = {ReachData.BoundingBox}';
            % X
            obj.X           = {ReachData.X}';
            % Y
            obj.Y           = {ReachData.Y}';
            % Assign mandatory input data reach code
            obj.ReachID     = [ReachData.ReachID]';
            % Starting node of the reach
            obj.FromNode    = [ReachData.FromNode]';
            % End node of the reach
            obj.ToNode      = [ReachData.ToNode]';
            % River mouth code
            obj.RiverMouth  = obj.ReachID([ReachData.RiverMouth]');
            % River type [mountain - 1, plane 0]
            obj.RiverType   = [ReachData.RiverType]';
            % Area [m^2]
            obj.A           = [ReachData.A]';
            % Reach length [m]
            obj.L           = [ReachData.L]';
            % River Discharge [m3/s]
            obj.Qr          = [ReachData.Qr]';
            % Wastewater Discharge [m3/s]
            obj.Qwwd        = [ReachData.Qwwd]';            
            % Total Discharge [m3/s]
            % Elevation [m.a.s.l.]
            obj.Z           = [ReachData.Z]';
            % Slope [m/m]
            obj.S           = [ReachData.S]';
             % Wigth [m]
            obj.W           = [ReachData.W]';
            % Depth [m]
            obj.H           = [ReachData.H]';
            % Flow velocity [m/s]
            obj.v           = [ReachData.v]';
            % Water temperature [°C]
            obj.T           = [ReachData.T]';

            % -------------------------------------------------------------
            % Check for optional input data
            % -------------------------------------------------------------
            ip = inputParser;    
            % Check - Solute velocity
            addParameter(ip, 'Vs',[],@ismatrix)
            % Check - Dispersive fraction
            addParameter(ip, 'DF',[],@ismatrix)
            % Check - Travel time
            addParameter(ip, 'tb',[],@ismatrix)
            % Check - Arrival or advection time
            addParameter(ip, 'tao',[],@ismatrix)   
            
            % Optional data check
            parse(ip,varargin{:})
                     
            % Dispersive fraction [Ad]
            DF  = ip.Results.DF;
            % Solute velocity [m/s]
            Vs  = ip.Results.Vs;
            % Travel time [d]
            tb  = ip.Results.tb;
            % Arrival or advection time [d]
            tao = ip.Results.tao;
         
            % Assign discharge status
            obj.Update_VerStatus;            
            
            % Initialization of loads of the determinants 
            % Total Coliform
            obj.W_X         = zeros(size(obj.ReachID));
            % Suspended Solids
            obj.W_SS        = zeros(size(obj.ReachID));
            % Organic Phosphorus
            obj.W_PO        = zeros(size(obj.ReachID));
            % Inorganic Phosphorus
            obj.W_PI        = zeros(size(obj.ReachID));
            % Organic Nitrogen
            obj.W_NO        = zeros(size(obj.ReachID));
            % Ammoniacal Nitrogen
            obj.W_NH4       = zeros(size(obj.ReachID));
            % Nitrate
            obj.W_NO3       = zeros(size(obj.ReachID));
            % Organic Matter
            obj.W_OM        = zeros(size(obj.ReachID));
            % Dissolved Oxygen Deficit
            obj.W_DOD       = zeros(size(obj.ReachID));
            % Temperature
            obj.W_T         = zeros(size(obj.ReachID));
            % Elemental Mercury
            obj.W_Hg0       = zeros(size(obj.ReachID));
            % Divalent Mercury
            obj.W_Hg2       = zeros(size(obj.ReachID));
            % Methyl Mercury
            obj.W_MeHg      = zeros(size(obj.ReachID));
            
            % Initialization of Concentrations of the determinants 
            % Total Coliform
            obj.C_X         = zeros(size(obj.ReachID));
            % Suspended Solids
            obj.C_SS        = zeros(size(obj.ReachID));
            % Turbidity
            obj.Turb        = zeros(size(obj.ReachID));
            % Organic Phosphorus
            obj.C_PO        = zeros(size(obj.ReachID));
            % Inorganic Phosphorus
            obj.C_PI        = zeros(size(obj.ReachID));
            % Organic Nitrogen
            obj.C_NO        = zeros(size(obj.ReachID));
            % Ammoniacal Nitrogen
            obj.C_NH4       = zeros(size(obj.ReachID));
            % Nitrate
            obj.C_NO3       = zeros(size(obj.ReachID));
            % Organic matter
            obj.C_OM        = zeros(size(obj.ReachID));
            % Dissolved Oxygen Deficit
            obj.C_DOD       = zeros(size(obj.ReachID));
            % Dissolved Oxygen
            obj.C_DO        = zeros(size(obj.ReachID));
            % Temperature
            obj.C_T         = zeros(size(obj.ReachID));
            % Elemental Mercury
            obj.C_Hg0       = zeros(size(obj.ReachID));
            % Divalent Mercury
            obj.C_Hg2       = zeros(size(obj.ReachID));
            % Methyl mercury
            obj.C_MeHg      = zeros(size(obj.ReachID));
                        
            % Initialization of assimilation factors of the determinants 
            % Total Coliform Factor in the boundary conditions [Ad]
            obj.AF_X        = zeros(size(obj.ReachID));
            % Total Suspended Solids Factor in Boundary Conditions [Ad]
            obj.AF_SS       = zeros(size(obj.ReachID));
            % Inorganic Phosphorus Factor in Boundary Conditions [Ad]
            obj.AF_PO       = zeros(size(obj.ReachID));
            % Organic Phosphorus Factor in Boundary Conditions [Ad]
            obj.AF_PI       = zeros(size(obj.ReachID));
            % Organic Nitrogen Factor in Boundary Conditions [Ad]
            obj.AF_NO       = zeros(size(obj.ReachID));
            % Ammoniacal Nitrogen Factor in Boundary Conditions [Ad]
            obj.AF_NH4      = zeros(size(obj.ReachID));
            % Nitrate factor in boundary conditions [Ad]
            obj.AF_NO3      = zeros(size(obj.ReachID));
            % Organic Matter Factor in Boundary Conditions [Ad]
            obj.AF_OM       = zeros(size(obj.ReachID));
            % Dissolved Oxygen Deficit Factor at Boundary Conditions [Ad]
            obj.AF_DOD      = zeros(size(obj.ReachID));
            % Temperature Factor in Boundary Conditions [Ad]
            obj.AF_T        = zeros(size(obj.ReachID));
            % Elemental Mercury Factor in Boundary Conditions [Ad]
            obj.AF_Hg0      = zeros(size(obj.ReachID));
            % Divalent Mercury Factor in Boundary Conditions [Ad]
            obj.AF_Hg2      = zeros(size(obj.ReachID));
            % Methyl Mercury factor in boundary conditions [Ad]
            obj.AF_MeHg     = zeros(size(obj.ReachID));
            
            % ------------------------------------------------------------- 
            % Parameter calculation 
            % ------------------------------------------------------------- 
            % Accumulate area and discharge
            obj.VariableAccumulation;                              
            
            % Oxygen Saturation Calculation [mg/l]
            obj.Cal_OS           
                                    
            % Calculation of Solute Velocity [m/s]
            if isempty(DF)
                obj.Beta_Gonzalez;
                obj.Vs_Lees;
            else
                obj.Vs = Vs;
            end         
            
            % Calculation of the Dispersive Fraction (Ad)
            if isempty(DF) 
                obj.DF_Gonzalez;
            else
                obj.DF = DF;
            end
            
            % Calculation of transport parameters
            if isempty(tao) && ~isempty(tb) 
                % Travel time of the current [d]
                obj.tb  = tb;
                % Arrival or advection time of the current [d]
                obj.tao = obj.tb.*(1 - obj.DF);
            elseif ~isempty(tao) && isempty(tb) 
                % Arrival or advection time of the current [d]
                obj.tao = tao;
                % Travel time of the current [d]
                obj.tb = obj.tao./(1 - obj.DF);
            else
                % Conversion Factor: Sec to day
                ConFac_1 = 1./(3600*24);
                % Travel time of the current [d]
                obj.tb = (obj.L./obj.Vs).*ConFac_1;  
                % Arrival or advection time of the current [d]
                obj.tao = obj.tb.*(1 - obj.DF);
            end                                           
            
            % Calculation of dead zone residence time [d]
            obj.Tr = obj.DF.*obj.tb;                        
        end
    end         
    
    % --------------------------------------------------------------------- 
    % Assimilation factor for determinants that do not depend on the 
    % concentration of other determinants 
    % ---------------------------------------------------------------------
    methods
        function AF = Cal_AF_CoIndependent(~,Q,tao,Tr,K)
             % tao  [d]     : Arrival or advection time 
             % Tr   [d]     : Dead zone residence time 
             % K    [/d]    : Decay constant of the determinant 
             % Q    [l/d]   : Flow rate of the river section 
             % AF   [l]     : Assimilation factor of the determinant
            
            Part_1  = (1 + (K.*Tr));
            Part_2  = exp(-K.*tao);
            AF  = Q.*(Part_1./Part_2);
        end
    end
    
    % --------------------------------------------------------------------- 
    % Assimilation factor for determinants that depend on the concentration 
    % of other determinants 
    % ---------------------------------------------------------------------
    methods
        function AF = Cal_AF_CoDependent(~,Q,tao,Tr,Kd,K,Factor,Factor2)
             % tao  [d]     : Arrival or advection time 
             % Tr   [d]     : Dead zone residence time 
             % Kd   [/d]    : Decay constant of the dependent determinant 
             % K    [/d]    : Decay constant of the determinant 
             % Q    [l/d]   : River section flow 
             % AF   [l]     : Assimilation factor of the determinant
            
            Part_1  = ((Kd.*Tr) + 1);
            Part_2  = exp(K.*tao) + (Tr.*(Factor2.*(K*Factor) + Kd));
            AF      = Q.*(Part_1./Part_2);
         end
    end
    
    % --------------------------------------------------------------------- 
    % Cumulative function for determinants that do not depend on the 
    % concentration of other determinants 
    % ---------------------------------------------------------------------
    methods
        function FunctionNetwork_CoIndependent(obj, Npre, Posi, varargin)
            % This function performs the cumulative calculation of the 
            % concentrations, loads and assimilation factors of the water 
            % quality determinants that are independent, in the cumulative 
            % scheme of the Functional Branch. Model [String]: Acronym of 
            % the model to be used. By default, it is assigned a
            
            % Check-in
            ip = inputParser;
            % Acronym of the model to be used. By default, total coliforms 
            % are assigned.
            addParameter(ip,'Model','X',@ischar)
            parse(ip,varargin{:})
            Model = ip.Results.Model;
            
            % Calculation of determinant loads
            if isempty(Posi)        
                % Calculation of loads for sections of order 1 (head) - Boundary conditions
                if obj.Load_Status(Npre)       
                    % Calculation of the determinant load upstream of the section, when there are discharges [mg]
                    obj.(['W_',Model])(Npre) = obj.(['BC_',Model,])(Npre) + obj.(['Load_',Model])(Npre);
                    % Change of status to no longer consider dumping in the accumulation
                    obj.Load_Status(Npre) = false;
                else
                    % Calculation of the determinant load upstream of the section, when there are no discharges [mg]
                    obj.(['W_',Model])(Npre)  = obj.(['BC_',Model])(Npre);
                end
            else  
                % Check for Posi greater than two
                if length(Posi)~=1
                    return
                end
                % Calculation of loads for order sections >1
                if obj.Load_Status(Npre)    
                    % Calculation of the determinant load upstream of the section, when there are discharges [mg]
                    obj.(['W_',Model])(Npre) = obj.(['W_',Model])(Npre) + (obj.(['C_',Model])(Posi).*obj.Q(Posi).*obj.ConFac_Q) + obj.(['Load_',Model])(Npre);
                    % Change of status to no longer consider the dumping in the accumulation
                    obj.Load_Status(Npre) = false;
                else
                    % Calculation of the determinant load upstream of the section, when there are no discharges [mg]
                    obj.(['W_',Model])(Npre) = obj.(['W_',Model])(Npre) + (obj.(['C_',Model])(Posi).*obj.Q(Posi).*obj.ConFac_Q);
                end  
            end

            % Calculation of the determinant concentration in the section [mg/l]
            obj.(['C_',Model])(Npre) = obj.(['W_',Model])(Npre)./obj.(['AF_',Model])(Npre);
        end
    end 
    
    % --------------------------------------------------------------------- 
    % Cumulative function for determinants that depend on the concentration 
    % of other determinants 
    % ---------------------------------------------------------------------
    methods
        function FunctionNetwork_CoDependent(obj, Npre, Posi,varargin)
            % This function performs the cumulative calculation of the 
            % concentrations, loads and assimilation factors of the water 
            % quality determinants that are dependent, in the cumulative 
            % scheme of the Functional Branch.
            
            % Check-in
            ip = inputParser;
            % Acronym of the model to be used. By default, 
            % Ammoniacal Nitrogen is assigned.
            addParameter(ip,'Model','NH4',@ischar)
            parse(ip,varargin{:})
            Model = ip.Results.Model;
            
            % Multiplying factor of the decay constant 1
            Factor  = 1;
            % Multiplying factor of the decay constant 2
            Factor2 = 1;
                
            % Parameters according to the determinant
            if strcmp(Model,'NH4')
                % Determinant on which Ammoniacal Nitrogen depends
                Model_d = 'NO';
            elseif strcmp(Model,'NO3')
                % Determinant on which Nitrates depend
                Model_d = 'NH4';
            elseif strcmp(Model,'PI')
                % Determinant on which inorganic phosphorus depends
                Model_d = 'PO';
            elseif strcmp(Model,'OM')
                % Determinant on which Organic Matter depends
                Model_d = 'NO3';
            elseif strcmp(Model,'Hg0')
                % Determinant on which Elemental Mercury
                Model_d = 'Hg2';
            elseif strcmp(Model,'Hg2')
                % Determinant on which Divalend Mercury
                Model_d = 'Hg0';
            elseif strcmp(Model,'MeHg')
                % Determinant on which Methyl Mercury
                Model_d = 'Hg2';
            end
            
            % Create a function handler map
            MethodMap1 = containers.Map(...
                        {'NH4', 'NO3', 'PI', 'OM', 'Hg0', 'Hg2', 'MeHg'},...
                        {@obj.Params_K__NH4, @obj.Params_K__NO3,...
                         @obj.Params_K__PI,  @obj.Params_K__OM,...
                         @obj.Params_K__Hg0, @obj.Params_K__Hg2,...
                         @obj.Params_K__MeHg});
            
            % Get the function handle from the mapper
            ParamMethod = MethodMap1(Model);

            % Calculation of determinant loads
            if isempty(Posi)        
                % Calculation of loads for sections of order 1 (head) - Boundary conditions
                if obj.Load_Status(Npre)   
                    % Calculation of the determinant load upstream of the section, when there are discharges [mg]
                    obj.(['W_',Model])(Npre)  = (obj.(['BC_',Model])(Npre)) + (obj.(['Load_',Model])(Npre));
                    % Change of status to no longer consider the dumping in the accumulation
                    obj.Load_Status(Npre)    = false;
                else
                    % Calculation of the determinant load upstream of the section, when there are no discharges [mg]
                    obj.(['W_',Model])(Npre)  = (obj.(['BC_',Model])(Npre));
                end  
                % Accumulation of flows for calculation of concentrations upstream of the section [m3/s]
                obj.Qaccum(Npre)        = obj.Q(Npre);
            else  
                % Check for Posi greater than two
                if length(Posi)~=1
                    return
                end
                % Calculation of loads for order sections >1
                if obj.Load_Status(Npre)  
                    % Calculation of the determinant load upstream of the section, when there are discharges [mg]
                    obj.(['W_',Model])(Npre) = obj.(['W_',Model])(Npre) + (obj.(['C_',Model])(Posi).*obj.Q(Posi).*obj.ConFac_Q) + (obj.(['Load_',Model])(Npre));
                    % Change of status to no longer consider the dumping in the accumulation
                    obj.Load_Status(Npre)    = false;
                else
                    % Calculation of the determinant load upstream of the section, when there are no discharges [mg]
                    obj.(['W_',Model])(Npre) = obj.(['W_',Model])(Npre) + (obj.(['C_',Model])(Posi).*obj.Q(Posi).*obj.ConFac_Q);
                end
                % Accumulation of flows for calculation of concentrations upstream of the section [m3/s]
                obj.Qaccum(Npre)        = obj.Qaccum(Npre) + obj.Q(Posi);
            end

            % Calculation of upstream concentration of the determinant [mg]
            Cu_Tmp = obj.(['W_',Model])(Npre)./(obj.Qaccum(Npre).*obj.ConFac_Q);
            
            % Calculation of decay constant 1 of determinant [1/d]
            obj.(['K__',Model])(Npre) = ParamMethod(obj.(['C_',Model_d])(Npre),Cu_Tmp,obj.(['K_',Model_d])(Npre),obj.(['K_',Model])(Npre));
                
            % Calculation of the assimilation factor of the determinant [Ad]
            obj.(['AF_',Model])(Npre) = obj.Cal_AF_CoDependent((obj.Q(Npre).*obj.ConFac_Q),obj.tao(Npre),obj.Tr(Npre),obj.(['K_',Model])(Npre),obj.(['K__',Model])(Npre),Factor,Factor2);

            % Calculation of the determinant concentration [mg/l]
            obj.(['C_',Model])(Npre) = obj.(['W_',Model])(Npre)./obj.(['AF_',Model])(Npre);
        end
    end        
    
    % --------------------------------------------------------------------- 
    % Cumulative function for determinants that do not depend on the 
    % concentration of other determinants 
    % ---------------------------------------------------------------------
    methods
        function FunctionNetwork_Temperature(obj, Npre, Posi)
            % This function performs the cumulative calculation of the 
            % concentrations, loads and assimilation factors of the water 
            % quality determinants that are independent, in the cumulative 
            % scheme of the Functional Branch. Model [String]: Acronym of 
            % the model to be used. By default, it is assigned a          
            
            % Calculation of determinant loads
            if isempty(Posi)        
                % Calculation of loads for sections of order 1 (head) - Boundary conditions
                if obj.Load_Status(Npre)&&(~isnan(obj.Load_T(Npre))&&(obj.Qwwd(Npre)>0)&&~isnan(obj.Qwwd(Npre)))
                    % Calculation of the determinant load upstream of the section, when there are discharges
                    obj.C_T(Npre) = ((obj.Qr(Npre)./obj.Qwwd(Npre)).*obj.C_T(Npre) + obj.Load_T(Npre))./(1+(obj.Qr(Npre)./obj.Qwwd(Npre)));
                    % Change of status to no longer consider dumping in the accumulation
                    obj.Load_Status(Npre) = false;
                else
                    % Calculation of the determinant load upstream of the section, when there are no discharges 
                    obj.C_T(Npre)  = obj.T(Npre);
                end                
            else  
                % Check for Posi greater than two
                if length(Posi)~=1
                    return
                end                
                IdT = (obj.ToNode == obj.FromNode(Npre))&obj.QtmpT;
                % Calculation of loads for order sections > 1
                if obj.Load_Status(Npre)&&(~isnan(obj.Load_T(Npre))&&(obj.Qwwd(Npre)>0)&&~isnan(obj.Qwwd(Npre)))  
                    % Calculation of the determinant load upstream of the section, when there are discharges
                    obj.C_T(Npre) = ((obj.Qr(Npre)./obj.Qwwd(Npre)).*obj.C_T(Npre) + obj.Load_T(Npre))./(1+(obj.Qr(Npre)./obj.Qwwd(Npre)));
                    % Calculation dischage
                    Qtmp = obj.Q(Npre) - sum(obj.Q(IdT));
                    % Calculation of the determinant load upstream of the section, when there are discharges                    
                    obj.C_T(Npre) = ((Qtmp./obj.Q(Posi)).*obj.C_T(Npre) + obj.C_T(Posi))./(1+(Qtmp./obj.Q(Posi)));
                    % Change of status to no longer consider the dumping in the accumulation
                    obj.Load_Status(Npre) = false;
                else
                    % Calculation dischage
                    Qtmp = obj.Q(Npre) - sum(obj.Q(IdT));
                    % Calculation of the determinant load upstream of the section, when there are no discharges
                    obj.C_T(Npre) = ((Qtmp./obj.Q(Posi)).*obj.C_T(Npre) + obj.C_T(Posi))./(1+(Qtmp./obj.Q(Posi)));
                end  
                obj.QtmpT(Posi) = false;
            end
        end
    end 

    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Temperature [SS] 
    % ---------------------------------------------------------------------
    methods                       
        function WQS_T(obj,varargin) 
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Temperature Load [°C]
            addParameter(ip, 'Load_T',zeros(size(obj.ReachID)),@ismatrix)
            % Temperature in Boundary Conditions [°C]
            addParameter(ip, 'BC_T',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Suspended Solids Load [mg/day]
            obj.Load_T      = ip.Results.Load_T;
            % Suspended Solids Load in Boundary Conditions [mg/day]
            obj.BC_T        = ip.Results.BC_T;
            % This function applies the total temperature model 
            % Initializing loads to 0
            obj.W_T         = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_T         = obj.T;
            % Initializing assimilation factors to 0
            obj.AF_T        = zeros(size(obj.ReachID));
            % Initializing accumulated flows to 0
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of assimilation factor
            obj.AF_T        = obj.Q.*obj.ConFac_Q;
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            obj.QtmpT       = true(size(obj.ReachID));
            % Function for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_Temperature(Npre, Posi(i));';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_Temperature(Npre, Posi), end;';            
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_T = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the temperature model',toc);
            fprintf(WQSMessage)
        end
    end

    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Suspended Solids [SS] 
    % ---------------------------------------------------------------------
    properties
        % Suspended solids decay rate by sedimentation [1/d]
        K_SS
    end
    
    methods               
        function Params_KSS(obj,Vss)
            % This function calculates the rate of decay of total suspended 
            % solids by sedimentation 
            % H         [m]     : Depth 
            % Vss       [m/d]   : Solids settling velocity 
            % K_SS      [1/d]   : Suspended solids decay constant
            
            obj.K_SS   = (Vss./obj.H);
        end
        
        function WQS_SS(obj,varargin) 
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Suspended solids sedimentation rate Vss [m/d] 
            % The default value is the value indicated in Qual2Kw [0.1]
            addParameter(ip, 'Vss',0.1,@isnumeric)
            % Suspended Solids Load [mg/day]
            addParameter(ip, 'Load_SS',zeros(size(obj.ReachID)),@ismatrix)
            % Suspended Solids Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_SS',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Suspended solids sedimentation rate Vss [m/d] 
            Vss             = ip.Results.Vss;
            % Suspended Solids Load [mg/day]
            obj.Load_SS     = ip.Results.Load_SS;
            % Suspended Solids Load in Boundary Conditions [mg/day]
            obj.BC_SS       = ip.Results.BC_SS;
            % This function applies the total suspended solids model 
            % Initializing loads to 0
            obj.W_SS        = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_SS        = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_SS       = zeros(size(obj.ReachID));
            % Initializing accumulated flows to 0
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Decay rate calculation [1/d]
            obj.Params_KSS(Vss)
            % Calculation of assimilation factor
            obj.AF_SS       = obj.Cal_AF_CoIndependent(obj.Q.*obj.ConFac_Q,obj.tao,obj.Tr,obj.K_SS);
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Function for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoIndependent(Npre, Posi(i),''Model'',''SS'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoIndependent(Npre, Posi,''Model'',''SS''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_SS = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Estimate turbidity
            obj.Cal_Turbidity;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the suspended solids model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Total Coliforms 
    % ---------------------------------------------------------------------
    properties
        % Total coliform decay rate by death and sedimentation [1/d]
        K_X 
    end
            
    methods                                    
        function Params_KX(obj,vx,Fpx,Kdx)
            % This function calculates the rate of decay of total coliforms by death 
            % and sedimentation 
            % H     [m]     : Depth 
            % Kd    [1/d]   : Reaction rates of total coliforms 
            % vx    [m/d]   : Sedimentation velocity of pathogens 
            % Fpx   [Ad]    : Fraction of bacteria adsorbed to solid particles 
            % Kdx   [1/d]   : Reaction rate of total coliforms 
            % K_X   [1/d]   : Decay rate of total coliforms

            % Total coliform decay rate by death and sedimentation [1/d] 
            % Coliforms outside the digestive tract show a decay expressed 
            % by 3 components (natural mortality, radiation and sedimentation) 
            % whose decay (kcoli) has been found between 5.66 and 10.32 d-1 
            % (Torres-Matallana, 2009). Accordingly, since decay by 
            % radiation is not being considered, the minimum value is 
            % considered to be 5.66 and the sedimentation and precipitation 
            % rates represent 50% of the death of the pathogens.
            obj.K_X    = 5.66 + (Kdx + ((Fpx.*vx)./obj.H));
        end
        
        function WQS_X(obj,varargin)
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Pathogen sedimentation rate [m/d] Default value indicated in Qual2Kw [1] is taken
            addParameter(ip,'vx',1,@isnumeric)
            % Adsorbed fraction of bacteria on solid particles Considering 
            % what was indicated by Rojas (2011) and in the document 
            % published by MADS and ANLA (2013), a typical value of 0.7 
            % can be adopted
            addParameter(ip, 'Fpx', 0.7, @isnumeric)
            % Calculation of total coliform reaction rate [1/d]
            % k = 0.8
            addParameter(ip, 'Kdx', obj.Kd_ArrheniusModel('NameParam','X','k',0.8), @ismatrix)
            % Total Coliform Load [UFC/day]
            addParameter(ip, 'Load_X',zeros(size(obj.ReachID)),@ismatrix)
            % Total Coliform Load in Boundary Conditions [UFC/day]
            addParameter(ip, 'BC_X',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Pathogen sedimentation rate [m/d]
            vx              = ip.Results.vx;
            % Fraction of bacteria adsorbed to solid particles
            Fpx             = ip.Results.Fpx;
            % Calculation of total coliform reaction rate [1/d]
            Kdx             = ip.Results.Kdx;
            % Total Coliform Load [UFC/day]
            obj.Load_X      = ip.Results.Load_X;
            % Total Coliform Load in Boundary Conditions [UFC/day]
            obj.BC_X        = ip.Results.BC_X;
            % This function applies the total coliform model Initializing loads to 0
            obj.W_X         = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_X         = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_X        = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Decay rate calculation [1/d]
            obj.Params_KX(vx,Fpx,Kdx);
            % Calculation of assimilation factor
            obj.AF_X        = obj.Cal_AF_CoIndependent(obj.Q.*obj.ConFac_Q,obj.tao,obj.Tr,obj.K_X);
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Function for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoIndependent(Npre, Posi(i),''Model'',''X'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoIndependent(Npre, Posi,''Model'',''X''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_X = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results            
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the total coliforms model',toc);
            fprintf(WQSMessage)
        end
        
    end   
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Organic Nitrogen 
    % ---------------------------------------------------------------------
    properties
        % Organic Nitrogen reaction rate [1/d]
        K_NO
        % Organic Nitrogen decay rate by hydrolysis and sedimentation [1/d]
        K__NO      
    end
    
    methods              
        function Params_KNO(obj,vno)
            % This function calculates the rate of decay of organic nitrogen 
            % by hydrolysis 
            % H     [m]   : Depth 
            % K_NO  [1/d] : Rate of reaction by hydrolysis 
            % vno   [m/d] : Sedimentation velocity of organic nitrogen 
            % K__NO [1/d] : Rate of decay of organic nitrogen by hydrolysis and sedimentation

            % Organic Nitrogen decay rate by hydrolysis and sedimentation [1/d]
            obj.K__NO   = obj.K_NO + (vno./obj.H);
        end
        
        function WQS_NO(obj,varargin)   
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Organic nitrogen sedimentation rate [m/d] The default value 
            % indicated in Qual2Kw [0.0005] is taken
            addParameter(ip, 'vno',0.0005 ,@isnumeric)
            % Hydrolysis reaction rate [1/d]. k = 0.02
            addParameter(ip, 'khno',obj.Kd_ArrheniusModel('NameParam','NO','k',0.02) ,@ismatrix)
            % Organic Nitrogen Load [mg/day]
            addParameter(ip, 'Load_NO',zeros(size(obj.ReachID)),@ismatrix)
            % Organic Nitrogen Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_NO',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Organic nitrogen sedimentation rate [m/d]
            vno             = ip.Results.vno;
            % Hydrolysis reaction rate [1/d]
            obj.K_NO        = ip.Results.khno;
            % Organic Nitrogen Load [mg/day]
            obj.Load_NO     = ip.Results.Load_NO;
            % Organic Nitrogen Load in Boundary Conditions [mg/day]
            obj.BC_NO       = ip.Results.BC_NO;
            % This function applies the organic nitrogen model Initializing loads to 0
            obj.W_NO        = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_NO        = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_NO       = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Decay rate calculation [1/d]
            obj.Params_KNO(vno);
            % Calculation of assimilation factor
            obj.AF_NO       = obj.Cal_AF_CoIndependent(obj.Q.*obj.ConFac_Q,obj.tao,obj.Tr,obj.K_NO);
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Function for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoIndependent(Npre, Posi(i),''Model'',''NO'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoIndependent(Npre, Posi,''Model'',''NO''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_NO = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the organic nitrogen model',toc);
            fprintf(WQSMessage)
        end
    end       
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Ammoniacal Nitrogen 
    % ---------------------------------------------------------------------
    properties        
        % Ammoniacal Nitrogen decay rate by nitrification [1/d]
        K_NH4
        % Ammoniacal Nitrogen decay rate combining Hydrolysis of organic 
        % nitrogen and Nitrification of ammoniacal nitrogen [1/d]
        K__NH4
    end
    
    methods               
        function K_NH4 = Params_KNH4(obj)
            % This function calculates the nitrification decay rate 
            % H     [m]     : Depth 
            % v     [m/s]   : Current velocity 
            % K_NH4 [1/d]   : Nitrification decay rate    

            % Robles and Camacho (2005) and Medina and Camacho (2008) 
            % demonstrated that the Bansal equation does not have a good 
            % performance for mountain rivers. In their studies, 
            % Robles and Camacho (2005), using the Couchaine method, 
            % proposed the following expression for mountain rivers:
            K_NH4_M  = ((0.4381.*(obj.v./obj.H)) + 0.5394);

            % For plain rivers, the nitrification decay rate is estimated 
            % according to the proposal of Bansal (1976) as:
            kv      = 0.0000010034; % water kinematic viscosity - 20 °C
            Part_1  = log(sqrt(9.81.*(obj.H.^3))./kv).^1.36;
            K_NH4_F = (10.*(-3.421 + Part_1)*kv)./(obj.H.^2);
            
            % Nitrification decay rate
            K_NH4   = (K_NH4_M.*obj.RiverType) + (K_NH4_F.*~obj.RiverType);
            % Correction for temperature
            K_NH4 = obj.Kd_ArrheniusModel('NameParam','NH4','k',K_NH4);
        end
        
        function K__NH4 = Params_K__NH4(~,C_NO,NH4u,K_NO,K_NH4)
             % This function calculates the decay rate of Ammoniacal 
             % Nitrogen combining Hydrolysis and Nitrification [1/d] 
             % C_NO     [mg/l]  : Organic Nitrogen concentration 
             % NH4u     [mg/l]  : Upstream Ammoniacal Nitrogen concentration 
             % K_NO     [1/d]   : Organic Nitrogen decay rate by Hydrolysis 
             % K_NH4    [1/d]   : Ammoniacal Nitrogen decay rate by nitrification 
             % K__NH4   [1/d]   : Ammoniacal Nitrogen decay rate combining Hydrolysis 
             %               of organic nitrogen and Nitrification of Ammoniacal Nitrogen
            
            % Controls for when concentrations are zero
            if NH4u == 0
                Var = 1;
            elseif C_NO == 0
                Var = 0;
            else
                Var = (C_NO./NH4u);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Ammoniacal Nitrogen decay rate combining Hydrolysis of organic 
            % nitrogen and Nitrification of ammoniacal nitrogen [1/d]
            K__NH4  = (K_NO*Var) - K_NH4;
        end
                        
        function WQS_NH4(obj,varargin)
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % nitrification decay rate [1/d]
            addParameter(ip, 'knh4',obj.Params_KNH4,@ismatrix)
            % Ammoniacal Nitrogen Load [mg/day]
            addParameter(ip, 'Load_NH4',zeros(size(obj.ReachID)),@ismatrix)
            % Ammoniacal Nitrogen Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_NH4',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Ammoniacal Nitrogen Load [mg/day]
            obj.Load_NH4    = ip.Results.Load_NH4;
            % Ammoniacal Nitrogen Load in Boundary Conditions [mg/day]
            obj.BC_NH4      = ip.Results.BC_NH4;
            % This function applies the ammonia nitrogen model Initializing loads to 0
            obj.W_NH4       = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_NH4       = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_NH4      = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of nitrification decay rate
            obj.K_NH4       = ip.Results.knh4;
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''NH4'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''NH4''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_NH4 = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the ammoniacal nitrogen model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Nitrates 
    % ---------------------------------------------------------------------
    properties        
        % Nitrate decay rate by denitrification [1/d]
        K_NO3
        % Nitrate decay rate combining Ammonium Nitrogen Nitrification and 
        % Nitrate Denitrification [1/d]
        K__NO3
    end
    
    methods               
        function Params_KNO3(obj,FoxdNO3,FdNO3)
            % This function calculates the rate of decay of nitrates 
            % by nitrification with oxygen effect 
            % Foxd_NO3 [Ad] : Effect of low oxygen on denitrification 
            % Fd_NO3 [1/d] : Denitrification rate 
            % K_NO3 [1/d] : Denitrification decay rate with oxygen effect
                        
            % Nitrate decay rate by oxygen-induced nitrification
            obj.K_NO3   = obj.ReachID*0;
            obj.K_NO3(:)= FoxdNO3.*FdNO3;
        end
        
        function K__NO3 = Params_K__NO3(~,C_NH4,NO3u,K_NH4,K_NO3)
            % This function calculates the nitrate decay rate combining 
            % Ammonium Nitrogen Nitrification and Nitrate Denitrification 
            % C_NH4 [mg/l]  : Ammoniacal Nitrogen Concentration 
            % NO3u  [mg/l]  : Upstream Nitrate Concentration 
            % K_NH4 [1/d]   : Ammoniacal Nitrogen Decay Rate by Nitrification 
            % K_NO3 [1/d]   : Nitrate Decay Rate by Denitrification 
            % K__NO3[1/d]   : Nitrate Decay Rate combining Ammonium Nitrogen 
            %                 Nitrification and Nitrate Denitrification
            
            % Controls for when concentrations are zero
            if NO3u == 0
                Var = 1;
            elseif C_NH4 == 0
                Var = 0;
            else
                Var = (C_NH4./NO3u);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Calculation of nitrate decay rate combining Nitrification of 
            % ammonium nitrogen and Denitrification of nitrates [1/d]
            K__NO3    = (K_NH4*Var) - K_NO3;
        end
        
        function WQS_NO3(obj,varargin)
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Effect of low oxygen on denitrification
            addParameter(ip, 'FoxdNO3',exp(-0.60) ,@isnumeric)
            % The denitrification rate kd_NO3 is the one indicated in Qual2Kw as the default value [0.1]
            addParameter(ip, 'FdNO3',obj.Kd_ArrheniusModel('NameParam','NO3','k',0.1) ,@ismatrix)
            % Nitrates Load [mg/day]
            addParameter(ip, 'Load_NO3',zeros(size(obj.ReachID)),@ismatrix)
            % Load Check Nitrates in boundary conditions [mg/day]
            addParameter(ip, 'BC_NO3',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Effect of low oxygen on denitrification
            FoxdNO3         = ip.Results.FoxdNO3;
            % Denitrification rate
            FdNO3           = ip.Results.FdNO3;
            % Nitrates Load [mg/day]
            obj.Load_NO3    = ip.Results.Load_NO3;
            % Load Check Nitrates in boundary conditions [mg/day]
            obj.BC_NO3      = ip.Results.BC_NO3;
            % This function applies the ammonia nitrogen model Initializing loads to 0
            obj.W_NO3       = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_NO3       = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_NO3      = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of nitrification decay rate
            obj.Params_KNO3(FoxdNO3,FdNO3)
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''NO3'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''NO3''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_NO3 = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the nitrates model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Organic Phosphorus 
    % ---------------------------------------------------------------------
    properties
        % Rate of decay of organic phosphorus by hydrolysis and sedimentation [1/d]
        K_PO
    end
    
    methods                
        function Params_KPO(obj,Khpo,vpo)
             % This function calculates the decay constant of organic phosphorus 
             % by hydrolysis and sedimentation 
             % H [m] : Depth 
             % Khpo [1/d] : Rate of reaction by hydrolysis 
             % vpo [m/d] : Sedimentation velocity of organic phosphorus 
             % K_PO [1/d] : Rate of decay of organic phosphorus by hydrolysis and sedimentation
            
            % Decay constant of organic phosphorus by hydrolysis and sedimentation[1/d]
            obj.K_PO    = Khpo + (vpo./obj.H);
        end
        
        function WQS_PO(obj,varargin) 
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Organic phosphorus sedimentation rate (m/d) 
            % The default value indicated in Qual2Kw [0.001] is taken
            addParameter(ip, 'vpo',0.001,@isnumeric)
            % Hydrolysis reaction rate [1/d] 
            % Default value indicated in Qual2Kw [0.03] is taken
            addParameter(ip, 'Khpo',0.03,@isnumeric)
            % Organic Phosphorus Load [mg/day]
            addParameter(ip, 'Load_PO',zeros(size(obj.ReachID)),@ismatrix)
            % Inorganic Phosphorus Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_PO',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Organic phosphorus sedimentation rate [m/d]
            vpo             = ip.Results.vpo;
            % Hydrolysis reaction rate [1/d]
            Khpo            = ip.Results.Khpo;
            % Organic Phosphorus Load [mg/day]
            obj.Load_PO     = ip.Results.Load_PO;
            % Inorganic Phosphorus Load in Boundary Conditions [mg/day]
            obj.BC_PO       = ip.Results.BC_PO;
            % This function applies the organic phosphorus model Initializing charges to 0
            obj.W_PO        = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_PO        = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_PO       = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Decay rate calculation [1/d]
            obj.Params_KPO(Khpo,vpo);
            % Calculation of assimilation factor
            obj.AF_PO       = obj.Cal_AF_CoIndependent(obj.Q.*obj.ConFac_Q,obj.tao,obj.Tr,obj.K_PO);
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Function for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoIndependent(Npre, Posi(i),''Model'',''PO'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoIndependent(Npre, Posi,''Model'',''PO''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_PO = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the organic phosphorus model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Inorganic Phosphorus 
    % ---------------------------------------------------------------------
    properties        
        % Decay rate of inorganic phosphorus by sedimentation [1/d]
        K_PI
        % Decay rate of inorganic phosphorus combining hydrolysis of organic 
        % phosphorus and sedimentation of inorganic phosphorus [1/d]
        K__PI
    end
    
    methods               
        function Params_KPI(obj,vpi)
             % This function calculates the rate of decay of inorganic phosphorus 
             % by sedimentation [1/d] 
             % H    [m]     : Depth 
             % vpi  [m/d]   : Sedimentation velocity of inorganic phosphorus 
             % K_PI [1/d]   : Decay rate of inorganic phosphorus by sedimentation
            
            % Decay rate of inorganic phosphorus by sedimentation [1/d]
            obj.K_PI    = (vpi./obj.H);
        end
        
        function K__PI = Params_K__PI(~,C_PO,PIu,K_PO,K_PI)                        
             % This function calculates the decay rate of inorganic phosphorus 
             % combining hydrolysis of organic phosphorus and sedimentation 
             % of inorganic phosphorus [1/d] 
             % C_PO [mg/l]  : Organic phosphorus concentration 
             % PIu  [mg/l]  : Upstream inorganic phosphorus concentration 
             % K_PO [1/d]   : Decay rates of organic phosphorus by hydrolysis 
             % K_PI [1/d]   : Decay rate of inorganic phosphorus by sedimentation 
             % K__PI[1/d]   : Decay rate of inorganic phosphorus combining 
             % hydrolysis of organic phosphorus and sedimentation of inorganic phosphorus
            
            % Controls for when concentrations are zero
            if PIu == 0
                Var = 1;
            elseif C_PO == 0
                Var = 0;
            else
                Var = (C_PO./PIu);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Calculation of inorganic phosphorus decay rate combining organic 
            % phosphorus hydrolysis and inorganic phosphorus sedimentation [1/d]
            K__PI    = (K_PO*Var) - K_PI;
        end
        
        function WQS_PI(obj,varargin)    
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Inorganic phosphorus sedimentation rate [m/d] Default value indicated in Qual2Kw [0.8] is taken
            addParameter(ip, 'vpi',0.8,@isnumeric)
            % Inorganic Phosphorus Load [mg/day]
            addParameter(ip, 'Load_PI',zeros(size(obj.ReachID)),@ismatrix)
            % Organic Phosphorus Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_PI',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Inorganic phosphorus sedimentation rate [m/d]
            vpi             = ip.Results.vpi;
            % Inorganic Phosphorus Load [mg/day]
            obj.Load_PI     = ip.Results.Load_PI;
            % Organic Phosphorus Load in Boundary Conditions [mg/day]
            obj.BC_PI       = ip.Results.BC_PI;
            % This function applies the inorganic phosphorus model Initializing charges to 0
            obj.W_PI        = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_PI        = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_PI       = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of nitrification decay rate
            obj.Params_KPI(vpi);
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''PI'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''PI''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_PI = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the inorganic phosphorus model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % Water Quality Simulation - Organic matter 
    % ---------------------------------------------------------------------
    properties
        % Rate of decay by decomposition of organic matter [1/d]
        KdOM
        % Nitrate nitrification decay rate considering oxygen [1/d]
        K_NO3_OM
        % Rate of decay of organic matter by oxidation [1/d]
        K_OM        
        % organic matter decay rate considering nitrification of nitrates 
        % and oxidation of organic matter [1/d]
        K__OM
        % Temporary variable to store data
        VarTmp        
    end
    
    methods   
        function KdOM = KdOM_WrightMcDonnell(obj)
            % This function calculates the rate of decay by decomposition 
            % of organic matter. The decomposition rate can be estimated 
            % as a function of the flow rate following the expression 
            % proposed by Wright and McDonnell (1979) 
            % Q     [m3/s]  : Flow rate 
            % KdMo  [1/d]   : Decay rate by decomposition
            
            % The range of application of this expression is for flow rates between 0.3 and 23 m3/s.
            id              = (obj.Q <= 23);
            % Rate of decay by decomposition of organic matter [1/d]
            KdOM        = obj.ReachID.*0;
            KdOM(id)    = 1.796.*(obj.Q(id).^-0.49);
            % Above the highest flow rate, Wright and McDonnell showed that 
            % kdMo is independent of flow rate and the rates found are 
            % consistent with laboratory rates. Above the range of the equation, 
            % for practical purposes, a constant value of 0.30 1/d may be assumed. 
            % The maximum value of kdMo is suggested to be 3.5 1/d.
            KdOM(~id)   = 0.3;  
            % Correction for temperature
            KdOM = obj.Kd_ArrheniusModel('NameParam','OM','k',KdOM);
            KdOM(KdOM > 3.5) = 3.5;            
        end
        
        function Cal_ro_OM(obj, FoxdOM)
             % This function calculates the decay constant of organic matter 
             % by oxidation FoxdOM [Ad] : Factor that considers the effect of oxygen. 
             % KdMo [1/d] : Decay rate by decomposition 
             % K_OM [1/d] : Decay rate by decomposition considering oxygen
             
            % Decay rate by decomposition considering oxygen
            obj.K_OM   = FoxdOM.*obj.KdOM;
        end
        
        function Cal_teta_OM(obj,FoxdNO3,FdNO3)    
            % This function calculates the nitrification decay rate 
            % considering oxygen 
            % FoxdNO3   [Ad]    : Factor that considers the effect of oxygen 
            % FdNO3     [1/d]   : Nitrate nitrification decay rate considering oxygen 
            % K_NO3_OM  [1/d]   : Nitration decay rate considering oxygen

            % nitrification decay rate considering oxygen [1/d]
            obj.K_NO3_OM    = obj.ReachID.*0;
            obj.K_NO3_OM(:) = (0.00286.*(1-FoxdNO3).*FdNO3);
        end
        
        function K__OM = Params_K__OM(~,C_NO3,OMu,K_NO3_OM,K_OM)
            % This function calculates the organic matter decay rate 
            % considering nitrate nitrification and organic matter oxidation 
            % C_NO3     [mg/l] : Nitrate concentration 
            % OMu       [mg/l] : Upstream organic matter concentration 
            % K_NO3_OM  [1/d]  : Nitrate decay rates by nitrification considering oxygen 
            % K_OM      [1/d]  : Organic matter decay rate by decomposition 
            % K__OM     [1/d]  : Inorganic phosphorus decay rate combining organic 
            %                    phosphorus hydrolysis and inorganic phosphorus sedimentation

            % Controls for when concentrations are zero
            if OMu == 0
                Var = 1;
            elseif C_NO3 == 0
                Var = 0;
            else
                Var = (C_NO3./OMu);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Organic matter decay rate considering nitrification of nitrates and oxidation of organic matter
            K__OM   = -(K_NO3_OM*Var) - K_OM;
        end
        
        function WQS_OM(obj,varargin)
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Calculation of decay rate by decomposition following the expression proposed by Wright and McDonnell (1979) [1/d]
            addParameter(ip, 'kdOM',obj.KdOM_WrightMcDonnell,@ismatrix)
            % Factor considering the effect of oxygen [Ad]
            addParameter(ip, 'FoxdOM',1 - exp(-0.6),@isnumeric)
            % Effect of low oxygen on denitrification
            addParameter(ip, 'FoxdNO3',exp(-0.60) ,@isnumeric)
            % The denitrification rate kd_NO3 is the one indicated in Qual2Kw as the default value [0.1]
            addParameter(ip, 'FdNO3',0.1 ,@isnumeric)
            % Organic Matter Load Check [mg/day]
            addParameter(ip, 'Load_OM',zeros(size(obj.ReachID)),@ismatrix)
            % Check of Organic Matter Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_OM',zeros(size(obj.ReachID)),@ismatrix)           
            % Checking input data
            parse(ip,varargin{:})
            % Inorganic phosphorus sedimentation rate [m/d]
            obj.KdOM        = ip.Results.kdOM;
            % Factor considering the effect of oxygen [Ad]
            FoxdOM          = ip.Results.FoxdOM;
            % Effect of low oxygen on denitrification
            FoxdNO3         = ip.Results.FoxdNO3;
            % Denitrification rate
            FdNO3           = ip.Results.FdNO3;
            % Organic Matter Load Check [mg/day]
            obj.Load_OM     = ip.Results.Load_OM;
            % Check of Organic Matter Load in Boundary Conditions [mg/day]
            obj.BC_OM       = ip.Results.BC_OM;
            % This function applies the inorganic phosphorus model Initializing charges to 0
            obj.W_OM        = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_OM        = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_OM       = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Combined decay rate initialization
            obj.K__OM       = zeros(size(obj.ReachID));
            % Assign discharge status
            obj.Update_VerStatus;
            % Calculation parameters
            obj.Cal_ro_OM(FoxdOM);
            obj.Cal_teta_OM(FoxdNO3,FdNO3);
            % Save original data of K_NO3
            obj.VarTmp      = obj.K_NO3;
            obj.K_NO3       = obj.K_NO3_OM;            
            % Activate function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''OM'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''OM''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable function scheme
            obj.StatusFun   = false;
            obj.StatusWQS_OM = true;
            % Assign discharge status
            obj.Update_VerStatus;
            % Assign original values ​​of KNO3
            obj.K_NO3 = obj.VarTmp;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the organic matter model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % WQS - Dissolved oxygen 
    % ---------------------------------------------------------------------
    properties
        % Reaeration rate [1/d]
        Ka_DO
        % Dissolved oxygen deficit decay rate [1/d]
        K_DOD
        % 
        K_NO3_DO
        % Dissolved oxygen deficit decay rate according to ammonia nitrogen and organic matter content
        K__DOD
    end
    
    methods            
        function Ka_DO = Ka_DO_Model(obj)
            % This function calculates the re-aeration rate according to 
            % the expression formulated by Tsivoglou and Neal (1976) 
            % for mountain rivers and those of O’Connor – Dubbins, 
            % Churchill and Owens - Gibbs for plain rivers 
            % H     [m]     : Depth of the water table of the stream section 
            % v     [m/s]   : Flow velocity of the stream section 
            % Q     [m3/s]  : Flow rate of the stream section 
            % S     [m/m]   : Slope of the stream section 
            % Ka_DO [1/d]   : Re-aeration rate
            
            Ka_DO       = obj.ReachID.*NaN;
            % Mountain rivers Formulated by Tsivoglou and Neal (1976).
            id              = obj.RiverType&...
                              (0.0283<=obj.Q)&...
                              (obj.Q < 0.4247);
            Ka_DO(id)   = 31.183.*obj.v(id).*obj.S(id);
            id              = obj.RiverType&...
                              (0.4247<=obj.Q)&...
                              (obj.Q < 84.938);
            Ka_DO(id)   = 15.308.*obj.v(id).*obj.S(id);

            % Rivers of the O’Connor-Dubbins Plain
            id              = ~obj.RiverType&...
                              (0.30<=obj.H)&(obj.H < 9.14)&...
                              (0.15<=obj.v)&(obj.v < 0.49);
            Ka_DO(id)   = 3.95.*((obj.v(id).^0.5)./(obj.H(id).^1.5));
            
            % Churchill
            id              = ~obj.RiverType&...
                             (0.31<=obj.H)&(obj.H < 3.35)&...
                             (0.55<=obj.v)&(obj.v < 1.52);
            Ka_DO(id)   = 5.026.*(obj.v(id)./(obj.H(id).^1.67));
            
            % Owens - Gibbs
            id              = ~obj.RiverType&...
                              (0.12<=obj.H)&(obj.H < 0.73)&...
                              (0.03<=obj.v)&(obj.v < 0.55);
            Ka_DO(id)   = 5.32.*((obj.v(id).^0.67)./(obj.H(id).^1.85));
            
            % Thackston & Dawson, 2001 Acceleration of gravity [m/s^2]
            g   = 9.81;
            % Hydraulic radius (assuming rectangular channel)
            Rn  = (obj.W.*obj.H)./(obj.W + (2*obj.H));
            % Cutting speed
            Vc = sqrt(g.*Rn.*obj.S);
            % Froude number
            Fr = (obj.v.^2)./(g.*obj.L);
            % Re-aeration rate calculation
            id = isnan(Ka_DO);
            Ka_DO(id) = 2.16.*( 1 + (9.*(Fr(id).^0.25))).*(Vc(id)./obj.H(id));
            % Correction for temperature
            Ka_DO = obj.Kd_ArrheniusModel('NameParam','DO','k',Ka_DO);
        end
        
        function K__DOD = Params_K__DOD(~,C_OM,C_NH4,DOu,KdOM,K_NH4,Ka)
            % This function calculates the rate of decay of dissolved oxygen 
            % deficit according to ammonia nitrogen and organic matter content 
            % C_OM  [mg/l]  : Concentration of organic matter 
            % C_NH4 [mg/l]  : Concentration of ammonia nitrogen 
            % DOu   [mg/l]  : Concentration of dissolved oxygen 
            % KdMo  [1/d]   : Rate of decay by decomposition of organic matter 
            % K_NH4 [1/d]   : Decay rates of ammonia nitrogen by nitrification considering oxygen 
            % Ka    [1/d]   : Reaeration rate 
            % K__DOD[1/d]   : Rate of decay of dissolved oxygen deficit according 
            %                 to ammonia nitrogen and organic matter content
            
            % Controls for when organic matter concentrations are zero
            if DOu == 0
                Var_1 = 1;
            elseif C_OM == 0
                Var_1 = 0;
            else
                Var_1 = (C_OM./DOu);
            end
            
            % Control so that no determinant matter is generated
            if Var_1 > 1
                Var_1 = 1;
            end
            
            % Controls for when ammonia nitrogen concentrations are zero
            if DOu == 0
                Var_2 = 1;
            elseif C_NH4 == 0
                Var_2 = 0;
            else
                Var_2 = (C_NH4./DOu);
            end
            
            % Control so that no determinant matter is generated
            if Var_2 > 1
                Var_2 = 1;
            end
            
            % Dissolved oxygen deficit decay rate according to ammonia nitrogen 
            % and organic matter content [1/d]
            K__DOD   = (KdOM*Var_1) + (4.57.*K_NH4*Var_2) - Ka;
        end
        
        function AF = Cal_AF_DOD(~,Q,C_DOu,C_OS,tao,Tr,Ka,Ko)
            % This function estimates the assimilation factor for dissolved 
            % oxygen deficit 
            % DOu   [mg/l]  : Dissolved oxygen upstream 
            % Os    [mg/l]  : Saturation oxygen 
            % tao   [d]     : Arrival or advection time 
            % Tr    [d]     : Dead zone residence time 
            % Ka    [/d]    : Reaeration constant 
            % Ko    [/d]    : Oxygen decay constant by consumption of organic 
            %                 matter and ammoniacal nitrogen 
            % Q     [lts]   : River section flow 
            % AF    [mg/l]  : Assimilation factor of the determinant
            
            Du      = C_OS - C_DOu;
            Part_1  = ((Tr.*Ka) + 1);
            Part_2  = (Du.*exp(Ko.*tao)) + (Tr.*C_DOu.*(Ka+Ko));
            AF      = Du.*Q.*(Part_1./Part_2);
        end
         
        function FunctionNetwork_CoDependent_DOD(obj, Npre, Posi,varargin)
            % This function performs the cumulative calculation of the 
            % concentration, load and assimilation factor of the dissolved 
            % oxygen deficit, in the cumulative scheme of the Functional Branch.
            
            % Calculation of dissolved oxygen loads
            if isempty(Posi)       
                % Saturation oxygen allocation in order 1 (header) sections - Boundary conditions [mg]
                obj.W_DO(Npre)  = obj.C_OS(Npre).*(obj.Q(Npre).*obj.ConFac_Q).*0.99;  
                % Accumulation of flows for calculation of concentrations upstream of the section [m3/s]
                obj.Qaccum(Npre)  = obj.Q(Npre);
            else  
                % Check for Posi greater than two
                if length(Posi)~=1
                    return
                end    
                % Calculation of loads for order sections >1 [mg]
                obj.W_DO(Npre) = obj.W_DO(Npre) + (obj.C_DO(Posi).*obj.Q(Posi).*obj.ConFac_Q);   
                % Accumulation of flows for calculation of concentrations upstream of the section [m3/s]
                obj.Qaccum(Npre) = obj.Qaccum(Npre) + obj.Q(Posi);
            end
            
            % Calculation of dissolved oxygen concentration [mg]
            obj.C_DO(Npre) = obj.W_DO(Npre)./(obj.Qaccum(Npre).*obj.ConFac_Q); 
            
            % Check if the dissolved oxygen concentration is higher than the saturation concentration [mg/l]
            if obj.C_DO(Npre) >= obj.C_OS(Npre)
                % C_DO can never be equal to C_OS
                obj.C_DO(Npre) = obj.C_OS(Npre) - 0.1; %[0.1 mg/l]
                % Recalculation of dissolved oxygen load
                obj.W_DO(Npre) = obj.C_DO(Npre).*(obj.Q(Npre).*obj.ConFac_Q);
            elseif obj.C_DO(Npre) <=0
                % [0.1 mg/l]
                obj.C_DO(Npre) = 0.1; 
                % Recalculation of dissolved oxygen load
                obj.W_DO(Npre) = obj.C_DO(Npre).*(obj.Q(Npre).*obj.ConFac_Q);
            end
            
            % Calculation of upstream dissolved oxygen concentration [mg]
            DOu                 = obj.C_DO(Npre);       
            
            % Calculation of decay constant [1/d]
            obj.K__DOD(Npre)     = obj.Params_K__DOD(obj.C_OM(Npre),obj.C_NH4(Npre),DOu,obj.KdOM(Npre),obj.K_NH4(Npre),obj.Ka_DO(Npre));
            
            % Calculation of the assimilation factor [Ad]
            obj.AF_DOD(Npre)  = obj.Cal_AF_DOD((obj.Q(Npre).*obj.ConFac_Q),DOu,obj.C_OS(Npre),obj.tao(Npre),obj.Tr(Npre),obj.Ka_DO(Npre),obj.K__DOD(Npre));            
            
            % Calculation of dissolved oxygen deficit concentration
            obj.C_DOD(Npre)      = obj.C_OS(Npre) - obj.C_DO(Npre);                        
            
            % Recalculation of dissolved oxygen deficit load [mg]
            obj.W_DOD(Npre)      = obj.C_DOD(Npre).*(obj.Qaccum(Npre).*obj.ConFac_Q);                        
            
            % Calculate dissolved oxygen deficit concentration [mg/l]
            obj.C_DOD(Npre)      = obj.W_DOD(Npre)./obj.AF_DOD(Npre);   
            
            % Check if the dissolved oxygen deficit is equal to zero or saturation oxygen
            if obj.C_DOD(Npre) <= 0
                % C_DOD should never be equal to zero [mg/l] [0.1 mg/l]
                obj.C_DOD(Npre) = 0.1;                 
            elseif obj.C_DOD(Npre) >= obj.C_OS(Npre)
                % C_DO should never be equal to zero [mg/l] [0.1 mg/l]
                obj.C_DOD(Npre) = obj.C_OS(Npre) - 0.1;           
            end
            
            % Calculate dissolved oxygen concentration [mg/l]
            obj.C_DO(Npre) = obj.C_OS(Npre) - obj.C_DOD(Npre);  
            
            % Calculate dissolved oxygen load
            obj.W_DO(Npre) = obj.C_DO(Npre).*(obj.Q(Npre).*obj.ConFac_Q);           
        end
        
        function WQS_DOD(obj,varargin)
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % nitrification decay rate [1/d]
            addParameter(ip, 'ka',obj.Ka_DO_Model,@ismatrix)
            % Dissolved Oxygen Deficit Load [mg/day]
            addParameter(ip, 'Load_DOD',zeros(size(obj.ReachID)),@ismatrix)
            % Dissolved Oxygen Deficit Load in boundary conditions [mg/day]
            addParameter(ip, 'BC_DOD',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Dissolved Oxygen Deficit Load [mg/day]
            obj.Load_DOD    = ip.Results.Load_DOD;
            % Dissolved Oxygen Deficit Load in boundary conditions [mg/day]
            obj.BC_DOD      = ip.Results.BC_DOD;

            % Initialize variables
            obj.Qaccum      = zeros(size(obj.ReachID));
            obj.W_DOD       = zeros(size(obj.ReachID));
            obj.C_DOD       = zeros(size(obj.ReachID));
            obj.W_DO        = zeros(size(obj.ReachID));
            obj.C_DO        = zeros(size(obj.ReachID));
            obj.AF_DOD      = zeros(size(obj.ReachID));
            obj.K__DOD      = zeros(size(obj.ReachID));
            % Assign discharge status
            obj.Update_VerStatus;
            % Re-aeration rate
            obj.Ka_DO       = ip.Results.ka;
            % Oxygen saturation
            obj.Cal_OS
            % Activate function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent_DOD(Npre, Posi(i));';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent_DOD(Npre, Posi), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable function scheme
            obj.StatusFun   = false;
            obj.StatusWQS_DO = true;
            % Correct sections so that they have saturation oxygen Assign discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the dissolved oxygen deficit model',toc);
            fprintf(WQSMessage)
        end
    end      
        
    % --------------------------------------------------------------------- 
    % WQS - Elemental Mercury 
    % ---------------------------------------------------------------------
    properties        
        % Decay rate of elemental mercury by volatilization and oxidation [1/d]
        K_Hg0
        % Elemental mercury decay rate combining reduction of divalent 
        % mercury and volatilization and oxidation of elemental mercury [1/d]
        K__Hg0
    end
    
    methods               
        function Params_KHg0_V1(obj,Vv,Kox,Fd)
            % Rate of decay of elemental mercury by volatilization and oxidation [1/d] 
            % H [m] : Depth 
            % Vv [m/d] : Volatilization rate of elemental mercury 
            % Kox [1/d] : Rate of reaction by oxidation of elemental mercury 
            % Fd [Ad] : Fraction of elemental mercury that is dissolved in water. 
            % K_Hg0 [1/d] : Rate of decay of elemental mercury by volatilization and oxidation

            % Decay rate of elemental mercury by volatilization and oxidation [1/d] [1/d]
            obj.K_Hg0    = (Fd.*(Vv./obj.H)) + Kox;
        end
        
        function Params_KHg2_V1(obj,Krx)
            % Decay rate of divalent mercury to elemental mercury by reduction 
            % Krx [1/d] : Decay rate by reduction 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury to elemental mercury by reduction

            obj.K_Hg2    = obj.ReachID.*0;
            obj.K_Hg2(:) = Krx;
        end
        
        function K__Hg0 = Params_K__Hg0(~,C_Hg2,Hg0u,K_Hg2,K_Hg0)                        
            % Elemental mercury decay rate combining reduction of divalent 
            % mercury and, volatilization and oxidation of elemental mercury 
            % C_Hg2 [mg/l] : Divalent mercury concentration 
            % Hg0u [mg/l] : Upstream elemental mercury concentration 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury to elemental mercury by reduction 
            % K_Hg0 [1/d] : Decay rate of elemental mercury by volatilization and oxidation 
            % K__Hg0[1/d] : Decay rate of elemental mercury combining reduction of divalent mercury and, volatilization and oxidation of elemental mercury
            
            % Controls for when concentrations are zero
            if Hg0u == 0
                Var = 1;
            elseif C_Hg2 == 0
                Var = 0;
            else
                Var = (C_Hg2./Hg0u);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Elemental mercury decay rate combining reduction of divalent 
            % mercury and volatilization and oxidation of elemental mercury [1/d]
            K__Hg0    = (K_Hg2*Var) - K_Hg0;
        end
        
        function WQS_Hg0(obj,varargin)  
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            addParameter(ip,'Krx',0.01,@isnumeric)   
            % Elemental mercury volatilization rate [m/d] Default value 
            % indicated in WASP [10] is taken
            addParameter(ip,'Vv',10,@isnumeric)            
            % Reaction rate for the oxidation of elemental mercury [1/d] 
            % Default value indicated in WASP [0.01] is taken
            addParameter(ip,'Kox',0.01,@isnumeric)
            % Fraction of elemental mercury found dissolved in water [Ad]. 
            % [0 log10 (L/Kg)] - Table 3. 
            % https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NERL&dirEntryId=135783
            Kd = 0;
            addParameter(ip,'Fd',1./(1 + ((10^Kd).*(obj.C_SS./1000000))),@ismatrix)
            % Elemental Mercury Load [mg/day]
            addParameter(ip, 'Load_Hg0',zeros(size(obj.ReachID)),@ismatrix)
            % Elemental Mercury Charge in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_Hg0',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})            
            Krx             = ip.Results.Krx;
            % Elemental mercury volatilization rate [m/d] Default value 
            Vv              = ip.Results.Vv;
            % Reaction rate for the oxidation of elemental mercury [1/d]
            Kox             = ip.Results.Kox;
            % Fraction of elemental mercury found dissolved in water [Ad].
            Fd              = ip.Results.Fd;
            % Elemental Mercury Load [mg/day]
            obj.Load_Hg0    = ip.Results.Load_Hg0;
            % Elemental Mercury Charge in Boundary Conditions [mg/day]
            obj.BC_Hg0      = ip.Results.BC_Hg0;
            
            % This function applies the elemental mercury model Initializing charges to 0
            obj.W_Hg0       = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_Hg0       = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_Hg0      = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of nitrification decay rate
            obj.Params_KHg0_V1(Vv,Kox,Fd)
            obj.Params_KHg2_V1(Krx)
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''Hg0'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''Hg0''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_Hg0 = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the elemental mercury model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % WQS - Divalent Mercury 
    % ---------------------------------------------------------------------
    properties        
        % Decay rate of divalent mercury combining, sedimentation, reduction 
        % and methylation of divalent mercury [1/d]
        K_Hg2
        % Decay rate of divalent mercury combining, oxidation of elemental 
        % mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
        K__Hg2
    end
    
    methods  
        function Params_KHg0_V2(obj,Kox)
            % Decay rate of elemental mercury to divalent mercury by oxidation 
            % Kox [1/d] : Decay rate by oxidation 
            % K_Hg0 [1/d] : Decay rate of elemental mercury to divalent mercury by oxidation

            obj.K_Hg0    = obj.ReachID.*0;
            obj.K_Hg0(:) = Kox;
        end
        
        function Params_KHg2_V2(obj,vs,Krx,Kme_d,Kme_p,Fp)
            % Decay rate of divalent mercury combining, sedimentation, 
            % reduction and methylation of divalent mercury [1/d] 
            % H [m] : Depth 
            % Vv [m/d] : Sedimentation velocity 
            % Krx [1/d] : Reduction rate 
            % Kme [1/d] : Methylation rate 
            % Fp [Ad] : Fraction of divalent mercury found as particulate in water. 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury combining, oxidation of elemental mercury, reduction and methylation of divalent mercury
            
            % Reaction rate for methylation of divalent mercury [1/d]
            Kme     = (Fp.*Kme_p) + ((1 - Fp).*Kme_d);
            
            % Decay rate of divalent mercury combining, sedimentation, reduction and methylation of divalent mercury [1/d]
            obj.K_Hg2    = (Fp.*(vs./obj.H)) + Krx + Kme;
        end
        
        function K__Hg2 = Params_K__Hg2(~,C_Hg0,Hg2u,K_Hg2,K_Hg0)                        
            % This function calculates the decay rate of divalent mercury 
            % combining, oxidation of elemental mercury, sedimentation, 
            % reduction and methylation of divalent mercury [1/d] 
            % C_Hg2 [mg/l] : Divalent mercury concentration 
            % Hg0u [mg/l] : Upstream elemental mercury concentration 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury combining, sedimentation, reduction and methylation of divalent mercury [1/d] 
            % K_Hg0 [1/d] : Decay rate of elemental mercury to divalent mercury by oxidation 
            % K__Hg2 [1/d] : Decay rate of divalent mercury combining, oxidation of elemental mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
            
            % Controls for when concentrations are zero
            if Hg2u == 0
                Var = 1;
            elseif C_Hg0 == 0
                Var = 0;
            else
                Var = (C_Hg0./Hg2u);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Decay rate of divalent mercury combining, oxidation of elemental mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
            K__Hg2    = (K_Hg0*Var) - K_Hg2;
        end
        
        function WQS_Hg2(obj,varargin)  
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Reaction rate for the oxidation of elemental mercury [1/d] 
            % Default value indicated in WASP [0.01] is taken
            addParameter(ip,'Kox',0.01,@isnumeric)
            % Divalent mercury sedimentation rate [m/d] 
            % Default value indicated in WASP [10] is taken
            addParameter(ip,'vs',0.6,@isnumeric)            
            % Reaction rate for reduction of divalent mercury [1/d] 
            % Default value indicated in WASP [0.01] is taken
            addParameter(ip,'Krx',0.01,@isnumeric)            
            % Dissolved HgII methylation rate [1/d] 
            % Default value indicated in WASP [0.001] is taken
            addParameter(ip,'Kme_d',0.001,@isnumeric)
            % Adsorbed HgII methylation rate [1/d] 
            % Default value indicated in WASP [0.01] is taken
            addParameter(ip,'Kme_p',0.01,@isnumeric) 
            % Fraction of divalent mercury found in particulate matter 
            % https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NERL&dirEntryId=135783 
            % [3.6 log10 (L/Kg)] - Table 3.
            Kd  = 3.6;
            addParameter(ip, 'Fp',((10^Kd).*(obj.C_SS./1000000))./(1 + ((10^Kd).*(obj.C_SS./1000000))),@ismatrix) 
            % Divalent Mercury Load [mg/day]
            addParameter(ip, 'Load_Hg2',zeros(size(obj.ReachID)),@ismatrix)
            % Divalent Mercury Load in Boundary Conditions [mg/day]
            addParameter(ip, 'BC_Hg2',zeros(size(obj.ReachID)),@ismatrix)
            % Checking input data
            parse(ip,varargin{:})
            % Reaction rate for oxidation of elemental mercury [1/d]
            Kox             = ip.Results.Kox;
            % Divalent mercury sedimentation rate [m/d]
            vs              = ip.Results.vs;
            % Reaction rate for reduction of divalent mercury [1/d]
            Krx             = ip.Results.Krx;
            % Dissolved HgII methylation rate [1/d]
            Kme_d           = ip.Results.Kme_d;
            % Adsorbed HgII methylation rate [1/d]
            Kme_p           = ip.Results.Kme_p;
            % Fraction of divalent mercury found in particulate matter 
            Fp              = ip.Results.Fp;
            % Divalent Mercury Load [mg/day]
            obj.Load_Hg2    = ip.Results.Load_Hg2;
            % Divalent Mercury Load in Boundary Conditions [mg/day]
            obj.BC_Hg2      = ip.Results.BC_Hg2;
            % This function applies the divalent mercury model Initializing charges to 0
            obj.W_Hg2       = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_Hg2       = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_Hg2      = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of nitrification decay rate
            obj.Params_KHg0_V2(Kox)
            obj.Params_KHg2_V2(vs,Krx,Kme_d,Kme_p,Fp)
            % Activate FunctionNetwork function mode
            obj.StatusFun       = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''Hg2'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''Hg2''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun   = false;
            obj.StatusWQS_Hg2 = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the divalent mercury model',toc);
            fprintf(WQSMessage)
        end
    end
    
    % --------------------------------------------------------------------- 
    % WQS - Methyl Mercury 
    % ---------------------------------------------------------------------
    properties        
        % Decay rate of divalent mercury combining, sedimentation, reduction 
        % and methylation of divalent mercury [1/d]
        K_MeHg
        % Decay rate of divalent mercury combining, oxidation of elemental 
        % mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
        K__MeHg
    end
    
    methods  
        function Params_KHg2_V3(obj,Kme_d,Kme_p,Fp)
            % Methylation rate of divalent mercury 
            % Kme [1/d] : Methylation rate 
            % K_Hg2 ​​[1/d] : Methylation rate of divalent mercury

            % Reaction rate for methylation of divalent mercury [1/d]
            Kme     = (Fp.*Kme_p) + ((1 - Fp).*Kme_d);
            
            obj.K_Hg2    = obj.ReachID.*0;
            obj.K_Hg2(:) = Kme;
        end
        
        function Params_KMeHg(obj,vs,Ku,Fp)
            % Rate of decay of methyl mercury by sedimentation and bioaccumulation [1/d] 
            % H [m] : Depth 
            % Vv [m/d] : Sedimentation velocity 
            % Krx [1/d] : Reduction rate 
            % Kme [1/d] : Methylation rate 
            % Fp [Ad] : Fraction of divalent mercury found in particulate form in water. 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury combining, oxidation of elemental mercury, reduction and methylation of divalent mercury
        
            obj.K_MeHg    = (Fp.*(vs./obj.H)) + Ku;
        end
        
        function K__MeHg = Params_K__MeHg(~,C_Hg2,MeHgu,K_Hg2,K_MeHg)                        
            % This function calculates the decay rate of divalent mercury 
            % combining, oxidation of elemental mercury, sedimentation, 
            % reduction and methylation of divalent mercury [1/d] 
            % C_Hg2 [mg/l] : Divalent mercury concentration 
            % Hg0u [mg/l] : Upstream elemental mercury concentration 
            % K_Hg2 ​​[1/d] : Decay rate of divalent mercury combining, sedimentation, reduction and methylation of divalent mercury [1/d] 
            % K_Hg0 [1/d] : Decay rate of elemental mercury to divalent mercury by oxidation 
            % K__Hg2 [1/d] : Decay rate of divalent mercury combining, oxidation of elemental mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
            
            % Controls for when concentrations are zero
            if MeHgu == 0 
                Var = 1;
            elseif C_Hg2 == 0
                Var = 0;
            else
                Var = (C_Hg2./MeHgu);
            end
            
            % Control so that no determinant matter is generated
            if Var > 1
                Var = 1;
            end
            
            % Decay rate of divalent mercury combining, oxidation of elemental mercury, sedimentation, reduction and methylation of divalent mercury [1/d]
            K__MeHg = (K_Hg2*Var) - K_MeHg;
        end
        
        function WQS_MeHg(obj,varargin) 
            % Start of stopwatch
            tic
            % Input parameters
            ip = inputParser;
            % Dissolved HgII methylation rate [1/d] 
            % Default value indicated in WASP [0.001] is taken
            addParameter(ip, 'Kme_d',0.001,@isnumeric)
            % Adsorbed HgII methylation rate [1/d] 
            % Default value indicated in WASP [0.01] is taken
            addParameter(ip, 'Kme_p',0.01,@isnumeric)            
            % Fraction of divalent mercury found in particulate matter The 
            % default value indicated in WASP is taken [0.01] 
            % https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NERL&dirEntryId=135783 
            % [2.7 log10 (L/Kg)] - Table 3.
            Kd  = 2.7;
            addParameter(ip, 'Fp',((10^Kd).*(obj.C_SS./1000000))./(1 + ((10^Kd).*(obj.C_SS./1000000))),@ismatrix)
            % Methyl mercury sedimentation rate [m/d] 
            % Default value indicated in WASP [10] is taken
            addParameter(ip, 'vs',0.5,@isnumeric)                     
            % Bioaccumulation rate of methyl mercury [1/d] 
            % Left at zero. That is, this process is not considered in the model.
            addParameter(ip, 'Ku',0,@isnumeric)          
            % Methyl Mercury Load [mg/day]
            addParameter(ip, 'Load_MeHg',zeros(size(obj.ReachID)),@ismatrix) 
            % Methyl Mercury Load in Boundary Conditions [mg]
            addParameter(ip, 'BC_MeHg',zeros(size(obj.ReachID)),@ismatrix) 

            % Checking input data
            parse(ip,varargin{:})
            % Dissolved HgII methylation rate [1/d]
            Kme_d           = ip.Results.Kme_d;
            % Adsorbed HgII methylation rate [1/d]
            Kme_p           = ip.Results.Kme_p;
            % Fraction of divalent mercury found in particulate matter 
            Fp              = ip.Results.Fp;
            % Methyl mercury sedimentation rate [m/d] 
            vs              = ip.Results.vs;
            % Bioaccumulation rate of methyl mercury [1/d] 
            Ku              = ip.Results.Ku;
            % Methyl Mercury Load [mg/day]
            obj.Load_MeHg   = ip.Results.Load_MeHg; 
            % Methyl Mercury Load in Boundary Conditions [mg]
            obj.BC_MeHg     = ip.Results.BC_MeHg;
            % This function applies the divalent mercury model Initializing charges to 0
            obj.W_MeHg      = zeros(size(obj.ReachID));
            % Initializing concentrations to 0
            obj.C_MeHg      = zeros(size(obj.ReachID));
            % Initializing assimilation factors to 0
            obj.AF_MeHg     = zeros(size(obj.ReachID));
            % Initialization of flow accumulator
            obj.Qaccum      = zeros(size(obj.ReachID));
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Calculation of mercury decay rate
            obj.Params_KMeHg(vs,Ku,Fp)
            obj.Params_KHg2_V3(Kme_d,Kme_p,Fp)
            % Activate FunctionNetwork function mode
            obj.StatusFun   = true;
            % Functions for accumulation in 1 and 2 sections
            obj.FunNetwork_1= 'obj.FunctionNetwork_CoDependent(Npre, Posi(i),''Model'',''MeHg'');';
            % Function for header section (when Posi = [])
            obj.FunNetwork_2= 'if isempty(Posi), obj.FunctionNetwork_CoDependent(Npre, Posi,''Model'',''MeHg''), end;';
            % Apply cumulative scheme
            obj.AnalysisNetwork_Obj;
            % Disable FunctionNetwork function mode
            obj.StatusFun       = false;
            obj.StatusWQS_MeHg = true;
            % Assignment of discharge status
            obj.Update_VerStatus;
            % Show results
            WQSMessage = sprintf('%60s | Time: %.4f sec\n','Execution of the methyl mercury model',toc);
            fprintf(WQSMessage)
        end
    end
        
    % -------------------------------------------------- ------------------- 
    % Functions to estimate reaction rates or physical parameters of the stream 
    % -------------------- -------------------------------------------------
    methods             
        % Reaction rate
        function k = Kd_ArrheniusModel(obj,varargin)
            % This function calculates the reaction rates of different determinants 
            % T [°C] : Temperature 
            % Theta [Ad] : Temperature correction factor 
            % k [1/d] : Interest rate calculated for a reference temperature of 20°C
            
            errorMsg = 'The argument must be matrix or numeric.';
            validationFcn = @(x) assert(isnumeric(x) || ismatrix(x),errorMsg);

            ip = inputParser;
            % Input parameter name
            addParameter(ip, 'NameParam','X',@ischar) 
            addParameter(ip, 'k',0.02,validationFcn) 
            % Check data
            parse(ip,varargin{:})
            NameParam   = ip.Results.NameParam;
            k           = ip.Results.k;

            % Internal parameters
            switch NameParam
                case 'OM'
                    % Temperature correction factor
                    Teta    = 1.047;
                case 'NO'
                    % Temperature correction factor
                    Teta    = 1.047;
                case 'NH4'
                    % Temperature correction factor
                    Teta    = 1.047;
                case 'NO3'
                    % Temperature correction factor
                    Teta    = 1.0698;
                case 'X'
                    % Temperature correction factor
                    Teta    = 1.07;
                case 'DO'
                    % Temperature correction factor
                    Teta    = 1.024;
            end

            % Reaction rates [1/d]
            k = k.*(Teta.^(obj.T - 20));   
        end
        
        
        % Update status of discharges
        function Update_VerStatus(obj)
            % This function determines that UHM has the presence of discharges in order to consider them in the cumulative scheme.
            % Tmp = obj.Load_X + obj.Load_PO + obj.Load_PI + obj.Load_NH4 + ...
            %       obj.Load_NO + obj.Load_OM + obj.Load_DOD + obj.Load_T + ... 
            %       obj.Load_Hg0 + obj.Load_Hg2 + obj.Load_MeHg + ...
            %       obj.Load_SS;
            Tmp = obj.ReachID;
            % Boolean Vector
            obj.Load_Status = (Tmp~=0);
        end
                
        % Dispersive fraction
        function DF_Gonzalez(obj)
            % This function assigns the dispersive fraction according to 
            % the type of river (Mountain or Plain). In previously carried 
            % out works (González, 2008) it has been found that the dispersive 
            % fraction (DF) for mountain rivers has a global magnitude of 0.27 
            % while for plain rivers it is 0.40. DF [Ad] : Dispersive fraction
            
            % Dispersive fraction for mountain rivers
            DF_MR   = 0.27;
            % Dispersive fraction for plain rivers
            DF_PR   = 0.40;
            % Assign dispersive fraction values
            obj.DF = obj.ReachID*0;
            obj.DF(obj.RiverType)   = DF_MR;
            obj.DF(~obj.RiverType)  = DF_PR;              
        end
                
        % Effective delay coefficient
        function Beta_Gonzalez(obj)
            % This function assigns the effective delay coefficient according 
            % to the type of river (Mountain or Plain). In previously carried 
            % out works (Gonzalez, 2008) it has been found that the effective 
            % delay coefficient for mountain rivers has a global magnitude of 
            % 1.10 while for plain rivers it is 2.0. 
            % Beta [Ad] : effective delay coefficient
            
            % Effective retardation coefficient for mountain rivers
            Beta_MR   = 1.1;
            % Effective retardation coefficient for plain rivers
            Beta_PR   = 2.0;
            % Assign values
            obj.Beta = obj.ReachID*0;
            obj.Beta(obj.RiverType)   = Beta_MR;
            obj.Beta(~obj.RiverType)  = Beta_PR;            
        end
                
        % Solute speed
        function Vs_Lees(obj)
            % This function estimates the solute velocity according to the 
            % relationship proposed by Lees et al. (2000). Lees et al. (2000) 
            % presented a relationship between flow velocity and solute velocity 
            % based on the moment matching technique between the ADZ-QUASAR 
            % and TS models. From this work it is concluded that the velocity 
            % relationship, both for plain rivers and for mountain rivers, 
            % is expressed by the effective delay coefficient ? 
            % (corresponding to the relationship between the cross-sectional 
            % areas of the main channel and the dead zones) (Camacho, 2000) 
            % (Rojas, 2011) 
            % v [m/s] : Flow velocity 
            % Vs [m/s] : Solute velocity
            
            % Solute velocity [m/s]
            obj.Vs = obj.v./(1 + obj.Beta);         
        end                
            
        % Variable Accumulation
        function VariableAccumulation(obj)  
            % This function accumulates the area of ​​the downstream sections 
            % A     [m2]    : Area of ​​the UHM
            % Qr    [m3/s]  : River discharge
            % Qwwd  [m3/s]  : Wastewater discharge
            
            % Only accumulate once
            if ~obj.StatusAccArea
                tic
                obj.AccumVar = [obj.A obj.Qwwd];
                obj.AnalysisNetwork_Obj;
                obj.BC_Status = (obj.A == obj.AccumVar(:,1));
                obj.A = obj.AccumVar(:,1);
                obj.Q = obj.Qr + obj.AccumVar(:,2);
                
                % State that the area is already accumulated
                obj.StatusAccArea = true;
                
                % Show results
                disp(['Acumulation Ok | Time: ',num2str(toc,'%.4f'),' Seg']);
            end
        end
                
        % Oxygen saturation
        function Cal_OS(obj)
            % This function estimates the saturation oxygen concentration 
            % of the stream according to the proposal of Zison et al. (1978) 
            % which considers correction for elevation 
            % T [°C] : Temperature 
            % Z [m.a.s.l] : Stream elevation 
            % C_OS [mg/l] : Dissolved oxygen concentration
            
            % Coversion Factor: Celsius to Kelvin
            ConFac = 273.15;
            % Only accumulate 1 time
            obj.C_OS = (1-(0.0001148.*obj.Z)).*...
                exp(-139.34411 + ...
                (1.575701E5./ ((obj.T+ConFac).^1)) -...
                (6.642308E7./ ((obj.T+ConFac).^2)) +...
                (1.243800E10./((obj.T+ConFac).^3))-...
                (8.621949E11./((obj.T+ConFac).^4)));
        end
        
        function Cal_Turbidity(obj)
            % This function estimates turbidity using an empirical relationship 
            % with data collected by Cornare in measurement campaigns carried 
            % out with variable periodicity between 2008 and 2018 in the 
            % Fe-Pantanillo basin. 
            % SS [mg/l] : Total suspended solids 
            % Turb [UNF] : Turbidity nephelometric turbidity units
            
            % Coversion Factor: [mg/l] to [Kg/m3]
            ConFac    = 1/1000;
            obj.Turb  = 304.90*((obj.C_SS.*ConFac).^0.81);
        end        
    end
end