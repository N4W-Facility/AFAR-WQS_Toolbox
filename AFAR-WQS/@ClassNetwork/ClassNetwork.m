classdef (HandleCompatible) ClassNetwork < handle
% -------------------------------------------------------------------------
% Matlab-R2023b
% -------------------------------------------------------------------------
%                              INFORMATION 
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
%--------------------------------------------------------------------------
%                               DESCRIPTION 
%--------------------------------------------------------------------------
% This object represent a topological network 
%
%--------------------------------------------------------------------------
%                               INPUT DATA 
%--------------------------------------------------------------------------
%   n = it is the section number that haves the topological network
%
%   ID          [n,1] = ID of the River Sections (Ad)
%   FromNode    [n,1] = Initial Node of the River Sections (Ad)
%   ToNode      [n,1] = End Node of the River Sections (Ad)
%   RiverMouth  [1,1] = ID of the River Sections corresponding to the River 
%                       Mouth (Ad)                    
%
%--------------------------------------------------------------------------
%                              OUTPUT DATA 
%--------------------------------------------------------------------------
% ClassNetwork [Object] = This object contain a topological network 

    %% Properties
    % Properties of the topological network
    properties        
        % ID of the River Sections (Ad)
        ReachID(:,1) double
        % Initial Node of the River Sections (Ad)
        FromNode(:,1) double
        % End Node of the River Sections (Ad)
        ToNode(:,1) double
        % ReachID of the River Sections corresponding to the River Mouth (Ad)
        RiverMouth(:,1) double
        % Barrier
        Barrier(:,1) double        
    end
    
    % Variables for object-oriented network analysis
    properties  
        % Variables to Propagate
        ProVar double
        % Variables to Accumulate
        AccumVar double
        % Variables to Accumulate with Clipping
        AccumClipVar double
        % Loss Rate (%)
        LossRate double
        % Variables to Accumulate with Losses
        AccumLossVar double
        % Variables to Accumulate with Losses and Clipping
        AccumClipLossVar double
        % Status function
        StatusFun(1,1) logical = false
        % Function        
        FunNetwork_1 char = ';'
        % Function        
        FunNetwork_2 char = ';'
    end
    
    methods      
        %% Analysis Network
        [FuncNetwork, varargout]    = AnalysisNetwork_Obj(obj, varargin)        
        [FuncNetwork,PoNet]         = FunctionalBranch(obj,ReachID_RM,ReachID_RM_i,CurrID,PoNet)
    end

end