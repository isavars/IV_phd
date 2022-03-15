classdef nexus < handle
    % nexus data class. Let's you load in neuronexus data into a class
    % object for analysis or data inspection
    %
    % Class constructor: scanpix.nexus
    % Construct class obj
    %
    % Usage:
    %       obj = scanpix.nexus;
    %       obj = scanpix.nexus(prmsMode);
    %       obj = scanpix.nexus(prmsMode, uiFlag);
    %
    % Inputs:
    %       prmsMode: 'default' - uses default parameter (default)
    %                 'ui'      - opens UI dialogue
    %                 'file'    - load params from file
    %
    %         uiFlag: true (default) or false
    %                 - if false will skip UI dialogue for data selection
    %                   (e.g. when you use constructor programmatically)
    %
    %  IV
    %% PROPERTIES
    properties % params %
        % Map object
        params                containers.Map
        chanMap               struct
    end
    
end 