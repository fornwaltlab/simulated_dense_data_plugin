classdef SimulatedDENSEDataPlugin < plugins.DENSEanalysisPlugin
    % SimulatedDENSEDataPlugin - A DENSEanalysis plugin
    %
    %   Generates a 2D cardiac displacement field with twist
    %
    % Copyright (c) 2016, Cardiac Imaging Technology Lab

    methods
        function validate(varargin)
            % validate - Check if the plugin can run.
            %
            %   Performs validation to ensure that the state of the program
            %   is correct to be able to run the plugin.
            %
            % USAGE:
            %   SimulatedDENSEDataPlugin.validate(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.

            % This plugin will always be able to be run
        end

        function run(self, data)
            % run - Method executed when user selects the plugin
            %
            % USAGE:
            %   SimulatedDENSEDataPlugin.run(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.

            params = getfield(self.Config, 'Parameters', struct());

            d = plugins.simulation.DNSCreator('Data', data, struct(params));
            L = addlistener(d, 'StateChanged', @(s,e)self.updateConfig(e));
            addlistener(d, 'ObjectBeingDestroyed', @(s,e)self.updateConfig(d))
        end

        function updateConfig(self, src)
            setfield(self.Config, 'Parameters', struct(src));
        end
    end
end
