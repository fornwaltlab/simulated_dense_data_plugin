classdef SimulatedDENSEDataPlugin < plugins.DENSEanalysisPlugin
    % SimulatedDENSEDataPlugin - A DENSEanalysis plugin
    %
    %   Generates a 2D cardiac displacement field with twist
    %
    % Copyright (c) 2016, Jonathan Suever

    methods
        function validate(self, data)
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

            % Assert that image data base been loaded
            assert(~isempty(data.seq), ...
                'You must load imaging data into DENSEanalysis first.')
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
        end
    end
end
