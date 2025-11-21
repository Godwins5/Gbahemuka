% Main class for hydrological catchment analysis
classdef Catchment < handle
    properties
        Name
        DEM                          % Digital Elevation Model
        FlowDirection
        FlowAccumulation
        StreamNetwork
        PourPoints
        Subcatchments
        LandUse
        SoilType
        Area                         % km²
        MeanSlope
        TimeOfConcentration          % hours
    end
    
    methods
        function obj = Catchment(name, dem)
            obj.Name = name;
            obj.DEM = dem;
            obj.Subcatchments = Subcatchment.empty;
        end
        
        function delineateCatchment(obj, pourPoint)
            fprintf('Delineating catchment: %s\n', obj.Name);
            
            % Fill sinks in DEM
            filledDEM = obj.fillSinks();
            
            % Calculate flow direction
            obj.calculateFlowDirection(filledDEM);
            
            % Calculate flow accumulation
            obj.calculateFlowAccumulation();
            
            % Extract stream network
            obj.extractStreamNetwork(1000); % Threshold of 1000 cells
            
            % Delineate catchment boundary
            obj.delineateWatershed(pourPoint);
            
            % Calculate catchment properties
            obj.calculateCatchmentProperties();
        end
        
        function filledDEM = fillSinks(obj)
            % Simple sink filling algorithm
            [rows, cols] = size(obj.DEM);
            filledDEM = obj.DEM;
            
            for i = 2:rows-1
                for j = 2:cols-1
                    neighbors = filledDEM(i-1:i+1, j-1:j+1);
                    minNeighbor = min(neighbors(:));
                    if filledDEM(i,j) < minNeighbor
                        filledDEM(i,j) = minNeighbor;
                    end
                end
            end
        end
        
        function calculateFlowDirection(obj, filledDEM)
            % D8 flow direction algorithm
            [rows, cols] = size(filledDEM);
            obj.FlowDirection = zeros(rows, cols);
            
            % Directions: 1=E, 2=SE, 3=S, 4=SW, 5=W, 6=NW, 7=N, 8=NE
            for i = 2:rows-1
                for j = 2:cols-1
                    if filledDEM(i,j) > 0 % Skip NoData cells
                        % Get elevation differences to neighbors
                        elevations = filledDEM(i-1:i+1, j-1:j+1);
                        center = elevations(2,2);
                        elevations(2,2) = inf; % Ignore center cell
                        
                        % Find steepest descent
                        [maxDrop, idx] = max(center - elevations(:));
                        
                        if maxDrop > 0
                            obj.FlowDirection(i,j) = idx;
                        end
                    end
                end
            end
        end
        
        function calculateFlowAccumulation(obj)
            [rows, cols] = size(obj.FlowDirection);
            obj.FlowAccumulation = ones(rows, cols);
            
            % Calculate flow accumulation based on flow direction
            for i = 1:rows
                for j = 1:cols
                    if obj.FlowDirection(i,j) > 0
                        % Accumulate flow from upstream cells
                        obj.accumulateFlow(i, j);
                    end
                end
            end
        end
        
        function accumulateFlow(obj, i, j)
            % Recursive function to accumulate flow
            [rows, cols] = size(obj.FlowDirection);
            
            % Define neighbor offsets for each direction
            offsets = [0,1; 1,1; 1,0; 1,-1; 0,-1; -1,-1; -1,0; -1,1];
            
            for dir = 1:8
                ni = i + offsets(dir,1);
                nj = j + offsets(dir,2);
                
                if ni >= 1 && ni <= rows && nj >= 1 && nj <= cols
                    if obj.FlowDirection(ni, nj) == dir
                        obj.FlowAccumulation(i,j) = obj.FlowAccumulation(i,j) + ...
                            obj.accumulateFlow(ni, nj);
                    end
                end
            end
        end
        
        function extractStreamNetwork(obj, threshold)
            obj.StreamNetwork = obj.FlowAccumulation >= threshold;
            fprintf('Stream network extracted with threshold: %d cells\n', threshold);
        end
        
        function delineateWatershed(obj, pourPoint)
            % Delineate watershed for given pour point
            [rows, cols] = size(obj.FlowDirection);
            watershed = false(rows, cols);
            
            queue = pourPoint;
            watershed(pourPoint(1), pourPoint(2)) = true;
            
            offsets = [0,1; 1,1; 1,0; 1,-1; 0,-1; -1,-1; -1,0; -1,1];
            
            while ~isempty(queue)
                current = queue(1,:);
                queue(1,:) = [];
                
                for dir = 1:8
                    ni = current(1) + offsets(dir,1);
                    nj = current(2) + offsets(dir,2);
                    
                    if ni >= 1 && ni <= rows && nj >= 1 && nj <= cols
                        if ~watershed(ni, nj) && obj.FlowDirection(ni, nj) == dir
                            watershed(ni, nj) = true;
                            queue = [queue; ni, nj];
                        end
                    end
                end
            end
            
            obj.PourPoints = watershed;
        end
        
        function calculateCatchmentProperties(obj)
            % Calculate basic catchment properties
            cellArea = 30 * 30 / 1e6; % Assuming 30m resolution, convert to km²
            obj.Area = sum(obj.PourPoints(:)) * cellArea;
            
            % Calculate mean slope
            [fx, fy] = gradient(obj.DEM);
            slope = sqrt(fx.^2 + fy.^2);
            obj.MeanSlope = mean(slope(obj.PourPoints), 'all');
            
            fprintf('Catchment Area: %.2f km²\n', obj.Area);
            fprintf('Mean Slope: %.2f degrees\n', rad2deg(atan(obj.MeanSlope)));
        end
        
        function plotCatchment(obj)
            figure;
            subplot(2,2,1);
            imagesc(obj.DEM);
            title('Digital Elevation Model');
            colorbar;
            axis equal;
            
            subplot(2,2,2);
            imagesc(obj.FlowAccumulation);
            title('Flow Accumulation');
            colorbar;
            axis equal;
            
            subplot(2,2,3);
            imagesc(obj.StreamNetwork);
            title('Stream Network');
            axis equal;
            
            subplot(2,2,4);
            imagesc(obj.PourPoints);
            title('Catchment Boundary');
            axis equal;
            
            colormap(jet);
        end
    end
end

% Subcatchment class for nested watersheds
classdef Subcatchment < handle
    properties
        ID
        Area
        Slope
        LandUse
        SoilType
        RunoffCoefficient
        Outlet
        UpstreamCatchments
    end
    
    methods
        function obj = Subcatchment(id, area, slope)
            obj.ID = id;
            obj.Area = area;
            obj.Slope = slope;
            obj.UpstreamCatchments = Subcatchment.empty;
        end
        
        function calculateRunoffCoefficient(obj)
            % Calculate runoff coefficient based on land use and soil type
            runoffCoefficients = containers.Map(...
                {'forest', 'agriculture', 'urban', 'pasture'}, ...
                [0.1, 0.3, 0.7, 0.2]);
            
            if isKey(runoffCoefficients, obj.LandUse)
                baseCoefficient = runoffCoefficients(obj.LandUse);
            else
                baseCoefficient = 0.3;
            end
            
            % Adjust for slope
            slopeAdjustment = min(1.0, 0.5 + obj.Slope * 10);
            obj.RunoffCoefficient = baseCoefficient * slopeAdjustment;
        end
    end
end

% Runoff Analysis Class
classdef RunoffAnalyzer < handle
    properties
        Catchment
        RainfallData
        TimeSeries
    end
    
    methods
        function obj = RunoffAnalyzer(catchment)
            obj.Catchment = catchment;
        end
        
        function calculateSCSCurveNumber(obj, landUse, soilType)
            % SCS Curve Number method
            cnTable = containers.Map();
            
            % Define Curve Numbers for different land use/soil combinations
            cnTable('forest_A') = 30;
            cnTable('forest_B') = 55;
            cnTable('forest_C') = 70;
            cnTable('forest_D') = 77;
            
            cnTable('agriculture_A') = 67;
            cnTable('agriculture_B') = 78;
            cnTable('agriculture_C') = 85;
            cnTable('agriculture_D') = 89;
            
            cnTable('urban_A') = 77;
            cnTable('urban_B') = 85;
            cnTable('urban_C') = 90;
            cnTable('urban_D') = 92;
            
            key = [landUse '_' soilType];
            if isKey(cnTable, key)
                curveNumber = cnTable(key);
            else
                curveNumber = 70; % Default
            end
            
            fprintf('Curve Number for %s: %d\n', key, curveNumber);
            return curveNumber;
        end
        
        function runoff = calculateRunoffSCS(obj, rainfall, curveNumber)
            % Calculate runoff using SCS method
            S = (1000 / curveNumber) - 10; % Potential maximum retention
            if rainfall > 0.2 * S
                runoff = (rainfall - 0.2 * S)^2 / (rainfall + 0.8 * S);
            else
                runoff = 0;
            end
        end
        
        function hydrograph = calculateUnitHydrograph(obj, duration, catchmentArea)
            % Calculate SCS Unit Hydrograph
            tp = 0.6 * obj.Catchment.TimeOfConcentration + duration/2; % Time to peak
            qp = 2.08 * catchmentArea / tp; % Peak discharge
            
            time = 0:0.1:tp*3;
            hydrograph = zeros(size(time));
            
            for i = 1:length(time)
                t = time(i);
                if t <= tp
                    hydrograph(i) = qp * (t / tp);
                else
                    hydrograph(i) = qp * exp(1 - (t / tp));
                end
            end
            
            figure;
            plot(time, hydrograph);
            title('SCS Unit Hydrograph');
            xlabel('Time (hours)');
            ylabel('Discharge (m³/s)');
            grid on;
        end
    end
end

% Hydrological Model Class
classdef HydrologicalModel < handle
    properties
        Name
        Catchments
        Rainfall
        TimeStep
        Results
    end
    
    methods
        function obj = HydrologicalModel(name)
            obj.Name = name;
            obj.Catchments = Catchment.empty;
        end
        
        function addCatchment(obj, catchment)
            obj.Catchments(end+1) = catchment;
        end
        
        function simulateRunoff(obj, rainfallDuration, rainfallIntensity)
            fprintf('Running hydrological simulation: %s\n', obj.Name);
            
            for i = 1:length(obj.Catchments)
                catchment = obj.Catchments(i);
                fprintf('\nAnalyzing catchment: %s\n', catchment.Name);
                
                % Calculate runoff for different scenarios
                analyzer = RunoffAnalyzer(catchment);
                
                % Example rainfall event
                rainfall = rainfallIntensity * rainfallDuration; % mm
                
                % Calculate runoff using SCS method
                cn = analyzer.calculateSCSCurveNumber('agriculture', 'B');
                runoffDepth = analyzer.calculateRunoffSCS(rainfall, cn);
                runoffVolume = runoffDepth * catchment.Area * 1000; % m³
                
                fprintf('Rainfall: %.1f mm\n', rainfall);
                fprintf('Runoff Depth: %.1f mm\n', runoffDepth);
                fprintf('Runoff Volume: %.0f m³\n', runoffVolume);
                
                % Generate unit hydrograph
                analyzer.calculateUnitHydrograph(rainfallDuration, catchment.Area);
            end
        end
        
        function plotResults(obj)
            figure;
            for i = 1:length(obj.Catchments)
                catchment = obj.Catchments(i);
                
                subplot(2,2,i);
                imagesc(catchment.PourPoints);
                title(sprintf('Catchment: %s (%.1f km²)', catchment.Name, catchment.Area));
                axis equal;
            end
        end
    end
end

% Terrain Analysis Class
classdef TerrainAnalyzer
    methods (Static)
        function slopeMap = calculateSlope(dem, cellSize)
            % Calculate slope from DEM
            [fx, fy] = gradient(dem, cellSize);
            slopeMap = atan(sqrt(fx.^2 + fy.^2));
        end
        
        function aspectMap = calculateAspect(dem, cellSize)
            % Calculate aspect from DEM
            [fx, fy] = gradient(dem, cellSize);
            aspectMap = atan2(-fy, fx);
            aspectMap(aspectMap < 0) = aspectMap(aspectMap < 0) + 2 * pi;
        end
        
        function wetnessIndex = calculateTopographicWetnessIndex(dem, flowAccumulation, cellSize)
            % Calculate topographic wetness index
            slope = TerrainAnalyzer.calculateSlope(dem, cellSize);
            slope(slope == 0) = 0.001; % Avoid division by zero
            
            specificCatchmentArea = flowAccumulation * cellSize^2;
            wetnessIndex = log(specificCatchmentArea ./ slope);
        end
    end
end

% Example usage and demonstration
function runHydrologicalExample()
    % Create a synthetic DEM for demonstration
    [X, Y] = meshgrid(1:100, 1:100);
    dem = 500 - sqrt((X-50).^2 + (Y-30).^2) + randn(100,100)*5;
    dem = max(dem, 0); % Ensure positive elevations
    
    % Create catchment
    mainCatchment = Catchment('Main River Basin', dem);
    
    % Delineate catchment with pour point at [80, 50]
    mainCatchment.delineateCatchment([80, 50]);
    
    % Calculate time of concentration (simplified)
    mainCatchment.TimeOfConcentration = 1.5; % hours
    
    % Create hydrological model
    model = HydrologicalModel('Regional Hydrological Model');
    model.addCatchment(mainCatchment);
    
    % Run simulation for 24-hour storm with 10 mm/hr intensity
    model.simulateRunoff(24, 10);
    
    % Plot results
    mainCatchment.plotCatchment();
    model.plotResults();
    
    % Demonstrate terrain analysis
    slope = TerrainAnalyzer.calculateSlope(dem, 30);
    aspect = TerrainAnalyzer.calculateAspect(dem, 30);
    wetnessIndex = TerrainAnalyzer.calculateTopographicWetnessIndex(...
        dem, mainCatchment.FlowAccumulation, 30);
    
    figure;
    subplot(1,3,1);
    imagesc(slope);
    title('Slope Map');
    colorbar;
    
    subplot(1,3,2);
    imagesc(aspect);
    title('Aspect Map');
    colorbar;
    
    subplot(1,3,3);
    imagesc(wetnessIndex);
    title('Topographic Wetness Index');
    colorbar;
    
    colormap(jet);
end

% Run the example
fprintf('=== HYDROLOGICAL AND CATCHMENT ANALYSIS ===\n',);
runHydrologicalExample();
