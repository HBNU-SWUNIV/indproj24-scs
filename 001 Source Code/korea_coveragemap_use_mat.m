%%%Coverage Maps for Satellite Constellation
%% Create and Visualize Scenario
startTime = datetime(2023,2,21,18,0,0);
stopTime = startTime + hours(1);
sampleTime = 60; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

%% Add Satellite Constellationnum
% 위성의 개수 설정
numOrbits=11; 
numSatellitesPerOrbitalPlane = 11;
numSatellites = numOrbits*numSatellitesPerOrbitalPlane; % Walker Star 모델에서는 총 위성의 개수를 설정

% RAAN (Right Ascension of Ascending Node) 설정 - 모든 위성이 동일한 궤도 공유 => 동일하게 설정
RAAN = zeros(1, numSatellites); % 동일한 궤도 평면

% True Anomaly 설정 - 균등한 각도 간격으로 분포
trueanomaly = zeros(1,numSatellites);

satIndex=1;
for i=1:numOrbits
    for j=1:numSatellitesPerOrbitalPlane
        RAAN(satIndex)=360*(i-1)/numOrbits;
        trueanomaly(satIndex)=360*(j-1)/numSatellitesPerOrbitalPlane;
        satIndex=satIndex+1;
    end
end

% 궤도 요소 설정
semimajoraxis = repmat((6371 + 780)*1e3, size(RAAN)); % 반장축 (고도 780km)
inclination = repmat(86.4, size(RAAN)); % 경사각
eccentricity = zeros(size(RAAN)); % 이심률 (원형 궤도)
argofperiapsis = zeros(size(RAAN)); % 근지점 인수

% 위성 생성
sats = satellite(sc, semimajoraxis, eccentricity, inclination, RAAN, argofperiapsis, trueanomaly, Name="Star " + string(1:numSatellites)');


%% Add Grid of Ground Stations Covering Korea
% 대한민국의 좌표계를 설정
latlim = [33, 39];  % 대한민국의 위도 범위
lonlim = [124, 131];  % 대한민국의 경도 범위

proj = projcrs(3857);
spacingInLatLon = 1; % degrees
%x-y 지도 좌표 계산
spacingInXY = deg2km(spacingInLatLon)*1000; % meters
[xlim,ylim] = projfwd(proj,latlim,lonlim);

R = maprefpostings(xlim,ylim,spacingInXY,spacingInXY);
[X,Y] = worldGrid(R);
[gridlat,gridlon] = projinv(proj,X,Y);
%korea추출
% 새로운 Shapefile의 내용을 확인
% Shapefile 읽기 시 좌표계를 지리적 좌표계로 지정
landareas = readgeotable('combined_landareas.shp', 'CoordinateSystemType', 'geographic');

korea = landareas(string(landareas.Name) == "South Korea", :);

% MAT 파일에서 한국 지리 정보 불러오기
load('korea_coverage_map_all.mat','gslat', 'gslon', 'koreab','gridpts','inregion');

gs = groundStation(sc,gslat,gslon);

%% Add Transmitters and Receivers
fq = 1625e6; % Hz
txpower = 20; % dBW
antennaType = "Gaussian";
halfBeamWidth = 62.7; % degrees

if antennaType == "Gaussian"
    lambda = physconst('lightspeed')/fq; % meters
    dishD = (70*lambda)/(2*halfBeamWidth); % meters
    tx = transmitter(sats, ...
        Frequency=fq, ...
        Power=txpower); 
    gaussianAntenna(tx,DishDiameter=dishD);
end

if antennaType == "Custom 48-Beam"
    antenna = helperCustom48BeamAntenna(fq);
    tx = transmitter(sats, ...
        Frequency=fq, ...
        MountingAngles=[0,-90,0], ... % [yaw, pitch, roll] with -90 using Phased Array System Toolbox convention
        Power=txpower, ...
        Antenna=antenna);  
end

isotropic = arrayConfig(Size=[1 1]);
rx = receiver(gs,Antenna=isotropic);

pattern(tx,Size=500000);

%% Compute Raster Coverage Map Data
delete(viewer)
maxsigstrength = satcoverage(gridpts,sc,startTime,inregion,halfBeamWidth);

%% Visualize Coverage on an axesm-Based Map
minpowerlevel = -120; % dBm
maxpowerlevel = -106; % dBm

figure
worldmap(latlim,lonlim)
mlabel north

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoshow(gridlat,gridlon,maxsigstrength,DisplayType="contour",Fill="on")
geoshow(korea,FaceColor="none")

cBar = contourcbar;
title(cBar,"dBm");
title("Signal Strength at " + string(startTime) + " UTC")

sydlatlon = [-32.13 151.21]; % Sydney
mellatlot = [-36.19 144.96]; % Melbourne
brislatlon = [-26.53 153.03]; % Brisbane

textm(sydlatlon(1),sydlatlon(2),"Sydney")
textm(mellatlot(1),mellatlot(2),"Melbourne")
textm(brislatlon(1),brislatlon(2),"Brisbane")

%% Compute Coverage Map Contours
levels = linspace(minpowerlevel,maxpowerlevel,8);
GT = contourDataGrid(gridlat,gridlon,maxsigstrength,levels,proj);
GT = sortrows(GT,"Power (dBm)");
disp(GT)

%% Visualize Coverage on a Map Axes
figure
newmap(proj)
hold on

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoplot(GT,ColorVariable="Power (dBm)",EdgeColor="none")
geoplot(korea,FaceColor="none")

cBar = colorbar;
title(cBar,"dBm");
title("Signal Strength at " + string(startTime) + " UTC")

text(sydlatlon(1),sydlatlon(2),"Sydney",HorizontalAlignment="center")
text(mellatlot(1),mellatlot(2),"Melbourne",HorizontalAlignment="center")
text(brislatlon(1),brislatlon(2),"Brisbane",HorizontalAlignment="center")

%% Compute and Visualize Coverage for a Different Time of Interest

secondTOI = startTime + minutes(2); % 2 minutes after the start of the scenario
maxsigstrength = satcoverage(gridpts,sc,secondTOI,inregion,halfBeamWidth);

GT2 = contourDataGrid(gridlat,gridlon,maxsigstrength,levels,proj);
GT2 = sortrows(GT2,"Power (dBm)");

figure
newmap(proj)
hold on

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoplot(GT2,ColorVariable="Power (dBm)",EdgeColor="none")
geoplot(korea,FaceColor="none")

cBar = colorbar;
title(cBar,"dBm");
title("Signal Strength at " + string(secondTOI) + " UTC")

text(sydlatlon(1),sydlatlon(2),"Sydney",HorizontalAlignment="center")
text(mellatlot(1),mellatlot(2),"Melbourne",HorizontalAlignment="center")
text(brislatlon(1),brislatlon(2),"Brisbane",HorizontalAlignment="center")

%% Compute Coverage Levels for Cities
covlevels1 = [coveragelevel(sydlatlon(1),sydlatlon(2),GT); ...
    coveragelevel(mellatlot(1),mellatlot(2),GT); ...
    coveragelevel(brislatlon(1),brislatlon(2),GT)];
covlevels2 = [coveragelevel(sydlatlon(1),sydlatlon(2),GT2); ...
    coveragelevel(mellatlot(1),mellatlot(2),GT2); ...
    coveragelevel(brislatlon(1),brislatlon(2),GT2)];

covlevels = table(covlevels1,covlevels2, ...
    RowNames=["Sydney" "Melbourne" "Brisbane"], ...
    VariableNames=["Signal Strength at T1 (dBm)" "Signal Strength T2 (dBm)"]);

%%
function coveragedata = satcoverage(gridpts,sc,timeIn,inregion,beamWidth)

    % Get satellites and ground station receivers
    sats = sc.Satellites;
    rxs = [sc.GroundStations.Receivers];

    % Compute the latitude, longitude, and altitude of all satellites at the input time
    lla = states(sats,timeIn,"CoordinateFrame","geographic");

    % Initialize coverage data
    coveragedata = NaN(size(gridpts));

    for satind = 1:numel(sats)
        % Create a geopolyshape for the satellite field-of-view
        fov = fieldOfViewShape(lla(:,1,satind),beamWidth);

        % Find grid and rx locations which are within the field-of-view
        gridInFOV = isinterior(fov,gridpts);
        rxInFOV = gridInFOV(inregion);

        % Compute sigstrength at grid locations using temporary link objects
        gridsigstrength = NaN(size(gridpts));
        if any(rxInFOV)
            tx = sats(satind).Transmitters;
            lnks = link(tx,rxs(rxInFOV));
            rxsigstrength = sigstrength(lnks,timeIn)+30; % Convert from dBW to dBm
            gridsigstrength(inregion & gridInFOV) = rxsigstrength;
            delete(lnks)
        end

        % Update coverage data with maximum signal strength found
        coveragedata = max(gridsigstrength,coveragedata);
    end
end

function satFOV = fieldOfViewShape(lla,beamViewAngle)

    % Find the Earth central angle using the beam view angle
    if isreal(acosd(sind(beamViewAngle)*(lla(3)+earthRadius)/earthRadius))
        % Case when Earth FOV is bigger than antenna FOV 
        earthCentralAngle = 90-acosd(sind(beamViewAngle)*(lla(3)+earthRadius)/earthRadius)-beamViewAngle;
    else
        % Case when antenna FOV is bigger than Earth FOV 
        earthCentralAngle = 90-beamViewAngle;
    end

    % Create a buffer zone centered at the position of the satellite with a buffer of width equaling the Earth central angle
    [latBuff,lonBuff] = bufferm(lla(1),lla(2),earthCentralAngle,"outPlusInterior",100);

    % Handle the buffer zone crossing over -180/180 degrees
    if sum(abs(lonBuff) == 180) > 2
        crossVal = find(abs(lonBuff)==180) + 1;
        lonBuff(crossVal(2):end) = lonBuff(crossVal(2):end) - 360 *sign(lonBuff(crossVal(2)));
    elseif sum(abs(lonBuff) == 180) == 2
        lonBuff = [lonBuff; lonBuff(end); lonBuff(1); lonBuff(1)];
        if latBuff(1) > 0
            latBuff = [latBuff; 90; 90; latBuff(1)];
        else
            latBuff = [latBuff; -90; -90; latBuff(1)];
        end
    end

    % Create geopolyshape from the resulting latitude and longitude buffer zone values
    satFOV = geopolyshape(latBuff,lonBuff);
end

function GT = contourDataGrid(latd,lond,data,levels,proj)

    % Pad each side of the grid to ensure no contours extend past the grid bounds
    lond = [2*lond(1,:)-lond(2,:); lond; 2*lond(end,:)-lond(end-1,:)];
    lond = [2*lond(:,1)-lond(:,2) lond 2*lond(:,end)-lond(:,end-1)];
    latd = [2*latd(1,:)-latd(2,:); latd; 2*latd(end,:)-latd(end-1,:)];
    latd = [2*latd(:,1)-latd(:,2) latd 2*latd(:,end)-latd(:,end-1)];
    data = [ones(size(data,1)+2,1)*NaN [ones(1,size(data,2))*NaN; data; ones(1,size(data,2))*NaN] ones(size(data,1)+2,1)*NaN];

    % Replace NaN values in power grid with a large negative number
    data(isnan(data)) = min(levels) - 1000;
    
    % Project the coordinates using an equal-area projection
    [xd,yd] = projfwd(proj,latd,lond);
    
    % Contour the projected data
    fig = figure("Visible","off");
    obj = onCleanup(@()close(fig));
    c = contourf(xd,yd,data,LevelList=levels);
    
    % Initialize variables
    x = c(1,:);
    y = c(2,:);
    n = length(y);
    k = 1;
    index = 1;
    levels = zeros(n,1);
    latc = cell(n,1);
    lonc = cell(n,1);
    
    % Calculate the area within each contour line. Remove areas that
    % correspond to holes and ignore negative areas.
    while k < n
        level = x(k);
        numVertices = y(k);
        yk = y(k+1:k+numVertices);
        xk = x(k+1:k+numVertices);
        k = k + numVertices + 1;
    
        [first,last] = findNanDelimitedParts(xk);
        jindex = 0;
        jy = {};
        jx = {};
        sumpart = 0;
    
        for j = 1:numel(first)
            % Process the j-th part of the k-th level
            s = first(j);
            e = last(j);
            cx = xk(s:e);
            cy = yk(s:e);
            if cx(end) ~= cx(1) && cy(end) ~= cy(1)
                cy(end+1) = cy(1); 
                cx(end+1) = cx(1);
            end

            if j == 1 && ~ispolycw(cx,cy)
                % First region must always be clockwise
                [cx,cy] = poly2cw(cx,cy);
            end
    
            jindex = jindex + 1;
            jy{jindex,1} = cy(:)';
            jx{jindex,1} = cx(:)';
    
            a = polyarea(cx,cy);
            if ~ispolycw(cx,cy)
                % Remove areas that correspond to holes
                a = -a;
            end

            sumpart = sumpart + a;
        end
    
        % Add a part when its area is greater than 0. Unproject the
        % coordinates.
        [jx,jy] = polyjoin(jx,jy);
        if length(jy) > 2 && length(jx) > 2 && sumpart > 0
            [jlat,jlon] = projinv(proj,jx,jy);
            latc{index,1} = jlat(:)';
            lonc{index,1} = jlon(:)';
            levels(index,1) = level;
            index = index + 1;
        end
    end
    
    % Create contour shapes from the unprojected coordinates. Include the
    % contour shapes, the areas, and the power levels in a geospatial
    % table.
    latc = latc(1:index-1);
    lonc = lonc(1:index-1);
    Shape = geopolyshape(latc,lonc);
    Levels = levels(1:index-1);
    allPartsGT = table(Shape,Levels);  

    % Condense the geospatial table into a new geospatial table with one
    % row per power level.
    GT = table.empty;
    levels = unique(allPartsGT.Levels);
    for k = 1:length(levels)
        gtLevel = allPartsGT(allPartsGT.Levels == levels(k),:);
        tLevel = geotable2table(gtLevel,["Latitude","Longitude"]);
        [lon,lat] = polyjoin(tLevel.Longitude,tLevel.Latitude);
        Shape = geopolyshape(lat,lon);
        Levels = levels(k);
        GT(end+1,:) = table(Shape,Levels);
    end

    maxLevelDiff = max(abs(diff(GT.Levels)));
    LevelRange = [GT.Levels GT.Levels + maxLevelDiff];
    GT.LevelRange = LevelRange;

    % Clean up the geospatial table
    GT.Properties.VariableNames = ...
        ["Shape","Power (dBm)","Power Range (dBm)"]; 
end

function powerLevels = coveragelevel(lat,lon,GT)

    % Determine whether points are within coverage
    inContour = false(length(GT.Shape),1);
    for k = 1:length(GT.Shape)
        inContour(k) = isinterior(GT.Shape(k),geopointshape(lat,lon));
    end

    % Find the highest coverage level containing the point
    powerLevels = max(GT.("Power (dBm)")(inContour));

    % Return -inf if the point is not contained within any coverage level
    if isempty(powerLevels)
        powerLevels = -inf;
    end
end


function [first,last] = findNanDelimitedParts(x)

    % Find indices of the first and last elements of each part in x. 
    % x can contain runs of multiple NaNs, and x can start or end with 
    % one or more NaNs.

    n = isnan(x(:));
    
    firstOrPrecededByNaN = [true; n(1:end-1)];
    first = find(~n & firstOrPrecededByNaN);
    
    lastOrFollowedByNaN = [n(2:end); true];
    last = find(~n & lastOrFollowedByNaN);
end

