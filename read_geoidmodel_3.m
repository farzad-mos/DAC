function deltaz=read_geoidmodel_3(lat_point,lon_point,lat_geo,lon_geo,z_geo)
%find geoidal height of interested point and apply to measured dynamic
%topography
%example: loading loksa EST-GEOID2017 to find estimated TG data near interesting
%point

%finding nearest point to lat+-0.005 deg (1/2 of mesh) & lon_geo+-0.01 deg (1/2
%mesh)
% interpolate the geoid hight around 9 nearest scatter data points instead of choosing nearest one

%inster mesh size*0.5 (1/2 of mesh size)
lat_m=abs((lat_geo(1,1)-lat_geo(2,1))/2); 
lon_m=abs((lon_geo(1,1)-lon_geo(2,1))/2);

[x,~]=find(lat_geo<=lat_point+lat_m & lat_geo>lat_point-lat_m);
[y,~]=find(lon_geo<=lon_point+lon_m & lon_geo>lon_point-lon_m);

if isempty(x)
    [x,~]=find(lat_geo<=lat_point+lat_m+(lat_m/2) & lat_geo>lat_point-lat_m-(lat_m/2));
end

if isempty(y)
    [y,~]=find(lon_geo<=lon_point+lon_m+(lon_m/2) & lon_geo>lon_point-lon_m-(lon_m/2));
end

if length(x)>1
    x=x(1);
end

if length(y)>1
    y=y(1);
end

if isempty(x) || isempty(y)
    deltaz=NaN;
else

% interpolate the z_geo with a squre of 6 around points
% find a solution for boundy situation
if x<length(lat_geo) && y<length(lon_geo) && x~=1 && y~=1
    f=scatteredInterpolant(repmat(lat_geo(x-1:x+1),3,1),repelem(lon_geo(y-1:y+1),3),reshape(z_geo(y-1:y+1,x-1:x+1),[],1));
        deltaz=f(lat_point,lon_point);
elseif x==1 || y==1 && x~=length(lat_geo) && y~=length(lon_geo)
    f=scatteredInterpolant(repmat(lat_geo(x:x+1),2,1),repelem(lon_geo(y:y+1),2),reshape(z_geo(y:y+1,x:x+1),[],1));
    f.Method = 'natural';
    deltaz=f(lat_point,lon_point);
elseif x==length(lat_geo) || y==length(lon_geo) && x~=1 && y~=1
    f=scatteredInterpolant(repmat(lat_geo(x-1:x),2,1),repelem(lon_geo(y-1:y),2),reshape(z_geo(y-1:y,x-1:x),[],1));
    deltaz=f(lat_point,lon_point);
else
    deltaz=z_geo(y,x);
end
end

