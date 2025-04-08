% This script read preprocessed SA data and compare it with Nemo-Nordic Hydrodynamic Model


% Inputs: You need to load below files:

% functions
% ConvertSerialYearToDate.m function to convert decilam year to date format
% ellipcor.m function for ellipsoidal vorrection (newer version is needed)
% grdread2.m function to read .grd file (NKG2015 format)
% read_geoidmodel_3.m to interpolate geoidal height at specific location

% data:
% pre processed SA data (.mat files) reading .nc file using related function for each data set e.g. balticseal_v02.m function to read Baltic+SEAL data
% EST_GEOID2017.nc geoid file to measure geoidal height
% NKG2015.grd geoid file to measure geoidal height
% HBM-EST hydrodynamic files folder location
% Nemo-Nordic files folder location
% BS.shp baltic sea coast shape file 
% TGdata.mat tide gaguges records
% tg_s3a.mat, tg_s3a.mat, tg_s3b.mat tidegauge number per each pass


% version 2:
% - interpolate grid data of Nemo-Nordic and HBM-EST HDM 
% - interpolate NKG2015 and EST2017 geoidal height at specific point
% - Topex/Posidon to WGS84 height difference
% - compare both HMD moldels by both geoid models
% - scatter plot of each SA at each year/overpass
% - correct NEMO-NORDIC model by TG (average of 6hrs) using weight (based
% on distance)
% - correct SA based on corrected HDM model
% - save figures of each pass/cycle

% Outputs
% 1-updated SA .mat file with added columns
% 2-res.mat for statisctical results
% 3-saved figures
 

 

%%%%%%%%%---------===========================================---------%%%%%%%%%
% Matlab R2019b used to write this code
% Majid Mostafavi, farzad.mostafavi@live.com, TalTech, Estonia 2021.March



clear
close all
clc


% Loading shoreline model

% change later

%% data processing

sa=load(uigetfile('*.mat',...
              'SA data'));         
storedvars = fieldnames(sa) ;
FirstVarName = storedvars{1};
FirstVarContent = sa.(FirstVarName);

          
%sort based on time
alldata=cell2mat(FirstVarContent);
alldata=sortrows(alldata,6); %sort data based on pass no.
alldata=sortrows(alldata,3); %sort data based on time

[m,~]=find(isnan(alldata(:,1)));
alldata(m,:)=[];
clear m

%remove NaN data
[m,~]=find(isnan(alldata(:,3)));
alldata(m,:)=[];
clear m

%remove negative data
[m,~]=find(alldata(:,4)<=0);
alldata(m,:)=[];
clear m



% convert to wgs84 ellipsoid
for i=1:length(alldata)
    alldata(i,4)=alldata(i,4)+ellipcor(alldata(i,1));
end

 
 clear storedvars
 clear FirstVarContent
 clear sa
 clear FirstVarName
     
%% Convert 2 Cell
% change track_all_cor to cell format (track_all_cor_cell) at each
% year/pass/cycle

    for i=min(alldata(:,6)):max(alldata(:,6)) %finding each pass no
        [m,~]=find(i==alldata(:,6));
            all_cor{i,:}=alldata(m,:);      % creat a cell of each pas per each pass
    end
    all_cor=all_cor(~cellfun('isempty',all_cor)); %remove empty cells (in case i did not start from 1)  
 
clear i
clear m


for j=1:length(all_cor)
    for i=min(floor(all_cor{j}(:,3))):max(floor(all_cor{j}(:,3))) %finding each year
        [m,~]=find(i==floor(all_cor{j}(:,3)));
        all_cor_p{i,:}=all_cor{j}(m,:);
    end
    all_cor_p=all_cor_p(~cellfun('isempty',all_cor_p)); %remove empty cells (in case i did not start from 1)
    all_cor_p_y{j,1}=all_cor_p; % creat a cell of each pas per each pass/year
    clear all_cor_p
end


clear i
clear j
clear m


for j=1:length(all_cor)
    for k=1:length(all_cor_p_y{j})
        for i=min(all_cor_p_y{j}{k}(:,6)):max(all_cor_p_y{j}{k}(:,6)) %finding each cycle
            [m,~]=find(i==all_cor_p_y{j}{k}(:,6));
            all_cor_p{i,:}=all_cor_p_y{j}{k}(m,:);
        end
        all_cor_p=all_cor_p(~cellfun('isempty',all_cor_p)); %remove empty cells (in case i did not start from 1)
        track_all_cor_cell{j,1}{k,1}=all_cor_p; % creat a cell of each pas per each pass/year/cycle
        clear all_cor_p
    end
end


clear i
clear j
clear m
clear all_cor
clear all_cor_p_y

%% Reading HDM data (Nemo-Nordic)
% this section will read Nemo-Nordic HDM at each SA point location and add
% to the 7th column of pre-readed SA matrix

%
% find the coast location form shape file
%     lat_coast=str2num(answer{1});
%     lon_coast=str2num(answer{2});
%

%save current folder directory for future use
oldFolder=pwd;

%Select the Model Data folder
name = uigetdir('./subfolder1/');
%Change the file directory to Model Data to read them
cd(name);

% Read all Model Data
nc=ls('*.nc');
% 	info=ncinfo(nc(1,:));% to find the content of nc files

base=datenum(1950,1,1);
bar=waitbar(0,'Extracting Nemo-Nordic HDM data...');

for k=1:length(track_all_cor_cell)
    
    Fraction= floor((k/length(track_all_cor_cell)*100));
    formatSpec = 'Reading model data...%d%% complete';
    str2=sprintf(formatSpec,Fraction);
    waitbar(Fraction/100,bar,str2)
    
    for i=1:length(track_all_cor_cell{k})
        
        for j=1:length(track_all_cor_cell{k}{i})
            t_sa=datevec(ConvertSerialYearToDate(track_all_cor_cell{k}{i}{j}(1,3)));
            for w=1:length(nc(:,1))  
                if str2double(nc(w,11:14))==t_sa(1) && str2double(nc(w,15:16))==t_sa(2) && str2double(nc(w,17:18))==t_sa(3)
                    
                    
                    lat=(ncread(nc(w,:),'latitude'));
                    lon=(ncread(nc(w,:),'longitude'));
                    elev=(ncread(nc(w,:),'SSH_inst'));
                    time=(ncread(nc(w,:),'time_instant'));
                    t=datevec((time/86400)+base);
                    
                    % clear HBM bias
                    [n,~]=find(t(:,1)==0);
                    time(n,:)=[];
                    t(n,:)=[];
                    elev(:,:,n)=[];
                    clear n
                    
                    [m,~]=find(t(:,4)==t_sa(4));
                    
                    
                    if isempty(m)
                        lat=(ncread(nc(w+1,:),'latitude'));
                        lon=(ncread(nc(w+1,:),'longitude'));
                        elev=(ncread(nc(w+1,:),'SSH_inst'));
                        time=(ncread(nc(w+1,:),'time_instant'));
                        t=datevec((time/86400)+base);
                        [m,~]=find(t(:,4)==t_sa(4));
                    end
                    
                    
                    if isempty(m)
                        lat=(ncread(nc(w+2,:),'latitude'));
                        lon=(ncread(nc(w+2,:),'longitude'));
                        elev=(ncread(nc(w+2,:),'SSH_inst'));
                        time=(ncread(nc(w+2,:),'time_instant'));
                        t=datevec((time/86400)+base);
                        [m,~]=find(t(:,4)==t_sa(4));
                    end
                    
                    
                    if length(m)>1
                        m=m(1);
                    end
                    
                    
                    lon2=repmat(lon,1,length(lat));
                    lat2=repmat(lat',length(lon),1);
                    
                    hdm_dt_sa=griddata(lon2,lat2,elev(:,:,m),track_all_cor_cell{k}{i}{j}(:,2),track_all_cor_cell{k}{i}{j}(:,1),'cubic');
                    track_all_cor_cell{k}{i}{j}(:,8)=hdm_dt_sa;
                    
                    clear lat
                    clear lat2
                    clear lon
                    clear lon2
                    clear time
                    clear elev
                    clear t
                    clear hdm_dt_sa
                    clear m
                    
%                 else
%                     track_all_cor_cell{k}{i}{j}(:,8)=NaN;
                end
                
            end
            clear w
            clear t_sa
        end
        clear j
    end
    clear i 
end
clear k
close(bar);
msgbox('Reading Completed!','Status');
clear bar
% change folder back to previous folder
cd(oldFolder)

clear base
clear oldFolder
clear name
clear nc
clear Fraction
clear formatSpec
clear str2

%% add geoid height (NKG2015)
% this section will add the geoid height to sea level of HDM to obtain HDM
% drived SSH

%  NKG 2015 geoid using grid read function
    [lon_geo, lat_geo, z_geo]=grdread2('NKG2015.grd');
lon_geo=lon_geo';
lat_geo=lat_geo';
z_geo=double(z_geo)';


for k=1:length(track_all_cor_cell)
    for j=1:length(track_all_cor_cell{k})
        for i=1:length(track_all_cor_cell{k}{j})
            for l=1:size(track_all_cor_cell{k}{j}{i},1)
                
                
                shift=read_geoidmodel_3(track_all_cor_cell{k}{j}{i}(l,1),track_all_cor_cell{k}{j}{i}(l,2),lat_geo,lon_geo,z_geo);
                track_all_cor_cell{k}{j}{i}(l,9)=track_all_cor_cell{k}{j}{i}(l,8)+shift;

                clear shift
            end
        end
    end
end

clear k
clear l
clear i
clear j
clear lat_geo
clear lon_geo
clear z_geo

    
    
    %% Model correction & SA bias removal  
load('tg_s3a.mat') %load TG stations number per pass 
load('TGdata.mat') %load TG records   
%     load('S3A.mat')
    
    [lon_geo, lat_geo, z_geo]=grdread2('NKG2015.grd');
    lon_geo=lon_geo';
    lat_geo=lat_geo';
    z_geo=double(z_geo)';
    
    
    f=1;
    %

    for k=1:length(track_all_cor_cell)
         for i=1:length(track_all_cor_cell{k})
             for j=1:length(track_all_cor_cell{k}{i})
                
                % Extract geoidal height to measure DT
                for q=1:size(track_all_cor_cell{k}{i}{j},1)
                    N=read_geoidmodel_3(track_all_cor_cell{k}{i}{j}(q,1),track_all_cor_cell{k}{i}{j}(q,2),lat_geo,lon_geo,z_geo);
                    
                    if ~isnan(N) && ~isempty(N)
                        track_all_cor_cell{k}{i}{j}(:,11)=N;
                    else
                        track_all_cor_cell{k}{i}{j}(:,11)=NaN;
                    end
                    clear N
                end
                
                
                ff(f)=figure;
                
                subplot(6,2,9:12)
                track_selected=track_all_cor_cell{k}{i}{j};
                
                
                % remove NaN value
                [n1,~]=find(isnan(track_selected(:,9)));
                loc=track_selected(:,1:2);
                loc(n1,:)=NaN;
                
                
                
                % find nearest Model points in track
                a=isnan(loc(:,1));
                J=1;
                for I=1:length(a)
                    if a(I)==0 && I~=1 && a(I-1)==1
                        index(J,1)=I;
                        J=J+1;
                    elseif a(I)==0 && I~=length(a) && a(I+1)==1
                        index(J,1)=I;
                        J=J+1;
                    end
                end
                
                clear J
                clear I
                clear a
                
                if exist('index')
                    [m,~]=find(datenum(TGdate(:,1))<(ConvertSerialYearToDate(track_all_cor_cell{k}{i}{j}(1,3))+0.125) & datenum(TGdate(:,1))>=(ConvertSerialYearToDate(track_all_cor_cell{k}{i}{j}(1,3))-0.125));
                    
                    x=1;
                    for y=2:17
                        if tg_s3a(k,y)~=0 % change tg_s3b
                            tg_selected(x,:)=[TGinfo.Lat(tg_s3a(k,y)),TGinfo.Lon(tg_s3a(k,y)),(mean(TGrecords(m,tg_s3a(k,y)))*0.01),tg_s3a(k,y)];
                            x=x+1;
                        end
                    end
                    
                    clear x
                    clear y
                    
                    for y=1:length(index)
                        
                        dis1=distance([loc(index(y),1),loc(index(y),2)],[tg_selected(:,1),tg_selected(:,2)],referenceEllipsoid('WGS84'))/1000; %measure TG to SA points disance
                        [m3,~]=find(dis1>130); %exclude tidegaguges more than 50km distance
                        dis1(m3)=nan;
                        p=min(dis1)./dis1; %measure weight based on distance
                        shift=read_geoidmodel_3(loc(index(y),1),loc(index(y),2),lat_geo,lon_geo,z_geo); % measure SA point geoid height
                        tg_sa(y,1)=(sum((p.*tg_selected(:,3)),'omitnan')/sum(p,'omitnan'))+shift; %mean of TGs
                        bias_m(y,1)=track_selected(index(y),9)-tg_sa(y,1); %TG-HDM bias
                        
                        % plot TG amount
                        for y2=1:length(dis1)
                            if ~isnan(dis1(y2))
                                h(1)=plot(loc(index(y),1),tg_selected(y2,3)+shift, '^k','MarkerFaceColor','k','MarkerSize',12,'DisplayName','TG');
                                hold on
                                %                 text(loc(index(y),1),(tg_selected(y2,3)+shift)-0.02, num2str(tg_selected(y2,4)),'FontSize',13,'Color',[0.7 0.7 0.7],'FontWeight','Bold');
                            end
                        end
                        
                        %         plot(loc(index(y),1),tg_sa(y,1),'or') %remove later
                        
                        clear dis1
                        clear m3
                        clear p
                        clear shift
                        clear y2
                    end
                    
                    [n2,~]=find(isnan(bias_m));
                    bias_m(n2)=[];
                    index(n2)=[];
                    
                    
                    clear y
                    
                    
                    h(2)=plot(loc(:,1),smooth(track_selected(:,9),'rlowess'),'-r','DisplayName','NEMO_N_K_G', 'LineWidth',2);
                    
                    
                    
                    if ~isempty(bias_m)
                        if length(bias_m)==1
                            bias_hdm=bias_m;
                        else
                            bias_hdm=interp1(track_selected(index,1),bias_m,track_selected(:,1));
                        end
                        
                        track_all_cor_cell{k}{i}{j}(:,10)=track_all_cor_cell{k}{i}{j}(:,9)-bias_hdm;
                        track_selected(:,10)=track_selected(:,9)-bias_hdm;
                        
                                              
                        
                        if ~isempty(~isnan(track_selected(:,10)))
                            
                            h(3)=plot(loc(:,1),smooth(track_selected(:,10),'rlowess'),'-','color','#77AC30','DisplayName','NEMO-corr_N_K_G','LineWidth',2);
                            
                            %measure SA bias
                            track_limited=track_selected;
                            [n2,~]=find(isnan(track_limited(:,10)));
                            track_limited(n2,:)=[]; %remove NaN values
                            
                            [m2,~]=find(abs(track_limited(:,4)-track_limited(:,10))<3.0); %consider points lower than 3m diff only
                            bias_sa=mean((track_limited(m2,4)-track_limited(m2,11))-(track_limited(m2,10)-track_limited(m2,11)));
                            
                            track_all_cor_cell{k}{i}{j}(:,12)= track_all_cor_cell{k}{i}{j}(:,4)-bias_sa;
                            track_selected(:,12)=track_selected(:,4)-bias_sa;
                            
                            
                            [m3,~]=find(abs(track_limited(:,5)-track_limited(:,10))<3.0); %consider points lower than 3m diff only
                            bias_sa2=mean((track_limited(m3,5)-track_limited(m3,11))-(track_limited(m3,10)-track_limited(m3,11)));
                            
                            track_all_cor_cell{k}{i}{j}(:,13)= track_all_cor_cell{k}{i}{j}(:,5)-bias_sa2;
                            track_selected(:,13)=track_selected(:,5)-bias_sa2;                           
                            
                            
                        else
                            
                            track_all_cor_cell{k}{i}{j}(:,12)=NaN;
                            track_all_cor_cell{k}{i}{j}(:,13)= NaN;
                            bias_sa=NaN;
                            bias_hdm=NaN;
                            
                            track_selected(:,12)=NaN;
                            track_selected(:,13)=NaN;
                        end
                    else
                        track_all_cor_cell{k}{i}{j}(:,12)=NaN;
                        track_all_cor_cell{k}{i}{j}(:,13)= NaN;
                        bias_sa=NaN;
                        bias_hdm=NaN;
                        
                        track_selected(:,13)=NaN;
                        track_selected(:,12)=NaN;
                    end
                    
                    
                    
                    h(4)=plot(track_selected(:,1),track_selected(:,4),'.','color',	'#D95319','DisplayName','ALES+', 'LineWidth',2);
                    
                    sd_sa_hdm=std(abs(track_selected(:,12)-track_selected(:,10)),'omitnan');
                    %
                    [m1,~]=find(abs(track_selected(:,12)-track_selected(:,10))<=sd_sa_hdm);
                    if ~isempty(m1)
                        h(5)=plot(track_selected(m1,1),track_selected(m1,12),'.','color','#0072BD','DisplayName','ALES+_c_o_r_r','LineWidth',2);
                    end
                    
                    if ~isempty(bias_m)
                        legend(h);
                    end
                    
                    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                    %     title(strcat('Cycle:' , num2str(track_all_cor_cell{k}{i}{j}(1,5)),' - ', datestr(ConvertSerialYearToDate(track_all_cor_cell{k}{j}{i}(1,3)),'dd mmm yyyy')),'FontSize',14,'FontWeight','bold');
                    
                    title(strcat('Cycle:' , num2str(track_selected(1,5)),' - ', datestr(ConvertSerialYearToDate(track_selected(1,3)),'dd mmm yyyy')),'FontSize',14,'FontWeight','bold');
                    
                    ylabel('SSH [m]','FontSize',14,'FontWeight','bold');
                    xlabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                    
                    
                    hold off
                    
                else
                    track_all_cor_cell{k}{i}{j}(:,13)=NaN;
                    track_all_cor_cell{k}{i}{j}(:,12)= NaN;
                    bias_sa=NaN;
                    bias_hdm=NaN;
                end
                
                clear n1
                clear ax
                clear track_selected
                clear h
                clear m
                clear index
                clear m1
                clear bias_m
                clear bias_hdm
                clear tg_sa
                clear hdm_sa
                clear loc
                clear Index
                clear n2
                clear track_limited
                clear m2
                clear m1
                clear tg_sa
                clear tg_selected
                clear index
                clear sd_sa_hdm
                clear q
                clear N
                
                subplot(6,2,[1,3,5,7])
                for y=2:17
                    if tg_s3a(k,y)~=0 % change tg_s3b
                        h(1)=plot(TGinfo.Lat(tg_s3a(k,y)),TGinfo.Lon(tg_s3a(k,y)), '^k','MarkerFaceColor','k','MarkerSize',9,'DisplayName','TG');
                        hold on
                        text(TGinfo.Lat(tg_s3a(k,y))-0.07,TGinfo.Lon(tg_s3a(k,y)), num2str(TGinfo.TGId(tg_s3a(k,y))),'FontSize',13,'Color','r','FontWeight','Bold');
                    end
                end
                
                
                track_selected=track_all_cor_cell{k}{i}{j};
                [m,~]=find(isnan(track_selected(:,13)));
                track_selected(m,1)=nan;
                
                
                sd_sa_hdm=std(abs(track_selected(:,12)-track_selected(:,10)),'omitnan');
                
                [m1,~]=find(abs(track_selected(:,12)-track_selected(:,10))>sd_sa_hdm);
                [m2,~]=find(abs(track_selected(:,12)-track_selected(:,10))<=sd_sa_hdm);
                
                
                scatter(track_selected(m2,1),track_selected(m2,2),25,(track_selected(m2,12)-track_selected(m2,10))*100,'filled');
                c=colorbar;
                caxis([-70 70]); %fix the bar
                c.Label.String = '[cm]';
                
                
                
                
                %
                %     if ~isempty(m1)
                %     h(2)=plot(track_selected(m1,1),track_selected(m1,2),'xr','MarkerSize',9,'DisplayName', strcat('No. outliers:',num2str(length(m1))), 'LineWidth',2);
                %     end
                
                rmse_sa_hdm=rms(abs(track_selected(m2,12)-track_selected(m2,10)),'omitnan');
                x=track_selected(:,10);
                
                
                leg=legend(h,'Box','off','Location','southwest');
                title(leg,strcat('RMSE=',num2str((rmse_sa_hdm*100),2),' cm','\newlineNo. outliers:',num2str(length(m1)),...
                    ' (',num2str((length(m1)/length(x(~isnan(x))))*100,2),'% )'),'FontSize',12,'FontWeight','bold');
                
                
                xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
                ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                
                
                if ~isnan(sd_sa_hdm)
                    ax2 = axes('Position',[0.3234375,0.750521711657832,0.094755595025878,0.155759444872576],'OuterPosition',...
                        [0.3075,0.7202,0.1223,0.2004],'InnerPosition',[0.3234375,0.750521711657832,0.094755595025878,0.155759444872576]);
                    box on;
                    hold(ax2, 'on');
                    ax3=gca; ax3.GridAlpha = 0.3; ax3.FontSize=14; ax3.FontWeight='Bold'; grid on
                    
                    histogram(track_selected(:,12)-track_selected(:,10),'FaceColor' ,[0 0.5 0.5]);
                    xline(sd_sa_hdm,'--','STD','color','r');
                    xline(-sd_sa_hdm,'--','STD','color','r');
                end
                
                
                %         plot(track_selected(m2,1),track_selected(m2,14)-(track_selected(m2,13)),'-','color','#0072BD','LineWidth',1.5)
                %         hold on
                %
                %         plot(track_selected(m1,1),track_selected(m1,14)-(track_selected(m1,13)),'x','color','#A2142F')
                %         yline(0,'--','HMD','color',	'#77AC30','LineWidth',1);
                
                
                % save result in res including: date, cycle, pass#, no out., no data, bias_hdm, bias_sa, rmse, std
                res{k,1}{i,1}(j,:)=[track_selected(1,3),track_selected(1,6),track_selected(1,7),length(m1),length(x(~isnan(x))),bias_sa,rmse_sa_hdm,sd_sa_hdm];
                
                %figure name
                str2=strcat('ja3','pass',num2str(track_selected(1,6)),'c',num2str(track_selected(1,5)));
                
                hold off
                clear m1
                clear m2
                clear m
                clear ax
                clear track_selected
                clear h
                clear ax2
                clear ax3
                clear sd_sa_hdm
                clear leg
                clear x
                clear bias_hdm
                clear bias_sa
                clear bias_sa2

                
                subplot(6,2,[2,4,6,8])
                geoplot(TGinfo.Lat,TGinfo.Lon, '^k','MarkerSize',8)
                hold on
                %     text(TGinfo.Lat,TGinfo.Lon+0.3, num2str(TGinfo.TGId),'FontSize',11,'Color',[0.7 0.7 0.7]);
                h(1)=geoplot(track_all_cor_cell{k}{i}{j}(:,1),track_all_cor_cell{k}{i}{j}(:,2),'o','color','#7E2F8E','MarkerFaceColor','#7E2F8E','DisplayName',strcat('Pass#',num2str(track_all_cor_cell{k}{i}{j}(1,7))));
                
                for y=2:17
                    if tg_s3a(k,y)~=0 % change tg_s3b
                        h(2)=geoplot(TGinfo.Lat(tg_s3a(k,y)),TGinfo.Lon(tg_s3a(k,y)), '^k','MarkerFaceColor','k','MarkerSize',8,'DisplayName','TG');
                        text(TGinfo.Lat(tg_s3a(k,y)),TGinfo.Lon(tg_s3a(k,y))+0.3, num2str(TGinfo.TGId(tg_s3a(k,y))),'FontSize',13,'Color','red','FontWeight','Bold');
                    end
                end
                
                %     plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast
                
                
                geobasemap landcover
                
                %
                %     xlabel('Lon [\circ]','FontSize',14);
                %     ylabel('Lat [\circ]','FontSize',14);
                %
                %     [leg,icons]=legend(h);
                %     icons(4).MarkerSize = 25;
                %     icons(6).MarkerSize = 11;
                %     icons(1).FontSize=14;
                %     icons(2).FontSize=14;
                leg=legend(h);
                set(leg,'Box','off','Location','northwest')
                %     legend boxoff
                %     leg.Location='northwest';
                ax=gca;  ax.FontSize=14; ax.FontWeight='Bold';
                
                %     pbaspect([0.9 1 1])
                
                %             clear icons
                clear h
                clear leg
                
                hold off
                %
                %     ax2 = axes('Position',[.7 .7 .2 .2]);
                %     box on;
                %     hold(ax2, 'on');
                % %     geolimits([59 60.8],[22.5 30.5]);
                %      scatter(track_all_cor_cell{k}{i}{j}(:,1),track_all_cor_cell{k}{i}{j}(:,2),25,(track_all_cor_cell{k}{i}{j}(:,4)-track_all_cor_cell{k}{i}{j}(:,9)),'filled')
                %      c=colorbar;
                
                
                
                %  pass(k,1)=track_all_cor_cell{k}{i}{j}(1,6);
                
                
                f=f+1
                %             if ~isnan(rmse_sa_hdm)
                %              savefig(ff,str2);
                %               clear ff
                % %             saveas(ff,str2,'fig')
                %             end
                %                           clear ff
                
                clear rmse_sa_hdm
                %             clear str2           
                 close all
             end
         end
    end
%


%% result comparision
    
% load res files and convert to table
s3a=load('res_s3a.mat');
s3a=s3a.res;


for k=1:length(s3a)
if k==1
    sa_all=cell2mat(s3a{k});
else
    sa_all=[sa_all;cell2mat(s3a{k})];
end
end
[m,~]=find(isnan(sa_all(:,6)));
sa_all(m,:)=[];


date_dec=sa_all(:,1);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,1))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,1)));
year=t_sa(:,1);
month=t_sa(:,2);


for i=1:length(t_sa)
    if t_sa(i,2)==3 || t_sa(i,2)==4 || t_sa(i,2)==5
season(i,1)=1; %spring
    elseif  t_sa(i,2)==6 || t_sa(i,2)==7 || t_sa(i,2)==8
season(i,1)=2; %summer
    elseif  t_sa(i,2)==9 || t_sa(i,2)==10 || t_sa(i,2)==11
season(i,1)=3; %autum
    elseif t_sa(i,2)==12 || t_sa(i,2)==1 || t_sa(i,2)==2
        season(i,1)=4; %spring
    end
end


cycle=sa_all(:,2);
pass=sa_all(:,3);
out=sa_all(:,4);
obs=sa_all(:,5);
bias=sa_all(:,6);
rmse=sa_all(:,7);
sd=sa_all(:,8);


S3A = table(date,date_dec,year,month,season,cycle,pass,out,obs,bias,rmse,sd);

clear k
clear m
clear date
clear date_dec
clear year
clear month
clear season
clear cycle
clear pass
clear out
clear obs
clear bias
clear rmse
clear sd
clear i
clear t_sa
clear sa_all


s3b=load('res_s3b.mat');
s3b=s3b.res;


for k=1:length(s3b)
if k==1
    sa_all=cell2mat(s3b{k});
else
    sa_all=[sa_all;cell2mat(s3b{k})];
end
end
[m,~]=find(isnan(sa_all(:,6)));
sa_all(m,:)=[];


date_dec=sa_all(:,1);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,1))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,1)));
year=t_sa(:,1);
month=t_sa(:,2);

for i=1:length(t_sa)
    if t_sa(i,2)==1 || t_sa(i,2)==2 || t_sa(i,2)==3
        season(i,1)=1;
    elseif  t_sa(i,2)==4 || t_sa(i,2)==5 || t_sa(i,2)==6
        season(i,1)=2;
    elseif  t_sa(i,2)==7 || t_sa(i,2)==8 || t_sa(i,2)==9
        season(i,1)=3;
    else
        season(i,1)=4;
    end
end
cycle=sa_all(:,2);
pass=sa_all(:,3);
out=sa_all(:,4);
obs=sa_all(:,5);
bias=sa_all(:,6);
rmse=sa_all(:,7);
sd=sa_all(:,8);


S3B = table(date,date_dec,year,month,season,cycle,pass,out,obs,bias,rmse,sd);

clear k
clear m
clear date
clear date_dec
clear year
clear month
clear season
clear cycle
clear pass
clear out
clear obs
clear bias
clear rmse
clear sd
clear i
clear t_sa
clear sa_all


ja3=load('res_ja3.mat');
ja3=ja3.res;


for k=1:length(ja3)
if k==1
    sa_all=cell2mat(ja3{k});
else
    sa_all=[sa_all;cell2mat(ja3{k})];
end
end
[m,~]=find(isnan(sa_all(:,6)));
sa_all(m,:)=[];




date_dec=sa_all(:,1);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,1))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,1)));
year=t_sa(:,1);
month=t_sa(:,2);

for i=1:length(t_sa)
    if t_sa(i,2)==1 || t_sa(i,2)==2 || t_sa(i,2)==3
season(i,1)=1;
    elseif  t_sa(i,2)==4 || t_sa(i,2)==5 || t_sa(i,2)==6
season(i,1)=2;
    elseif  t_sa(i,2)==7 || t_sa(i,2)==8 || t_sa(i,2)==9
season(i,1)=3;
    else
        season(i,1)=4;
    end
end
cycle=sa_all(:,2);
pass=sa_all(:,3);
out=sa_all(:,4);
obs=sa_all(:,5);
bias=sa_all(:,6);
rmse=sa_all(:,7);
sd=sa_all(:,8);


JA3 = table(date,date_dec,year,month,season,cycle,pass,out,obs,bias,rmse,sd);

clear k
clear m
clear date
clear date_dec
clear year
clear month
clear season
clear cycle
clear pass
clear out
clear obs
clear bias
clear rmse
clear sd
clear i
clear t_sa
clear sa_all

clear s3a
clear s3b
clear ja3




%% compare results
% S3B.year(536) = 2018;
% S3B.month(536) = 0;
% S3B.season(536) = 0;
% S3B.cycle(536) = 0;
% S3B.month(536) = 1;
% S3B.pass(536) = 0;
% S3B.out(536) = 0;
% S3B.obs(536) = 0;
% S3B.bias(536) = 0;
% S3B.rmse(536) = 0;
% S3B.sd(536) = 0;
% S3B.Var13(536,1) = 0;
% S3B.month(537) = 2;
% S3B.month(538) = 3;
% S3B.month(539) = 4;
% S3B.month(540) = 5;
% S3B.year(537) = 2018;
% S3B.year(538) = 2018;
% S3B.year(539) = 2018;
% S3B.year(540) = 2018;

% RMSE

y1=splitapply(@mean, S3A.rmse(S3A.year==2017),S3A.month(S3A.year==2017))*100;
y2=splitapply(@mean, S3A.rmse(S3A.year==2018),S3A.month(S3A.year==2018))*100;
y3=splitapply(@mean, S3A.rmse(S3A.year==2019),S3A.month(S3A.year==2019))*100;
[~,ge] = findgroups(S3A.month(S3A.year==2019));



y1=splitapply(@mean, JA3.rmse(JA3.year==2017),JA3.month(JA3.year==2017))*100;
y2=splitapply(@mean, JA3.rmse(JA3.year==2018),JA3.month(JA3.year==2018))*100;
y3=splitapply(@mean, JA3.rmse(JA3.year==2019),JA3.month(JA3.year==2019))*100;

y2=splitapply(@mean, S3B.rmse(S3B.year==2018),S3B.month(S3B.year==2018))*100;
y3=splitapply(@mean, S3B.rmse(S3B.year==2019),S3B.month(S3B.year==2019))*100;


% BIAS

y1=splitapply(@mean, S3A.bias(S3A.year==2017),S3A.month(S3A.year==2017))*100;
y2=splitapply(@mean, S3A.bias(S3A.year==2018),S3A.month(S3A.year==2018))*100;
y3=splitapply(@mean, S3A.bias(S3A.year==2019),S3A.month(S3A.year==2019))*100;


y1=splitapply(@mean, JA3.bias(JA3.year==2017),JA3.month(JA3.year==2017))*100;
y2=splitapply(@mean, JA3.bias(JA3.year==2018),JA3.month(JA3.year==2018))*100;
y3=splitapply(@mean, JA3.bias(JA3.year==2019),JA3.month(JA3.year==2019))*100;

y2=splitapply(@mean, S3B.bias(S3B.year==2018),S3B.month(S3B.year==2018))*100;
y3=splitapply(@mean, S3B.bias(S3B.year==2019),S3B.month(S3B.year==2019))*100;




 % No. out
y1=splitapply(@sum, S3A.out(S3A.year==2017),S3A.month(S3A.year==2017));
y2=splitapply(@sum, S3A.out(S3A.year==2018),S3A.month(S3A.year==2018));
y3=splitapply(@sum, S3A.out(S3A.year==2019),S3A.month(S3A.year==2019));


y4=splitapply(@sum, JA3.out(JA3.year==2017),JA3.month(JA3.year==2017));
y5=splitapply(@sum, JA3.out(JA3.year==2018),JA3.month(JA3.year==2018));
y6=splitapply(@sum, JA3.out(JA3.year==2019),JA3.month(JA3.year==2019));


y7=splitapply(@sum, S3B.out(S3B.year==2018),S3B.month(S3B.year==2018));
y8=splitapply(@sum, S3B.out(S3B.year==2019),S3B.month(S3B.year==2019));


y=[y1,y2];
b = bar(y,1.5);
b(3).FaceColor = [0.4660 0.6740 0.1880];
b(2).FaceColor = [0.6350 0.0780 0.1840];
b(1).FaceColor = [0 0.4470 0.7410]	;



gscatter(S3A.month(S3A.year==2017),S3A.bias(S3A.year==2017),S3A.pass(S3A.year==2017))
hold on
gscatter(S3A.month(S3A.year==2018),S3A.bias(S3A.year==2018),S3A.pass(S3A.year==2018))
gscatter(S3A.month(S3A.year==2019),S3A.bias(S3A.year==2019),S3A.pass(S3A.year==2019))



gscatter(S3A.season,S3A.bias,S3A.pass)
%% plot gross errors


f=1;
figure(1)
            latlim = [53 67 ];
            lonlim = [10 31];

            ax = usamap(latlim, lonlim);
            geoshow('landareas.shp','FaceColor',[0.5 0.7 0.5])
            
            for k=1:length(track_all_cor_cell)
                for i=1:length(track_all_cor_cell{k})
                    for j=1:length(track_all_cor_cell{k}{i})
                        sd_sa_hdm=JA3.sd(JA3.pass==track_all_cor_cell{k}{i}{j}(1,6) & JA3.cycle==track_all_cor_cell{k}{i}{j}(1,5));
                        if ~isempty(sd_sa_hdm)
                            [m1,~]=find(abs(track_all_cor_cell{k}{i}{j}(:,4)-track_all_cor_cell{k}{i}{j}(:,13))>3);
                            for l=1:length(m1)
                                geoshow(track_all_cor_cell{k}{i}{j}(m1(l),1),track_all_cor_cell{k}{i}{j}(m1(l),2),'Marker', 'o','Color','r','MarkerFaceColor','r','MarkerSize',2)
                            end
                            
                            hold on
                            clear m1
                            
                        end
                        clear sd_sa_hdm
                        f=f+1
                    end
                end
            end

setm(gca,'FLineWidth',5,'Grid','on','FontSize',14,'fontweight','bold')


xLoc =1.3001e+05;
yLoc =7.0876e+06;

scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc);




%% scatter plot (before and after correction)

% S = shaperead('BS.shp');
[lat2,lon2,z] = read_kml('coast.kml');


% dt=0.074;
dt=0.0272;

% for JA3: 0.0272
% for S3: 0.074

str='JA3';
% str='S3A';
% str='S3B';


f1=1;
f2=1;


for k=2:length(track_all_cor_cell{2})
    for i=1:length(track_all_cor_cell)
        if i==1
            track_all=cell2mat(track_all_cor_cell{i}{k});
        end
        track_all=[track_all;cell2mat(track_all_cor_cell{i}{k})];
    end
    
    
    mi=min(track_all(:,3));
    [le,~]=find(~isnan(track_all(:,9)));
    l=1;
    
    
    while mi<=max(track_all(:,3))
        [m,~]=find(track_all(:,3)>=mi & track_all(:,3)<mi+dt);
        if ~isempty(m)
            track_selected=track_all(m,:);
            
            [m,~]=find(isnan(track_selected(:,9)));
            track_selected(m,:)=[];
            
            avg=mean(abs(track_selected(:,14)-(track_selected(:,13))));
            sd=std(abs(track_selected(:,14)-(track_selected(:,13))));
            
            
            [n1,~]=find(abs(track_selected(:,14)-(track_selected(:,13)))<avg+3*sd);
            [n2,~]=find(abs(track_selected(:,14)-(track_selected(:,13)))<avg+2*sd);
            [n3,~]=find(abs(track_selected(:,14)-(track_selected(:,13)))<avg+sd);
            
            
            
            [o,~]=find(abs(track_selected(:,13)-(track_selected(:,14)))>=3);
            
            [m2,~]=find(abs(track_selected(:,4)-(track_selected(:,9)))<0.9);
            [m1,~]=find(abs(track_selected(:,14)-(track_selected(:,13)))<0.9);
            
            %             [m22,~]=find(abs(track_selected(:,4)-(track_selected(:,9)))>=0.2);
            %             [m11,~]=find(abs(track_selected(:,14)-(track_selected(:,13)))>=0.2);
            
            if ~isempty(m2)
                
                figure(f1)
                subplot(2,2,1)
                scatter(track_selected(m2,2),track_selected(m2,1),[],(track_selected(m2,4)-(track_selected(m2,9)))*100,'filled');
                colormap(bluewhitered(10));
                
                c=colorbar;
                hold on
                plot(lon2,lat2,'-k')
                %                 plot(track_selected(m22,2),track_selected(m22,1),'.r');
                
                
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                
                title(strcat(str,' vs Nemo-Nordic : ',datestr(ConvertSerialYearToDate(min(track_selected(:,3))),'dd mmm yyyy'),' - ',...
                    datestr(ConvertSerialYearToDate(max(track_selected(:,3))),'dd mmm yyyy')),'FontSize',14,'FontWeight','bold');
                
                ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
                c.Label.String = '[cm]';
                box on
                pbaspect([1,0.986038663089498,0.986038663089498]);
                set(gca,'fontname','Times New Roman');
                
                
                hold off
                
                
                subplot(2,2,2)
                scatter(track_selected(m1,2),track_selected(m1,1),[],(track_selected(m1,14)-(track_selected(m1,13)))*100,'filled');
                colormap(bluewhitered(10));
                
                c=colorbar;
                hold on
                plot(lon2,lat2,'-k')
                %                 plot(track_selected(m11,2),track_selected(m11,1),'.r')
                
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                
                title(strcat(str,' vs Nemo-Nordic : ',datestr(ConvertSerialYearToDate(min(track_selected(:,3))),'dd mmm yyyy'),' - ',...
                    datestr(ConvertSerialYearToDate(max(track_selected(:,3))),'dd mmm yyyy')),'FontSize',14,'FontWeight','bold');
                
                ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
                c.Label.String = '[cm]';
                box on
                pbaspect([1,0.986038663089498,0.986038663089498]);
                set(gca,'fontname','Times New Roman');
                set(gca,'yticklabel',[])

                hold off
                
                %                 subplot(2,2,3)
                %                  plot(track_selected(o,2),track_selected(o,1),'+r');
                %                 hold on
                % plot(lon2,lat2,'-k')
                %
                %                 hold off
                %                subplot(2,2,4)
                %
                %                 plot(track_selected(:,2),track_selected(:,1),'.k');
                %                 hold on
                %                                 plot(track_selected(n1,2),track_selected(n1,1),'or');
                %                                 plot(track_selected(n2,2),track_selected(n2,1),'oy');
                %                                 plot(track_selected(n3,2),track_selected(n3,1),'og');
                %
                %
                % plot(lon2,lat2,'-k')
                %
                %                 ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                %
                %                 title(strcat(str,' vs Nemo-Nordic : ',datestr(ConvertSerialYearToDate(min(track_selected(:,3))),'dd mmm yyyy'),' - ',...
                %                     datestr(ConvertSerialYearToDate(max(track_selected(:,3))),'dd mmm yyyy')),'FontSize',14,'FontWeight','bold');
                %
                %                 ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                %                 xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
                %
                %                 hold off
                %
                %
                
                
                
                mi=mi+dt;
                clear m1
                clear m2
                clear m22
                clear m11
                clear n
                clear n1
                clear n2
                clear n3
                clear track_selected
                
            end
            
            
            clear c
            clear ax
            
            f1=f1+1;
            %             f2=f2+1;
            
        end
    end
    clear track_all
    clear mi
    
end


%% plot bias
% 
% for i=1:length(tg_s3b)
%     subplot (3,1,1)
% plot(S3B.date(S3B.pass==tg_s3b(i,1) & abs(S3B.bias)<=0.6),smooth(S3B.bias(S3B.pass==tg_s3b(i,1) & abs(S3B.bias)<=0.6),'rlowess'),'-')
% hold on
% end
col=jet(length(tg_s3a));
for i=1:length(tg_s3a)
    subplot (2,1,1)
    if ~isempty(S3A.date(S3A.pass==tg_s3a(i,1)))
plo(i)=plot(S3A.date(S3A.pass==tg_s3a(i,1)),smooth(S3A.bias(S3A.pass==tg_s3a(i,1)),'rlowess'),'-','LineWidth',1.5,'DisplayName',num2str(tg_s3a(i,1)),'Color',col(i,:));
hold on
    end
end
ylabel('Bias_p_a_s_s [cm]','FontSize',12);
title(strcat('SA-HDM bias'),'FontSize',12,'FontWeight','bold')
h=legend('show');
h.NumColumns=20;
datetick('x','mmm yy','keepticks')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
pbaspect([1,0.25,1]);
hold off
clear i
clear plo
clear col

col=jet(length(tg_ja3));
for i=1:length(tg_ja3)
    subplot (2,1,2)
     if ~isempty(JA3.date(JA3.pass==tg_ja3(i,1)))
plo(i)=plot(JA3.date(JA3.pass==tg_ja3(i,1)),smooth(JA3.bias(JA3.pass==tg_ja3(i,1)),'rlowess'),'-','LineWidth',1.5,'DisplayName',num2str(tg_s3a(i,1)),'Color',col(i,:));
     end
hold on
end

ylabel('Bias_p_a_s_s [cm]','FontSize',12);
h=legend('show');
h.NumColumns=20;
datetick('x','mmm yy','keepticks')
clear i
clear plo    
 ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
 pbaspect([1,0.25,1]);
 box on

 hold off

 
 %% HDM bias
 
S = shaperead('BS.shp');
f=1;
l=1;

for k=1:length(track_all_cor_cell)
    %     for i=1:length(track_all_cor_cell)
    i=2;
    %         for j=1:length(track_all_cor_cell{k}{i})
    j=4;
    figure(f)
    scatter(track_all_cor_cell{k}{i}{j}(:,2),track_all_cor_cell{k}{i}{j}(:,1),25,(track_all_cor_cell{k}{i}{j}(:,13)-(track_all_cor_cell{k}{i}{j}(:,9))),'filled');
    c=colorbar;
    hold on
    l=l+1
end
%     end
% end



plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast
                
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
                
                title('Nemo-Nordic correction','FontSize',14,'FontWeight','bold');
                
                ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
                xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
                c.Label.String = '[m]';
                
%                 hold off


%% corrected file comparision
    
% load res files and convert to table
s3a=load('S3A.mat');
s3a=s3a.track_all_cor_cell;


for k=1:length(s3a)
    for j=1:length(s3a{k})
        if j==1 && k==1
            sa_all=cell2mat(s3a{k}{j});
        else
            sa_all=[sa_all;cell2mat(s3a{k}{j})];
        end
    end
end
% [m,~]=find(isnan(sa_all(:,9)));
% sa_all(m,:)=[];


time=sa_all(:,3);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,3))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,3)));
year=t_sa(:,1);
month=t_sa(:,2);

for i=1:length(t_sa)
    if t_sa(i,2)==1 || t_sa(i,2)==2 || t_sa(i,2)==3
season(i,1)=1;
    elseif  t_sa(i,2)==4 || t_sa(i,2)==5 || t_sa(i,2)==6
season(i,1)=2;
    elseif  t_sa(i,2)==7 || t_sa(i,2)==8 || t_sa(i,2)==9
season(i,1)=3;
    else
        season(i,1)=4;
    end
end
cycle=sa_all(:,5);
pass=sa_all(:,6);
lat=sa_all(:,1);
lon=sa_all(:,2);
sa_orig=sa_all(:,4);
hdm_orig=sa_all(:,9);
N=sa_all(:,9)-sa_all(:,7);
sa=sa_all(:,14);
hdm=sa_all(:,13);

S3A = table(lat,lon,date,time,year,month,season,cycle,pass,N,sa_orig,hdm_orig,sa,hdm);

clear k
clear m
clear date
clear time
clear year
clear month
clear season
clear cycle
clear pass
clear sa_orig
clear hdm_orig
clear hdm
clear sa
clear lat
clear lon
clear i
clear t_sa
clear sa_all
clear j
clear N
clear s3a

% load res files and convert to table
s3b=load('S3B.mat');
s3b=s3b.track_all_cor_cell;


for k=1:length(s3b)
    for j=1:length(s3b{k})
        if j==1 && k==1
            sa_all=cell2mat(s3b{k}{j});
        else
            sa_all=[sa_all;cell2mat(s3b{k}{j})];
        end
    end
end
% [m,~]=find(isnan(sa_all(:,9)));
% sa_all(m,:)=[];


time=sa_all(:,3);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,3))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,3)));
year=t_sa(:,1);
month=t_sa(:,2);

for i=1:length(t_sa)
    if t_sa(i,2)==1 || t_sa(i,2)==2 || t_sa(i,2)==3
season(i,1)=1;
    elseif  t_sa(i,2)==4 || t_sa(i,2)==5 || t_sa(i,2)==6
season(i,1)=2;
    elseif  t_sa(i,2)==7 || t_sa(i,2)==8 || t_sa(i,2)==9
season(i,1)=3;
    else
        season(i,1)=4;
    end
end
cycle=sa_all(:,5);
pass=sa_all(:,6);
lat=sa_all(:,1);
lon=sa_all(:,2);
sa_orig=sa_all(:,4);
hdm_orig=sa_all(:,9);
N=sa_all(:,9)-sa_all(:,7);
sa=sa_all(:,14);
hdm=sa_all(:,13);

S3B = table(lat,lon,date,time,year,month,season,cycle,pass,N,sa_orig,hdm_orig,sa,hdm);

clear k
clear m
clear date
clear time
clear year
clear month
clear season
clear cycle
clear pass
clear sa_orig
clear hdm_orig
clear hdm
clear sa
clear lat
clear lon
clear i
clear t_sa
clear sa_all
clear j
clear N
clear s3b



% load res files and convert to table
ja3=load('JA3.mat');
ja3=ja3.track_all_cor_cell;


for k=1:length(ja3)
    for j=1:length(ja3{k})
        if j==1 && k==1
            sa_all=cell2mat(ja3{k}{j});
        else
            sa_all=[sa_all;cell2mat(ja3{k}{j})];
        end
    end
end
% [m,~]=find(isnan(sa_all(:,9)));
% sa_all(m,:)=[];


time=sa_all(:,3);
date=datetime(datestr(ConvertSerialYearToDate(sa_all(:,3))));
t_sa=datevec(ConvertSerialYearToDate(sa_all(:,3)));
year=t_sa(:,1);
month=t_sa(:,2);

for i=1:length(t_sa)
    if t_sa(i,2)==1 || t_sa(i,2)==2 || t_sa(i,2)==3
season(i,1)=1;
    elseif  t_sa(i,2)==4 || t_sa(i,2)==5 || t_sa(i,2)==6
season(i,1)=2;
    elseif  t_sa(i,2)==7 || t_sa(i,2)==8 || t_sa(i,2)==9
season(i,1)=3;
    else
        season(i,1)=4;
    end
end
cycle=sa_all(:,5);
pass=sa_all(:,6);
lat=sa_all(:,1);
lon=sa_all(:,2);
sa_orig=sa_all(:,4);
hdm_orig=sa_all(:,9);
N=sa_all(:,9)-sa_all(:,7);
sa=sa_all(:,14);
hdm=sa_all(:,13);

JA3 = table(lat,lon,date,time,year,month,season,cycle,pass,N,sa_orig,hdm_orig,sa,hdm);

clear k
clear m
clear date
clear time
clear year
clear month
clear season
clear cycle
clear pass
clear sa_orig
clear hdm_orig
clear hdm
clear sa
clear lat
clear lon
clear i
clear t_sa
clear sa_all
clear j
clear N
clear ja3

%% 
ye=2019;
for i=1:12
z=(JA3.hdm_orig(~isnan(JA3.N) & JA3.year==ye & JA3.month==i)-JA3.hdm(~isnan(JA3.N) & JA3.year==ye & JA3.month==i))*100;
[m,~]=find(isnan(z));
z(m,:)=[];
re_ja3(i,:)=[mean(z),rms(z),std(z)];
clear z
clear m
end

for i=1:12
z=(S3A.hdm_orig(~isnan(S3A.N) & S3A.year==ye & S3A.month==i)-S3A.hdm(~isnan(S3A.N) & S3A.year==ye & S3A.month==i))*100;
[m,~]=find(isnan(z));
z(m,:)=[];
re_s3a(i,:)=[mean(z),rms(z),std(z)];
clear z
clear m
end


for i=1:12
z=(S3B.hdm_orig(~isnan(S3B.N) & S3B.year==ye & S3B.month==i)-S3B.hdm(~isnan(S3B.N) & S3B.year==ye & S3B.month==i))*100;
[m,~]=find(isnan(z));
z(m,:)=[];
re_s3b(i,:)=[mean(z),rms(z),std(z)];
clear z
clear m
end



% plot
load('LATLON.mat')
S = shaperead('BS.shp');



i=2;
ye=2019;

x1=JA3.lat(~isnan(JA3.N) & JA3.year==ye & JA3.month==i);
y1=JA3.lon(~isnan(JA3.N) & JA3.year==ye & JA3.month==i);
zzz1=(JA3.hdm_orig(~isnan(JA3.N) & JA3.year==ye & JA3.month==i)-JA3.hdm(~isnan(JA3.N) & JA3.year==ye & JA3.month==i))*100;
zz1=(JA3.hdm_orig(~isnan(JA3.N) & JA3.year==ye & JA3.month==i));
z1=JA3.hdm(~isnan(JA3.N) & JA3.year==ye & JA3.month==i);

[m1,~]=find(isnan(z1));
x1(m1,:)=[];
y1(m1,:)=[];
z1(m1,:)=[];
zz1(m1,:)=[];
zzz1(m1,:)=[];


x2=S3A.lat(~isnan(S3A.N) & S3A.year==ye & S3A.month==i);
y2=S3A.lon(~isnan(S3A.N) & S3A.year==ye & S3A.month==i);
zzz2=(S3A.hdm_orig(~isnan(S3A.N) & S3A.year==ye & S3A.month==i)-S3A.hdm(~isnan(S3A.N) & S3A.year==ye & S3A.month==i))*100;
zz2=(S3A.hdm_orig(~isnan(S3A.N) & S3A.year==ye & S3A.month==i));
z2=S3A.hdm(~isnan(S3A.N) & S3A.year==ye & S3A.month==i);

[m2,~]=find(isnan(z2));
x2(m2,:)=[];
y2(m2,:)=[];
z2(m2,:)=[];
zz2(m2,:)=[];
zzz2(m2,:)=[];


x3=S3B.lat(~isnan(S3B.N) & S3B.year==ye & S3B.month==i);
y3=S3B.lon(~isnan(S3B.N) & S3B.year==ye & S3B.month==i);
zzz3=(S3B.hdm_orig(~isnan(S3B.N) & S3B.year==ye & S3B.month==i)-S3B.hdm(~isnan(S3B.N) & S3B.year==ye & S3B.month==i))*100;
zz3=(S3B.hdm_orig(~isnan(S3B.N) & S3B.year==ye & S3B.month==i));
z3=S3B.hdm(~isnan(S3B.N) & S3B.year==ye & S3B.month==i);

[m3,~]=find(isnan(z3));
x3(m3,:)=[];
y3(m3,:)=[];
z3(m3,:)=[];
zz3(m3,:)=[];
zzz3(m3,:)=[];


x=[x1;x2;x3];
y=[y1;y2;y3];
z=[z1;z2;z3];
zz=[zz1;zz2;zz3];
zzz=[zzz1;zzz2;zzz3];


F1 = scatteredInterpolant(x,y,z,'nearest');
F2 = scatteredInterpolant(x,y,zz,'nearest');
F3 = scatteredInterpolant(x,y,zzz,'nearest');
% nearest

Z=nan(728,733);
ZZ=nan(728,733);
ZZZ=nan(728,733);


for i=1:728
    for j=1:733
        if ~isnan(double(Lat(i,j)))
            Z(i,j)=F1(double(Lat(i,j)),double(Lon(i,j)));
        end
    end
end

for i=1:728
    for j=1:733
        if ~isnan(double(Lat(i,j)))
            ZZ(i,j)=F2(double(Lat(i,j)),double(Lon(i,j)));
        end
    end
end

for i=1:728
    for j=1:733
        if ~isnan(double(Lat(i,j)))
            ZZZ(i,j)=F3(double(Lat(i,j)),double(Lon(i,j)));
        end
    end
end

figure(1)
subplot(1,3,1)
contourf(Lon,Lat,ZZ,15)
hold on
c=colorbar; 
pbaspect([1 .91 .91])
plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast           
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
title('Nemo-Nordic Original Surface','FontSize',14,'FontWeight','bold');
ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
c.Label.String = '[m]';
hold off

subplot(1,3,2)
[~,M] =contourf(Lon,Lat,ZZZ,20);
hold on
c=colorbar; 
set(M,'LineColor','none')
pbaspect([1 .91 .91])
plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast           
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
title('Correction Surface','FontSize',14,'FontWeight','bold');
% ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
c.Label.String = '[cm]';
hold off


subplot(1,3,3)
contourf(Lon,Lat,ZZ,15)
hold on
c=colorbar; 
pbaspect([1 .91 .91])
plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast           
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold'; grid on
title('Nemo-Nordic Corrected Surface','FontSize',14,'FontWeight','bold');
% ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
c.Label.String = '[m]';
hold off


figure(1)
[~,M] =contourf(Lon,Lat,ZZZ,20);
hold on
c=colorbar; 
set(M,'LineColor','none')
pbaspect([1 .91 .91])
plot(S.X,S.Y,'Color', [0.7 0.7 0.7]); % plot the coast           
ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; %grid on;  ax.GridAlpha = 0.3;
title('Correction Surface','FontSize',14,'FontWeight','bold');
ylabel('Latitude [°]','FontSize',14,'FontWeight','bold');
xlabel('Longitude [°]','FontSize',14,'FontWeight','bold');
c.Label.String = '[cm]';
caxis([-20 90]);
hold off


%% stats

% table 4: monthly avarage bias_sa
% figure 7: monthly avg bias
splitapply(@mean, S3A.bias(S3A.year==2017),S3A.month(S3A.year==2017))*100
splitapply(@mean, S3A.bias(S3A.year==2018),S3A.month(S3A.year==2018))*100
splitapply(@mean, S3A.bias(S3A.year==2019),S3A.month(S3A.year==2019))*100



splitapply(@mean,  JA3.bias( JA3.year==2017), JA3.month( JA3.year==2017))*100
splitapply(@mean,  JA3.bias( JA3.year==2018), JA3.month( JA3.year==2018))*100
splitapply(@mean,  JA3.bias( JA3.year==2019), JA3.month( JA3.year==2019))*100

splitapply(@mean, S3B.bias(S3B.year==2018),S3B.month(S3B.year==2018))*100
splitapply(@mean, S3B.bias(S3B.year==2019),S3B.month(S3B.year==2019))*100



% table 5,stats
sum(S3A.obs)
sum(S3A.out)
mean(S3A.rmse,'omitnan')*100

sum(JA3.obs)
sum(JA3.out)
mean(JA3.rmse,'omitnan')*100


sum(S3B.obs)
sum(S3B.out)
mean(S3B.rmse,'omitnan')*100


% figure 10: RMSE monthly average 
splitapply(@mean, S3A.rmse(S3A.year==2017),S3A.month(S3A.year==2017))*100
splitapply(@mean, S3A.rmse(S3A.year==2018),S3A.month(S3A.year==2018))*100
splitapply(@mean, S3A.rmse(S3A.year==2019),S3A.month(S3A.year==2019))*100


splitapply(@mean, JA3.rmse(JA3.year==2017),JA3.month(JA3.year==2017))*100
splitapply(@mean, JA3.rmse(JA3.year==2018),JA3.month(JA3.year==2018))*100
splitapply(@mean, JA3.rmse(JA3.year==2019),JA3.month(JA3.year==2019))*100

splitapply(@mean, S3B.rmse(S3B.year==2018),S3B.month(S3B.year==2018))*100
splitapply(@mean, S3B.rmse(S3B.year==2019),S3B.month(S3B.year==2019))*100

%% filter data

% JA3
stat=load('ja3_table.mat');
tableja3=stat.JA3;
% data=load('ja3_data.mat');
% dataja3=data.JA3;
for i=1:height(tableja3)
    [m,~]=find(abs(dataja3.sa(dataja3.pass==tableja3.pass(i)&dataja3.cycle==tableja3.cycle(i))-dataja3.hdm(dataja3.pass==tableja3.pass(i)&dataja3.cycle==tableja3.cycle(i)))>tableja3.sd(i));
    

    dataja3.sa(m)=nan;
    dataja3.hdm(m)=nan;

    a{i}=m;
    clear m

end
    [n1,~]=find(dataja3.sa>100);
    [n2,~]=find(dataja3.hdm>100); 
    dataja3.sa(n1)=nan;
    dataja3.hdm(n1)=nan;
    dataja3.sa(n2)=nan;
    dataja3.hdm(n2)=nan;
    clear n1
    clear n2
    clear i
    clear stat
    clear data
    
%  S3A   
stat=load('s3a_table.mat');
tables3a=stat.S3A;
% data=load('s3a_data.mat');
% datas3a=data.S3A;
for i=1:height(tables3a)
    [m,~]=find(abs(datas3a.sa(datas3a.pass==tables3a.pass(i)&datas3a.cycle==tables3a.cycle(i))-datas3a.hdm(datas3a.pass==tables3a.pass(i)&datas3a.cycle==tables3a.cycle(i)))>tables3a.sd(i));
    

    datas3a.sa(m)=nan;
    datas3a.hdm(m)=nan;

    a{i}=m;
    clear m

end
    [n1,~]=find(datas3a.sa>100);
    [n2,~]=find(datas3a.hdm>100); 
    datas3a.sa(n1)=nan;
    datas3a.hdm(n1)=nan;
    datas3a.sa(n2)=nan;
    datas3a.hdm(n2)=nan;
    clear n1
    clear n2
    clear i   
    clear stat
    clear data
    
% S3B    
stat=load('s3b_table.mat');
tables3b=stat.S3B;
% data=load('s3b_data.mat');
% datas3b=data.S3B;
for i=1:height(tables3b)
    [m,~]=find(abs(datas3b.sa(datas3b.pass==tables3b.pass(i)&datas3b.cycle==tables3b.cycle(i))-datas3b.hdm(datas3b.pass==tables3b.pass(i)&datas3b.cycle==tables3b.cycle(i)))>tables3b.sd(i));
    

    datas3b.sa(m)=nan;
    datas3b.hdm(m)=nan;

    a{i}=m;
    clear m

end
    [n1,~]=find(datas3b.sa>100);
    [n2,~]=find(datas3b.hdm>100); 
    datas3b.sa(n1)=nan;
    datas3b.hdm(n1)=nan;
    datas3b.sa(n2)=nan;
    datas3b.hdm(n2)=nan;
    clear n1
    clear n2
    clear i
	clear stat
    clear data
      
%% sub basin analysis

[lat,lon,z] = read_kml('gof.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));

figure(1)
subplot(3,3,1)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1)+1,max(y1)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+1,max(y1)-1,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+1,max(y1)-2,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+1,max(y1)-3,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,2)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2)+1,max(y2)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+1,max(y2)-1,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+1,max(y2)-2,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+1,max(y2)-3,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,3)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3)+1,max(y3)-1,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+1,max(y3)-2,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+1,max(y3)-3,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+1,max(y3)-4,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


clear lat
clear lon
clear in1
clear in2
clear in3
clear x1
clear x2
clear x3
clear y1
clear y2
clear y3
clear m1
clear m2
clear m3



[lat,lon,z] = read_kml('gor.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));


subplot(3,3,4)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1),max(y1)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-0.5,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-1,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-1.5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,5)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2),max(y2)-0,sprintf('y = %0.2fx + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-0.5,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-1,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-1.5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,6)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3),max(y3)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-0.5,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-1,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-1.5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o



clear lat
clear lon
clear in1
clear in2
clear in3
clear x1
clear x2
clear x3
clear y1
clear y2
clear y3
clear m1
clear m2
clear m3


[lat,lon,z] = read_kml('baltic.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));


subplot(3,3,7)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1)+0.2,max(y1)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-4,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-6,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-8,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,8)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2)+0.2,max(y2)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-4,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-6,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-8,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,9)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3)+0.2,max(y3)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-4,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-6,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-8,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o



%figure2


[lat,lon,z] = read_kml('br.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));

figure(2)
subplot(3,3,1)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1),max(y1)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-1,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-2,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-3,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,2)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2),max(y2)-0,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-1,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-2,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-3,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,3)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3),max(y3)-1,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-2,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-3,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-4,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


clear lat
clear lon
clear in1
clear in2
clear in3
clear x1
clear x2
clear x3
clear y1
clear y2
clear y3
clear m1
clear m2
clear m3



[lat,lon,z] = read_kml('bb.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));


subplot(3,3,4)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1),max(y1)-1,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-2,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-3,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1),max(y1)-4,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,5)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2),max(y2)-1,sprintf('y = %0.2fx + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-2,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-3,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2),max(y2)-4,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,6)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3),max(y3)-1,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-2,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-3,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3),max(y3)-4,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
% xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o



clear lat
clear lon
clear in1
clear in2
clear in3
clear x1
clear x2
clear x3
clear y1
clear y2
clear y3
clear m1
clear m2
clear m3


[lat,lon,z] = read_kml('arkona.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);


x1=dataja3.hdm(in1);
y1=dataja3.sa(in1);
x2=datas3a.hdm(in2);
y2=datas3a.sa(in2);
x3=datas3b.hdm(in3);
y3=datas3b.sa(in3);


[m1,~]=find(isnan(x1));
x1(m1,:)=[];
y1(m1,:)=[];

[m2,~]=find(isnan(x2));
x2(m2,:)=[];
y2(m2,:)=[];

[m3,~]=find(isnan(x3));
x3(m3,:)=[];
y3(m3,:)=[];




[n,~]=find(abs(x1-y1)>3);
x1(n,:)=[];
y1(n,:)=[];

[m,~]=find(abs(x1-y1)>3*std(x1-y1));
[o,~]=find(abs(x1-y1)<3*std(x1-y1));


subplot(3,3,7)
plot(x1(m),y1(m),'.r')
hold on
plot(x1(o),y1(o),'+','Color',[0.6350 0.0780 0.1840])

b = [ones(size(x1,1),1) x1]\y1;
RegressionLine = [ones(size(x1,1),1) x1]*b;
plot(x1,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y1-mean(y1)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x1)+0.2,max(y1)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-3,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-4,strcat('STD= ',num2str(std(x1-y1)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x1)+0.2,max(y1)-5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x1))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman')



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o


[n,~]=find(abs(x2-y2)>3);
x2(n,:)=[];
y2(n,:)=[];

[m,~]=find(abs(x2-y2)>3*std(x2-y2));
[o,~]=find(abs(x2-y2)<3*std(x2-y2));

subplot(3,3,8)
plot(x2(m),y2(m),'.','Color',[0.3010 0.7450 0.9330])
hold on
plot(x2(o),y2(o),'+','Color',[0 0.4470 0.7410])
b = [ones(size(x2,1),1) x2]\y2;
RegressionLine = [ones(size(x2,1),1) x2]*b;
plot(x2,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y2-mean(y2)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y2-mean(y2)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x2)+0.2,max(y2)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-3,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-4,strcat('STD= ',num2str(std(x2-y2)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x2)+0.2,max(y2)-5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x2))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');



clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o

[n,~]=find(abs(x3-y3)>3);
x3(n,:)=[];
y3(n,:)=[];

[m,~]=find(abs(x3-y3)>3*std(x3-y3));
[o,~]=find(abs(x3-y3)<3*std(x3-y3));


subplot(3,3,9)
plot(x3(m),y3(m),'.g')
hold on
plot(x3(o),y3(o),'+','Color',[0.4660 0.6740 0.1880])
b = [ones(size(x3,1),1) x3]\y3;
RegressionLine = [ones(size(x3,1),1) x3]*b;
plot(x3,RegressionLine,'k','LineWidth',2.5,'displayname',sprintf('Regression line (y = %0.2fx + %0.2f)',b(2),b(1)))
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y3-mean(y3)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y3-mean(y3)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
text(min(x3)+0.2,max(y3)-2,sprintf('y = %0.2f*x + %0.2f',b(2),b(1)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-3,strcat('R^{2}= ','  ',num2str(R_squared)),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-4,strcat('STD= ',num2str(std(x3-y3)*100,3),' cm '),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
text(min(x3)+0.2,max(y3)-5,strcat('No. Out= ',num2str(length(m)),'(',num2str((length(m)/length(x3))*100,2),'%)'),'FontName','Times New Roman','FontWeight','Bold','FontSize',12)
% ylabel('SSH_S_A','FontSize',14,'FontWeight','bold');
xlabel('SSH_H_D_M','FontSize',14,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=14; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');


clear SS_X
clear SS_Y
clear SS_XY
clear b
clear R_squared
clear RegressionLine
clear m
clear o
%% plot subbasin location

[lat,lon,z] = read_kml('gof.kml'); %1 %4
geoplot([lat(1) lat(4)],[lon(1) lon(4)],'k-','LineWidth',4)
hold on
[lat1,lon1,z] = read_kml('gor.kml'); %3 %4
geoplot([lat1(4) lat(4)],[lon1(4) lon(4)],'k-','LineWidth',4)


% [lat2,lon2,z] = read_kml('br.kml'); %2 %3
geoplot([lat(1) 59.80],[lon(1) 19],'k-','LineWidth',4)

[lat,lon,z] = read_kml('baltic.kml'); %13 %2 %3
geoplot([lat(2) lat(3)],[lon(2) lon(3)],'k-','LineWidth',4)
geoplot([lat(2) lat(13)],[lon(2) lon(13)],'k-','LineWidth',4)

[lat,lon,z] = read_kml('bb.kml'); %1 %2
geoplot([63.12 63.82],[21.663 20.47],'k-','LineWidth',4)



ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');

geobasemap bluegreen

clear lat
clear lon

%% remove NaN data from datble 

dataja3=dataja3(~any(ismissing(dataja3),2),:);
datas3a=datas3a(~any(ismissing(datas3a),2),:);
datas3b=datas3b(~any(ismissing(datas3b),2),:);

%% dist2coast
% [lat,lon,z] = read_kml('coast2.kml');
% coast=[lat lon];


load('bs.mat')
% 
% lat=bs(:,1);
% lon=bs(:,2);

[k,~] = dsearchn(bs,[dataja3.lat dataja3.lon]);
for i=1:length(k)
dis(i,:)=distance([dataja3.lat(i) dataja3.lon(i)],[bs(k(i),1),bs(k(i),2)],referenceEllipsoid('WGS84'))/1000; %measure distance to the bs
i
end
dataja3= [dataja3 table(dis, 'VariableNames', {'dist2'})];
clear k
clear dis

[k,~] = dsearchn(bs,[datas3a.lat datas3a.lon]);
for i=1:length(k)
dis(i,:)=distance([datas3a.lat(i) datas3a.lon(i)],[bs(k(i),1),bs(k(i),2)],referenceEllipsoid('WGS84'))/1000; %measure distance to the bs
i
end
datas3a= [datas3a table(dis, 'VariableNames', {'dist2'})];
clear k
clear dis

[k,~] = dsearchn(bs,[datas3b.lat datas3b.lon]);
for i=1:length(k)
dis(i,:)=distance([datas3b.lat(i) datas3b.lon(i)],[bs(k(i),1),bs(k(i),2)],referenceEllipsoid('WGS84'))/1000; %measure distance to the bs
i
end
datas3b= [datas3b table(dis, 'VariableNames', {'dist2'})];
clear k
clear dis

clear lat
clear lon
clear z
%% plot histogram
%complete later

x1=[dataja3.hdm-dataja3.sa dataja3.hdm_orig-dataja3.sa_orig dataja3.dist];
[m1,~]=find(isnan(x1(:,1)));
x1(m1,:)=[];
[n1,~]=find(abs(x1(:,1))>3);
x1(n1,:)=[];

subplot(2,3,1)
h2=histogram(x1(:,2));
hold on
h1=histogram(x1(:,1));

h1.FaceColor = 'r';
h1.EdgeColor  = 'r';
h2.FaceColor = 'b';
h2.EdgeColor  = 'b';
pbaspect([1 1 1])

ylabel('No. Data','FontSize',18,'FontWeight','bold');
xlabel('\DeltaSSH[m]','FontSize',18,'FontWeight','bold');

ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
box on
xlim([-1 3])

subplot(2,3,4)
plot(x1(:,3),x1(:,2),'.b')
hold on
plot(x1(:,3),x1(:,1),'.r')
pbaspect([1 1 1])


clear x1

x1=[datas3a.hdm-datas3a.sa datas3a.hdm_orig-datas3a.sa_orig datas3a.dist];
[m1,~]=find(isnan(x1(:,1)));
x1(m1,:)=[];
[n1,~]=find(abs(x1(:,1))>3);
x1(n1,:)=[];
subplot(2,3,2)
h2=histogram(x1(:,2));
hold on
h1=histogram(x1(:,1));
h1.FaceColor = 'r';
h1.EdgeColor  = 'r';
h2.FaceColor = 'b';
h2.EdgeColor  = 'b';
pbaspect([1 1 1])
ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
box on
% ylabel('No. Data','FontSize',18,'FontWeight','bold');
xlabel('\DeltaSSH[m]','FontSize',18,'FontWeight','bold');
xlim([-1 3])


subplot(2,3,5)
plot(x1(:,3),x1(:,2),'.b')
hold on
plot(x1(:,3),x1(:,1),'.r')
pbaspect([1 1 1])


clear x1


x1=[datas3b.hdm-datas3b.sa datas3b.hdm_orig-datas3b.sa_orig datas3b.dist];
[m1,~]=find(isnan(x1(:,1)));
x1(m1,:)=[];
[n1,~]=find(abs(x1(:,1))>3);
x1(n1,:)=[];
subplot(2,3,3)
h2=histogram(x1(:,2));
hold on
h1=histogram(x1(:,1));
h1.FaceColor = 'r';
h1.EdgeColor  = 'r';
h2.FaceColor = 'b';
h2.EdgeColor  = 'b';
pbaspect([1 1 1])
ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
box on
% ylabel('No. Data','FontSize',18,'FontWeight','bold');
xlabel('\DeltaSSH[m]','FontSize',18,'FontWeight','bold');
xlim([-1 3])


subplot(2,3,6)
plot(x1(:,3),x1(:,2),'.b')
hold on
plot(x1(:,3),x1(:,1),'.r')
pbaspect([1 1 1])
clear x1







%% plot dssh-dist2coast

[lat,lon,z] = read_kml('baltic.kml');
in1=inpolygon(dataja3.lat,dataja3.lon,lat,lon);
in2=inpolygon(datas3a.lat,datas3a.lon,lat,lon);
in3=inpolygon(datas3b.lat,datas3b.lon,lat,lon);

 x1=[dataja3.sa(in1)-dataja3.hdm(in1),dataja3.dist(in1)];
 x2=[datas3a.sa(in2)-datas3a.hdm(in2),datas3a.dist(in2)];
 x3=[datas3b.sa(in3)-datas3b.hdm(in3),datas3b.dist(in3)];
% 
% 
% 
% x1=[(dataja3.sa(in1)-dataja3.N(in1))-(dataja3.hdm(in1)-dataja3.N(in1)),dataja3.dist(in1)];
% x2=[(datas3a.sa(in2)-datas3a.N(in2))-(datas3a.hdm(in2)-datas3a.N(in2)),datas3a.dist(in2)];
% x3=[(datas3b.sa(in3)-datas3b.N(in3))-(datas3b.hdm(in3)-datas3b.N(in3)),datas3b.dist(in3)];

dis1(:,1)=min(x1(:,2))-0.5:0.5:max(x1(:,2))+0.5;
dis2(:,1)=min(x2(:,2))-0.5:0.5:max(x2(:,2))+0.5;
dis3(:,1)=min(x3(:,2))-0.5:0.5:max(x3(:,2))+0.5;


 [n1,~]=find(abs(x1(:,1))>3);
 [n2,~]=find(abs(x2(:,1))>3);
 [n3,~]=find(abs(x3(:,1))>3);

 
 
%  [n1,~]=find(abs((dataja3.sa-dataja3.N)-(dataja3.hdm-dataja3.N))<1);
%  [n2,~]=find(abs((datas3a.sa-datas3a.N)-(datas3a.hdm-datas3a.N))<1);
%  [n3,~]=find(abs((datas3b.sa-datas3b.N)-(datas3b.hdm-datas3b.N))<1);
% 
% % figure(1)
% plot(dataja3.dist(n1,:),(dataja3.sa(n1,:)-dataja3.hdm(n1,:)),'.','LineWidth',2.5,'Color',[0.6350 0.0780 0.1840]);
% hold on
% plot(datas3a.dist(n2,:),(datas3a.sa(n2,:)-datas3a.hdm(n2,:)),'.','LineWidth',2.5,'Color',[0.3010 0.7450 0.9330]);
% plot(datas3b.dist(n3,:),(datas3b.sa(n3,:)-datas3b.hdm(n3,:)),'.','LineWidth',2.5,'Color',[0.4660 0.6740 0.1880]);


% 
% figure(1)
% plot(dataja3.dist(n1,:),((dataja3.sa(n1,:)-dataja3.N(n1,:))-(dataja3.hdm(n1,:)-dataja3.N(n1,:)))*100,'.','LineWidth',2.5,'Color',[0.6350 0.0780 0.1840]);
% hold on
% plot(datas3a.dist(n2,:),((datas3a.sa(n2,:)-datas3a.N(n2,:))-(datas3a.hdm(n2,:)-datas3a.N(n2,:)))*100,'.','LineWidth',2.5,'Color',[0.3010 0.7450 0.9330]);
% plot(datas3b.dist(n3,:),((datas3b.sa(n3,:)-datas3b.N(n3,:))-(datas3b.hdm(n3,:)-datas3b.N(n3,:)))*100,'.','LineWidth',2.5,'Color',[0.4660 0.6740 0.1880]);
% 
% 
% figure(1)
% plot(dataja3.dist,(dataja3.sa-dataja3.hdm),'.','LineWidth',2.5,'Color',[0.6350 0.0780 0.1840]);
%  hold on
%  plot(datas3a.dist,(datas3a.sa-datas3a.hdm),'.','LineWidth',2.5,'Color',[0.3010 0.7450 0.9330]);
%  plot(datas3b.dist,(datas3b.sa-datas3b.hdm),'.','LineWidth',2.5,'Color',[0.4660 0.6740 0.1880]);



x1(n1,:)=[];
x2(n2,:)=[];
x3(n3,:)=[];


%JA3
for i=1:length(dis1)
[m1,~]=find(x1(:,2)>=dis1(i)&x1(:,2)<dis1(i)+0.5);
rmse1(i,1)=rms(x1(m1,1))*100;
end

%S3A
for i=1:length(dis2)
[m2,~]=find(x2(:,2)>=dis2(i)&x2(:,2)<dis2(i)+0.5);
rmse2(i,1)=rms(x2(m2,1))*100;  
end

%S3B
for i=1:length(dis3)
[m3,~]=find(x3(:,2)>=dis3(i)&x3(:,2)<dis3(i)+0.5);
rmse3(i,1)=rms(x3(m3,1))*100;  
end


figure(2)
plot(dis1,rmse1,'-','LineWidth',2.5,'Color',[0.6350 0.0780 0.1840])
hold on
plot(dis2,rmse2,'-','LineWidth',2.5,'Color',[0.3010 0.7450 0.9330])
plot(dis3,rmse3,'-','LineWidth',2.5,'Color',[0.4660 0.6740 0.1880])

yline(std([x1(:,1);x2(:,1);x3(:,1)])*100,'--','STD','LineWidth',3 ,'FontName','Times New Roman','FontWeight','Bold','FontSize',18);
xline(10,':r','LineWidth',3);


patch('XData',[0 0 60 60],'YData',[0 std([x1(:,1);x2(:,1);x3(:,1)])*100  std([x1(:,1);x2(:,1);x3(:,1)])*100 0],'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.3)

xlabel('Lat','FontSize',18,'FontWeight','bold');
ylabel('Diff_S_A_H_D_M [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman');
% legend({'JA3','S3A','S3B'})
xlim([0 50])
pbaspect([1,0.125,0.125])

clear dis1
clear dis2
clear dis3
clear rmse1
clear rmse2
clear rmse3
clear x1
clear x2
clear x3
clear in1
clear in2
clear in3



% 
% d=splitapply(@std, (dataja3.hdm-dataja3.sa),dataja3.dist(dataja3.dist>1.5));


clear n1
clear n2
clear n3

