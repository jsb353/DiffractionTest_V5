%% Get the XRD Pattern
% Run through StructureFactor with the same parameters except the
% Lattice.Normal to obtain multiple 2theta angles.
% diffraction geometry of the (102) peak (noncoplanar)
% It needs the Lattice and Probe to get information for StructureFactor,
% Threshold to limit the peaks to only the important ones (it is
% intensity), Range to define how close two peaks can be together before
% they are interpreted as a single one from the same family of planes and
% span/IndexMax for what Miller indeces should be analyzed 
% e.g. span=1 100 010 001 110 101 110 111 
% e.g. span=2 100 010 001 110 101 110 111 102 120 201 210 012 021 and so
% on.
% 
%
%
% Last updated 4-12-2017 Cosmin Popescu

% function [PlotIntensity,Lattice, Probe]=Generate_Intensity_2theta(Lattice, Probe,Threshold, Range, IndexMax,FigNum)
%% Get the XRD Pattern
% Run through StructureFactor with the same parameters except the
% Lattice.Normal to obtain multiple 2theta angles.
% diffraction geometry of the (102) peak (noncoplanar)
% Last updated 4-10-2017 Cosmin Popescu

addpath(genpath('C:\Users\Cosmin\Desktop\Grand Diffraction Master'))
addpath(genpath('C:\Users\Cosmin\Desktop\Cr2AlC'))
addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\StructureLibrary'))
addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\TestScripts'))
%addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master'))
%clear;
% Load your material
%load graphite.mat

% DEFINE YOUR X-RAYS
%Probe.Type = 'x-ray';
%Probe.Energy = 8047; % [eV] % define either Energy or lambda
Probe.Polarization = 's'; % s (perpendicular) or p (parallel)

Threshold=0.001;
Range=0.1;
IndexMax=9;
FigNum=[];

%% Loop over different hkl values individually. 
% There will be repetitions
% like 001 and 002
DataOfIntensity=zeros(10,4);
countofdata=1;


for h=0:IndexMax
    for k=0:IndexMax
        for l=0:IndexMax
            if h==0 && k==0 && l==0
                % This is to eliminate the 000
            else
            Lattice.Reflection=[h k l];
            DataOfIntensity(countofdata,1)=h;
            DataOfIntensity(countofdata,2)=k;
            DataOfIntensity(countofdata,3)=l;
            [Result,Lattice,Probe]= StructureFactor(Lattice,Probe);
            DataOfIntensity(countofdata,4)=Result.Intensity;
            DataOfIntensity(countofdata,5)=2*Result.BraggAngle;
            countofdata=countofdata+1;
            end
        end
    end
end

%% Eliminate the values with intensity below a threshold.
countofplot=1;
PlotIntensity=zeros(1,5);
for i=1:countofdata-1
    if(DataOfIntensity(i,4)>Threshold && isreal(DataOfIntensity(i,5)))
        PlotIntensity(countofplot,:)=DataOfIntensity(i,:);
        countofplot=countofplot+1;
    end
end
%% Sort everything by the Bragg angle.
sorted=false;
while sorted==false
    n=1;
   while n<countofplot-1
        if PlotIntensity(n,5)>PlotIntensity(n+1,5)
            
            store=PlotIntensity(n+1,:);
            PlotIntensity(n+1,:)=PlotIntensity(n,:);
            PlotIntensity(n,:)=store;
            n=1;
        else 
            n=n+1;
        end
   
   end
   sorted=true;
    
end
    

%% Make Plot
if isempty(FigNum)
    figure
else
    figure(FigNum)
end
plot(PlotIntensity(:,5),PlotIntensity(:,4),'ob','linewidth',2);
xlabel('2\theta');
ylabel('Intensity');
Title=strcat(Lattice.Symbol,'-',Lattice.Type);
title(Title);
axis([0 180 0 1.1*max(PlotIntensity(:,4))]);
hold on;

%% Add the Miller indeces to the plot
xseparation=2;
yseparation=max(PlotIntensity(:,4))*0.04;
% This is a way to separate multiple combinations of hkl that give the same
% theta. 
PlotIntensity(:,6)=ones(countofplot-1,1);
for i=2:countofplot-1
    if PlotIntensity(i,5)==PlotIntensity(i-1,5)
        PlotIntensity(i,6)=PlotIntensity(i-1,6)+1;
    end
end
% Write text at a separation from the point. If multiple point in the same
% location, stack the hkl in order defined by previous for.
for i=1:countofplot-1
    h=PlotIntensity(i,1);
    k=PlotIntensity(i,2);
    l=PlotIntensity(i,3);
    text( PlotIntensity(i,5)+xseparation,PlotIntensity(i,4)+yseparation*PlotIntensity(i,6),strcat([num2str(h), num2str(k), num2str(l)]))

end
%% Intensities on the same angle are added up

PlotIntensity(:,7)=zeros(countofplot-1,1);
PlotIntensity(:,8)=PlotIntensity(:,4);
[length,~]=size(PlotIntensity);
for i=1:length
    n=i;
    currentposition=i;
    if(n<length-1)
        while(PlotIntensity(n+1,6)>1)
            PlotIntensity(currentposition,8)=PlotIntensity(currentposition,8);%+PlotIntensity(n+1,8);
            PlotIntensity(n+1,8)=0;
            if(n+1<length)
            n=n+1;
            else
                break
            end
        end
    else
        if(PlotIntensity(n,6)==1)
            PlotIntensity(n,4)=PlotIntensity(n,4);
        end
    end
    
end

%% Generating XRD plot with resolution based on range.
if isempty(FigNum)
    figure
else
    figure(FigNum+1)
end
countAngle=1;
XRD_plot=zeros(ceil(180/Range+length),5);
indexPlot=1;
for angle=0:Range:180
    if (angle<PlotIntensity(countAngle,5))&& (angle+Range>PlotIntensity(countAngle,5)&&countAngle<length-1)
        XRD_plot(indexPlot,1)=PlotIntensity(countAngle,1); %h
        XRD_plot(indexPlot,2)=PlotIntensity(countAngle,2); %k
        XRD_plot(indexPlot,3)=PlotIntensity(countAngle,3); %l
        XRD_plot(indexPlot,4)=PlotIntensity(countAngle,5); %2theta
        XRD_plot(indexPlot,5)=PlotIntensity(countAngle,8); %Intensity of the actual point
        indexPlot=indexPlot+1;
        if countAngle<length-1
            while(PlotIntensity(countAngle+1,6)>1&&(countAngle<(length-2)))
                countAngle=countAngle+1;
            end
        end
        countAngle=countAngle+1;
    else
        XRD_plot(indexPlot,4)=angle; %Looping through the angles.
        XRD_plot(indexPlot,5)=0; %No intensity
        indexPlot=indexPlot+1;
        
    end       
end
plot(XRD_plot(:,4),XRD_plot(:,5),'-b');
xlabel('2\theta');
ylabel('Intensity');
Title=strcat(Lattice.Symbol,'-',Lattice.Type);
title(Title);
axis([0 180 0 1.1*max(XRD_plot(:,5))]);
hold on;
maxx=max(XRD_plot(:,5));
XRD_plot(:,5) = 100*XRD_plot(:,5)/maxx;

%% A new plot similar to actual XRD is generated.
% Xtheta=zeros(5,180+length);
% count=1;
% % Make 0 in intensity everything that is not in the data.
% for i=1:180+length
%     if count<=length
%         if((PlotIntensity(count,5)-i+count)<Range)
%             Xtheta(1,i)=PlotIntensity(count,5);
%             Xtheta(2,i)=PlotIntensity(count,8);
%             Xtheta(3,i)=PlotIntensity(count,1);
%             Xtheta(4,i)=PlotIntensity(count,2);
%             Xtheta(5,i)=PlotIntensity(count,3);
%             count=count+1;
%         else
%             Xtheta(1,i)=i-count;
%         end
%     else
%         Xtheta(1,i)=i-count;
%         
%     end
% end
% Xtheta=Xtheta';
% figure;
% plot(Xtheta(:,1),Xtheta(:,2),'-b');
% xlabel('2\theta');
% ylabel('Intensity');
% title(Lattice.Symbol);
% axis([0 180 0 1.1*max(Xtheta(:,2))]);
% hold on;
%% Make name for the text file

    Filename=strcat(Lattice.Symbol,'_',Lattice.Type,'.txt');


%% Make table with Miller indeces.
fid=fopen(Filename,'a');
fprintf(fid,'\r\n');
fprintf(fid,'Start of a new data set \r\n');
fprintf(fid,strcat('The energy is ',num2str(Probe.Energy)));
fprintf(fid,'\r\n %s %s %s %s %s\r\n','h','k','l','2theta','Intensity');
for i=1:length
    if(PlotIntensity(i,6)==1)
        fprintf (fid,'%d %d %d %f %f\r\n',PlotIntensity(i,1),PlotIntensity(i,2),PlotIntensity(i,3),PlotIntensity(i,5),PlotIntensity(i,8));
        
%     else
%         fprintf(fid,'%d %d %d\r\n',PlotIntensity(i,1),PlotIntensity(i,2),PlotIntensity(i,3));
        
    end
end
fclose(fid);

