
%% Get the XRD Pattern
% Run through StructureFactor with the same parameters except the
% Lattice.Normal to obtain multiple 2theta angles.
% diffraction geometry of the (102) peak (noncoplanar)
% It needs the Lattice and Probe to get information for StructureFactor,
% Threshold to limit the peaks to only the important ones (it is
% intensity), Resolution to define how close two peaks can be together before
% they are interpreted as a single one from the same family of planes and
% span/IndexMax for what planes should be analyzed
%
%
%
% Last updated 5-8-2017 Cosmin Popescu

function [Table]=Generate_Intensity_2theta(Lattice, Probe,Threshold, Resolution, hkl,FigNum)
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

% % Load your material
% load Cr2AlC.mat
%
% % DEFINE YOUR X-RAYS
% Probe.Type = 'x-ray';
% Probe.Energy = 8047; % [eV] % define either Energy or lambda
% Probe.Polarization = 's'; % s (perpendicular) or p (parallel)
%
% Threshold=1;
% Resolution=0.1;
% IndexMax=9;
% FigNum=[];

%% Loop over different hkl values individually.
% There will be repetitions
% like 001 and 002
% The software goes through all combinations given by user of Miller
% indeces and returns the resulting information from Structure Factor. If
% the Bragg angle is real and the intensity is above a threshold, these
% values are save.
TypeOfFile=input('What kind of file extension do you want?\nOptions: nothing (press Enter), txt, xlsx, xls, dat, csv \n','s');

MainData=zeros(10,4);
countofdata=1;


for h=hkl:-1:-hkl
    for k=hkl:-1:-hkl
        for l=hkl:-1:-hkl
            if h==0 && k==0 && l==0
                % This is to eliminate the 000
            else
                
                Lattice.Reflection=[h k l];
                [Result,Lattice,Probe]= StructureFactor(Lattice,Probe);
                if Result.Intensity>Threshold && isreal(Result.BraggAngle)
                    
                    MainData(countofdata,1)=h;
                    MainData(countofdata,2)=k;
                    MainData(countofdata,3)=l;
                    
                    MainData(countofdata,4)=Result.Intensity;
                    MainData(countofdata,5)=2*Result.BraggAngle;
                    MainData(countofdata,6)=Result.Distance;
                    countofdata=countofdata+1;
                end
            end
        end
    end
end


length=size(MainData,1);
%% Sort everything by the Bragg angle.
% It uses a bubble sort to rearrange all rows based on the fifth column
% which is the two theta.
sorted=false;
while sorted==false
    n=1;
    while n<length
        if MainData(n,5)>MainData(n+1,5)
            % If it finds a place where it is not sorted, it rearranges the
            % two rows and restarts going through all lines. It needs
            % changed to speed up the process because this is the slowest
            % algorithm for sorting and the only one I had in mind then.
            store=MainData(n+1,:);
            MainData(n+1,:)=MainData(n,:);
            MainData(n,:)=store;
            n=1;
        else
            n=n+1;
        end
    end
    sorted=true;
end
% SortedData=zeros(size(MainData,1),size(MainData,2));
% SortedData=TopDownMergeSort(MainData,SortedData,size(MainData,2),5);
% MainData=SortedData;

%% Make Plot

if nargin>5
    figure(FigNum)
    plot(MainData(:,5),MainData(:,4),'ob','linewidth',2);
    xlabel('2\theta');
    ylabel('Intensity');
    Title=strcat(Lattice.Symbol,'-',Lattice.Type);
    title(Title);
    axis([0 180 0 1.1*max(MainData(:,4))]);
    hold on;
end

%% Add the Miller indeces to the plot
xseparation=-2;
yseparation=max(MainData(:,4))*0.05;
% This is a way to separate multiple combinations of hkl that give the same
% theta.
MainData(:,7)=ones(length,1);
for i=2:length
    if MainData(i,5)==MainData(i-1,5)
        MainData(i,7)=MainData(i-1,7)+1;
    end
end
% Write text at a separation from the point. If multiple point in the same
% location, stack the hkl in order defined by previous for.
if nargin >5
    for i=1:length
        if MainData(i,7)==1
            h=MainData(i,1);
            k=MainData(i,2);
            l=MainData(i,3);
            text( MainData(i,5)+xseparation,MainData(i,4)+yseparation,strcat([num2str(h), num2str(k), num2str(l)]))
        end
    end
end
%% Intensities on the same angle are added up


MainData(:,8)=MainData(:,4);

for i=1:length
    n=i;
    currentposition=i;
    if(n<length-1)
        while(MainData(n+1,7)>1)
            MainData(currentposition,8)=MainData(currentposition,8)+MainData(n+1,8);
            MainData(n+1,8)=0;
            if(n+1<length)
                n=n+1;
            else
                break
            end
        end
    else
        if(MainData(n,7)==1)
            MainData(n,8)=MainData(n,4);
        end
    end
    
end

%% Generating XRD plot with resolution based on range.
if nargin>5
    figure(FigNum+1)
else
    figure
end
countAngle=1;
XRD_plot=zeros(ceil(180/Resolution+length),5);
indexPlot=1;
ArrayIndex=1;
maxI=max(MainData(:,8));
for angle=0:Resolution:180
    if (angle<MainData(countAngle,5))&& (angle+Resolution>MainData(countAngle,5)&&countAngle<length-1)
        Table.h=MainData(countAngle,1); %h
        Table.k=MainData(countAngle,2); %k
        Table.l=MainData(countAngle,3); %l
        Table.BraggAngle=MainData(countAngle,5)/2;
        Table.TwoTheta=MainData(countAngle,5);
        Table.Intensity=MainData(countAngle,8);
        Table.RelativeIntensity=Table.Intensity/maxI*100;
        Table.d=MainData(countAngle,6);
        Table.TwoPi_Distance=2*pi/Table.d;
        Table.d_r=1/Table.d;
        
        XRD_plot(indexPlot,1)=MainData(countAngle,5); %2theta
        XRD_plot(indexPlot,2)=MainData(countAngle,8); %Intensity of the actual point
        indexPlot=indexPlot+1;
        if countAngle<length-1
            while(MainData(countAngle+1,7)>1&&(countAngle<(length-2)))
                countAngle=countAngle+1;
            end
        end
        Table.m=MainData(countAngle,7);
        Array(ArrayIndex)=Table;
        ArrayIndex=ArrayIndex+1;
        countAngle=countAngle+1;
        
    else
        XRD_plot(indexPlot,1)=angle; %Looping through the angles.
        XRD_plot(indexPlot,2)=0; %No intensity
        indexPlot=indexPlot+1;
        
    end
end
% subplot(2,1,1);
hold on;
plot(XRD_plot(:,1),XRD_plot(:,2)*100/max(XRD_plot(:,2)),'-b','linewidth',1.5);
xlabel('2\theta');
ylabel('Intensity percent I/I_M_a_x');
Title=strcat(Lattice.Symbol,'-',Lattice.Type,'- ',num2str(Probe.Energy),'eV');
title(Title);
axis([0 180 0 105]);
hold on;

legend('MATLAB');
%% A new plot similar to actual XRD is generated.
% Xtheta=zeros(5,180+length);
% count=1;
% % Make 0 in intensity everything that is not in the data.
% for i=1:180+length
%     if count<=length
%         if((MainData(count,5)-i+count)<Resolution)
%             Xtheta(1,i)=MainData(count,5);
%             Xtheta(2,i)=MainData(count,8);
%             Xtheta(3,i)=MainData(count,1);
%             Xtheta(4,i)=MainData(count,2);
%             Xtheta(5,i)=MainData(count,3);
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

%% Make table with Miller indeces.
Table=struct2table(Array);
time=datestr(datetime('now'));
time(15)='_';
time(18)='_';
time(12)='_';
NameOFile=strcat(Lattice.Symbol,'_',Lattice.Type,'_',num2str(Probe.Energy),'eV','_',time,'.');

if size(TypeOfFile,2)==4
    if TypeOfFile=='xlsx'
        writetable(Table,strcat(NameOFile,'xlsx'));
    end
elseif size(TypeOfFile,2)==3
    if TypeOfFile== 'txt'
        writetable(Table,strcat(NameOFile,'txt'));
    elseif TypeOfFile=='dat'
        writetable(Table,strcat(NameOFile,'dat'));
    elseif TypeOfFile == 'csv'
        writetable(Table,strcat(NameOFile,'csv'));
    elseif TypeOfFile == 'xls'
        writetable(Table,strcat(NameOFile,'xls'));
    end
else %default do nothing
end




% Former text file generation on the basis of the data generated. Updated
% version of data output above in table generation with different
% extension.
% fid=fopen(Filename,'a');
% fprintf(fid,'\r\n');
% fprintf(fid,'Start of a new data set \r\n');
% fprintf(fid,strcat('The energy is  ',num2str(Probe.Energy)));
% fprintf(fid,'\r\n %s %s %s %s %s\r\n','h','k','l','2theta','Intensity');
% for i=1:length
%     if(MainData(i,6)==1)
%         fprintf (fid,'%d %d %d %f %f\r\n',MainData(i,1),MainData(i,2),MainData(i,3),MainData(i,5),MainData(i,8));
%
% %     else
% %         fprintf(fid,'%d %d %d\r\n',MainData(i,1),MainData(i,2),MainData(i,3));
%
%     end
% end
% fclose(fid);

