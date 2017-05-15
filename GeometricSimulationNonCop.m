% Generate a table of values for the angles for the incident and reflected
% beam. Needs commenting.

function Table=GeometricSimulationNonCop(Lattice, Probe, hkl_space)
% Check what kind of file does the user want.

TypeOfFile=input('What kind of file extension do you want? Options: none, txt, xlsx, xls, dat, csv \n','s');

index=1;
for h=hkl_space
    for k=hkl_space
        for l=hkl_space
            if h==0 && k==0 && l==0 
            else
                Lattice.Reflection=[h k l];
                [SF,Lattice, Probe]=StructureFactor(Lattice,Probe);
                SF.AssymAngle = SF.AssymAngle
               % if isreal(SF.BraggAngle)==1 && SF.Intensity>10^-10
                 if abs(SF.BraggAngle)<=90    && SF.Intensity > 1e-10
                    TEST = [h k l]
                    Result=NonCoplanarDiffraction(SF.BraggAngle, Probe.psi,SF.CrystalNormal,SF.ReflectionNormal );
                    % Make struct variable to store the values from
                    % NonCoplanarDiffraction
                    if isreal(Result.Incident)&& isreal(Result.IncidentSpherical) && ...
                            isreal(Result.Reflected) && isreal(Result.ReflectedSpherical)
                        Table.h=h;
                        Table.k=k;
                        Table.l=l;
                        Table.Incident_Theta=Result.IncidentSpherical(1);
                        Table.Incident_Psi=Result.IncidentSpherical(2);
                        Table.Reflected_Theta=Result.ReflectedSpherical(1);
                        Table.Reflected_Psi=Result.ReflectedSpherical(2);
                        Table.IncidentXYZ=Result.Incident;
                        Table.ReflectedXYZ=Result.Reflected;
                        Table.Intensity = SF.Intensity;
                      %  Table.BraggAngle = Result.BraggAngle;
                      %  Table.AssymAngle = Result.AssymAngle;
                      %  Table.YES = 1;
                      %  if (Result.BraggAngle > Result.AssymAngle)
                      %      Table.YES = 0;
                     %   else
                        %    Table.YES = 0;
                      %  end
                        
                        Array(index)=Table;
                        
                        index=index+1;
                    end
                end
            end
        end
    end
end
Table=struct2table(Array);
if size(TypeOfFile,2)==4
    if TypeOfFile=='xlsx'
        Name=strcat(Lattice.Symbol,'_',Lattice.Type,'_','XRD_spherical_angles.xlsx');
        writetable(Table,Name);
    end
elseif size(TypeOfFile,2)==3
    if TypeOfFile== 'txt'
        Name=strcat(Lattice.Symbol,'_',Lattice.Type,'_','XRD_spherical_angles.txt');
        writetable(Table,Name);
    elseif TypeOfFile=='dat'
        Name=strcat(Lattice.Symbol,'_',Lattice.Type,'_','XRD_spherical_angles.dat');
        writetable(Table,Name);
    elseif TypeOfFile == 'csv'
        Name=strcat(Lattice.Symbol,'_',Lattice.Type,'_','XRD_spherical_angles.csv');
        writetable(Table,Name);
    elseif TypeOfFile == 'xls'
        Name=strcat(Lattice.Symbol,'_',Lattice.Type,'_','XRD_spherical_angles.xls');
        writetable(Table,Name);
    end
else %default
end




