close all
clear all
clc

%% 
% *10 MV FFF versus 6 MV FF*
%% 
% *Field Size: 6x6 cm2*
%% 
% *Vassiliev Normalization*

%Importing raw data
raw_data = importdata('fs6x6_6X_10FFF_corrected.xlsx');

%6 MV profile at different depths; field size = 2x2
FF6MV = raw_data.data.x6X;

%10 MV FFF profile at different depths; field size = 2x2
FFF10MV = raw_data.data.x10FFF;

%Voxel dimensions and locations are the same at all depths; therefore, only
%one column of voxel data is necessary.
voxels = FF6MV(:,1);


off_axis_distance = FF6MV(2:end,2);


loop = 0;

DoseData = 3;

depthvar = 1; 

while loop <= 12
       
%Dose Profile for the 6 MV beam
FFProfile = FF6MV(2:end,DoseData);

%Dose Profile for the 10 MV beam
FFFProfile = FFF10MV(2:end,DoseData);

%The first step is to normalize the 6MV profile; this is done by dividing
%all doses by dose on the central axis.

%Finding the location of the central axis 

%Generating a new x-axis to account for dose at the central axis
xaxis = zeros(length(off_axis_distance)+1,1);

negative_index = find(off_axis_distance<0);

xaxis(negative_index)=off_axis_distance(negative_index);

xaxis((negative_index(end)+2):end) = off_axis_distance((negative_index(end)+1):end);

FFcenter_point_dose = interp1(off_axis_distance,FFProfile,0);

FF_placehold=zeros(length(FFProfile)+1,1);

FF_placehold(negative_index)=FFProfile(negative_index);

FF_placehold(negative_index(end)+1)=FFcenter_point_dose;

FF_placehold((negative_index(end)+2):end)=FFProfile((negative_index(end)+1):end);

FFProfile=zeros(size(FF_placehold));

FFProfile=FF_placehold;

%Doing all that for the FFF profile

FFFcenter_point_dose = interp1(off_axis_distance,FFFProfile,0);

FFF_placehold=zeros(length(FFFProfile)+1,1);

FFF_placehold(negative_index)=FFFProfile(negative_index);

FFF_placehold(negative_index(end)+1)=FFFcenter_point_dose;

FFF_placehold((negative_index(end)+2):end)=FFFProfile((negative_index(end)+1):end);

FFFProfile=zeros(size(FFF_placehold));

FFFProfile=FFF_placehold;

zero_index = negative_index(end)+1;

NormValue = FFProfile(zero_index)/1.05;

%Normalizing the FF profile s.t. the dose on the central axis is 105%
FFProfile = (FFProfile/NormValue);

%Finding where PDD = 100% occurs on the FF profile.
[uniqueFF,old_FF_index,new_FF_index]=unique(FFProfile(negative_index));

%Getting the off-axis distance of the 100% point from the FF profile
PDD100FF = interp1(uniqueFF,xaxis(old_FF_index),1);

%Finding the dose at the location of 100% PDD of the FF profile
[uniqueFFF,old_FFF_index,new_FFF_index]=unique(FFFProfile(negative_index));

NormValueFFF = interp1(xaxis(old_FFF_index),uniqueFFF,PDD100FF);

FFFNormalized = FFFProfile/NormValueFFF;

FFFNormalized = FFFNormalized(:)*100;

FFProfile = FFProfile(:)*100;

Renorm_Point = FFFNormalized(zero_index);

NominalValues = transpose(0:10:100);

%Finding the penumbral width of the FFF Profile

FFFNormalized(isnan(FFFNormalized))=0;

[CoarseFFF,FFFindex,CoarseFFFIndex]=unique(FFFNormalized(1:zero_index));

FFFNormInterp = interp1(CoarseFFF,off_axis_distance(FFFindex),NominalValues);

FFF80 = find(NominalValues == 80);

FFF20 = find(NominalValues == 20);

FFFPenumbraWidth = abs(FFFNormInterp(FFF80)-FFFNormInterp(FFF20));

%Finding the penumbral width of the FF Profile

FFProfile(isnan(FFProfile))=0;

[CoarseFF,FFindex,CoarseFFIndex]=unique(FFProfile(1:zero_index));

FFNormInterp = interp1(CoarseFF,off_axis_distance(FFindex),NominalValues);

FF80 = find(NominalValues == 80);

FF20 = find(NominalValues == 20);

FFPenumbraWidth = abs(FFNormInterp(FFF80)-FFNormInterp(FFF20));

%Comparison of Penumbral Widths

WidthTable = table;

WidthTable.FF_Penumbra_cm = FFPenumbraWidth;

WidthTable.FFF_Penumbra_cm = FFFPenumbraWidth;

%Examining the doses for the FFF beam at 2mm, 5mm, 10mm, 20 mm, and 50 mm from the beam edge

FFF50 = FFFNormInterp(find(NominalValues == 50));

OutofFieldFFF = zeros(5,1);

OutofFieldFFF(1) = FFF50 - 0.2;

OutofFieldFFF(2) = FFF50 - 0.5;

OutofFieldFFF(3) = FFF50 - 1;

OutofFieldFFF(4) = FFF50 - 2;

OutofFieldFFF(5) = FFF50 - 5;

OutofFieldFFFDose = interp1(off_axis_distance(FFFindex(:)),CoarseFFF,OutofFieldFFF);

%Examining the doses for the FF beam at 2mm, 5mm, 10mm, 20 mm, and 50 mm from the beam edge

FF50 = FFNormInterp(find(NominalValues == 50));

OutofFieldFF = zeros(5,1);

OutofFieldFF(1) = FFF50 - 0.2;

OutofFieldFF(2) = FFF50 - 0.5;

OutofFieldFF(3) = FFF50 - 1;

OutofFieldFF(4) = FFF50 - 2;

OutofFieldFF(5) = FFF50 - 5;

OutofFieldFFDose = interp1(off_axis_distance(FFindex(:)),CoarseFF,OutofFieldFF);

%Creating a table to display relative out of field doses

FieldEdgeDistances = [2,5,10,20,50];

DisplayTable = table;

DisplayTable.DistanceFromFieldEdge_mm = FieldEdgeDistances';

DisplayTable.FF_Dose = OutofFieldFFDose;

DisplayTable.FFF_Dose = OutofFieldFFFDose;
%%

depth = raw_data.textdata.x10FFF(1,depthvar);

disp(depth)
disp(DisplayTable)
disp(WidthTable)


%Plotting the normalized profiles
figure
plot(xaxis,FFProfile,xaxis,FFFNormalized,xaxis(zero_index),Renorm_Point,'.','MarkerSize',20,'MarkerEdgeColor','r','MarkerFaceColor','r')
xlabel('Off-Axis Distance (cm)')
ylabel('% Dose')
legend('6 MV FF Profile', '10 MV FFF Profile')
text(xaxis(zero_index),Renorm_Point+4,num2str(Renorm_Point));

depthvar = depthvar + 5; 

DoseData = DoseData + 5; 

loop = loop + 1;

end
