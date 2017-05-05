close all
clear all
clc

%% 10 MV FFF versus 6 MV FF
%% Ponisch Normalization
%% Field Size: 2x2 cm2

%Importing raw data
raw_data = importdata('fs10x10_6X_10FFF.xlsx');

%6 MV profile at different depths; field size = 2x2
FF6MV = raw_data.data.x6X;

%10 MV FFF profile at different depths; field size = 2x2
FFF10MV = raw_data.data.x10FFF;

%Voxel dimensions and locations are the same at all depths; therefore, only
%one column of voxel data is necessary.
voxels = FF6MV(:,1);


off_axis_distance = FF6MV(2:end,2);


loop = 1;

DoseData = 3;

penumbra_table = zeros(12,3);

penumbra_table(1:10,1) = [1:1:10];

penumbra_table(11,1) = 20;

penumbra_table(12,1) = 28;


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

% %Normalizing the FF profile
% FFProfile = (FFProfile(:)/FFcenter_point_dose);

zero_index = negative_index(end)+1;

%First step in Ponisch normalization is to find the inflection points of
%the flattened profile as well as the unflattened profile

g1_FF = gradient(FFProfile);
g2_FF = gradient(g1_FF);

g1_FFF = gradient(FFFProfile);
g2_FFF = gradient(g1_FFF);

%Need to set the second derivative equal to zero to find locations of the
%inflection points

[rough_FF,original_FF_index,brand_new_FF_index] = unique(g2_FF);

[rough_FFF,original_FFF_index,brand_new_FFF_index] = unique(g2_FFF);

inflection_point_FF = interp1(rough_FF,xaxis(original_FF_index),0);

inflection_point_FFF = interp1(rough_FFF,xaxis(original_FFF_index),0);

%Doses at inflection points

[rough_FF,original_FF_index,brand_new_FF_index] = unique(FFProfile);

[rough_FFF,original_FFF_index,brand_new_FFF_index] = unique(FFFProfile);

FF_inflection_dose = interp1(xaxis(original_FF_index),rough_FF,inflection_point_FF);

FFF_inflection_dose = interp1(xaxis(original_FFF_index),rough_FFF,inflection_point_FFF);

%Using the Ponisch normalization equation:

normalized_FFF_dose = (FFF_inflection_dose/FF_inflection_dose)*FFcenter_point_dose;

%The penumbra can then be found as being the distance between the points
%where the dose is equal to 20% of the normalized_FFF_dose and 80% of the
%normalized_FFF_dose

FFF20 = interp1(rough_FFF,xaxis(original_FFF_index),.2*normalized_FFF_dose);

FFF80 = interp1(rough_FFF,xaxis(original_FFF_index),.8*normalized_FFF_dose);

%penumbral width of the FFF beam can now be found

FFF_penumbral_width = abs(FFF80-FFF20);

%Comparing it to the FF beam:

FF20 = interp1(rough_FF,xaxis(original_FF_index),.2*FFcenter_point_dose);

FF80 = interp1(rough_FF,xaxis(original_FF_index),.8*FFcenter_point_dose);

FF_penumbral_width = abs(FF80-FF20);

penumbra_table(loop,2) = FFF_penumbral_width;

penumbra_table(loop,3) = FF_penumbral_width;

DoseData = DoseData + 5; 

loop = loop + 1;

end






