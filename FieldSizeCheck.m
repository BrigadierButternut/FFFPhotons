%% Field Size check


%Importing raw data
raw_data = importdata('fs6x6_6X_10FFF.xlsx');

%6 MV profile at different depths; field size = 2x2
FF6MV = raw_data.data.x6X;

%10 MV FFF profile at different depths; field size = 2x2
FFF10MV = raw_data.data.x10FFF;

center_points = FF6MV(2:end,2);

Depth = cell(12,1);

FF_Field_Size = zeros(12,1);

FFF_Field_Size = zeros(12,1);

Difference = zeros(12,1);

loop = 1;

DoseData = 3;

depthvar = 1; 

while loop <= 12
    
    %Dose Profile for the 6 MV beam
    FFProfile = FF6MV(2:end,DoseData);

    %Dose Profile for the 10 MV beam
    FFFProfile = FFF10MV(2:end,DoseData);
    
    xaxis = zeros(length(center_points)+1,1);
    
    
    FF_zero_dose = interp1(center_points,FFProfile,0);
    
    FFF_zero_dose = interp1(center_points,FFFProfile,0);
    
    FFProfile = FFProfile/FF_zero_dose;
    
    FFFProfile = FFFProfile/FFF_zero_dose;
    
    
    [coarseFF,old_FF_index,new_FF_index] = unique(FFProfile);
    
    [coarseFFF,old_FFF_index,new_FFF_index] = unique(FFFProfile);
    
    FF_PDD50 = interp1(coarseFF,center_points(old_FF_index),0.5);
    
    FFF_PDD50 = interp1(coarseFFF,center_points(old_FFF_index),0.5);
    
    FF_field_size = abs(FF_PDD50*2);
    
    FFF_field_size = abs(FFF_PDD50*2);
    
    field_size_diff = abs(FF_PDD50-FFF_PDD50);
    
    
    Depth(loop)=raw_data.textdata.x10FFF(1,DoseData-2);
    
    FF_Field_Size(loop) = FF_field_size;
    
    FFF_Field_Size(loop) = FFF_field_size;
    
    Difference(loop) = field_size_diff;
    
    
    DoseData = DoseData +5;
    
    loop = loop + 1;
    
end

Ouput = table(Depth,FF_Field_Size,FFF_Field_Size,Difference)
  
    
    
    
    