close all
clear all
clc

%% 6 MV FFF versus 4 MV FF 
%% (Vassiliev Normalization)

FFFdata = importdata('6FFF_Beam_Data.xlsx');

FFdata = importdata('4MV_Beam_Data.xlsx');

field_size = 2; 

%Collecting all the field names for FFF and FF data
FFF_field_names = fieldnames(FFFdata.data);

FF_field_names = fieldnames(FFdata.data);

%Finding common field names
field_names = intersect(FFF_field_names,FF_field_names);

depths = strfind(field_names,'OpenFieldProfiles');

field_name_index = find(~cellfun(@isempty,depths));

depth_data = 1;

%This allows me to cycle through all of the different depths that match
%between the two beams
while depth_data <= length(field_name_index)
    
    FFF_Profile_Data = FFFdata.data.(char(field_names(field_name_index(depth_data))));
    
    FF_Profile_Data = FFdata.data.(char(field_names(field_name_index(depth_data))));
    
    %Given a particular depths, it is necessary to cycle through all of the
    %field sizes; that is the purpose of the second loop
    field_size = 2;
    
    
    
   
    
    while field_size <= 6
        
        FFF_Profile = FFF_Profile_Data(:,field_size);
        
        FF_Profile = FF_Profile_Data(:,field_size);
        
        FFFx = FFF_Profile_Data(:,1);

        FFx = FF_Profile_Data(:,1);
        
        %I really only want to match up the data that is aligned for the
        %two profiles; that is the purpose of intersect
        
        [x,FFx_Match_Index,FFFx_Match_Index]=intersect(FFx,FFFx);
        
        FFx=FFx(FFx_Match_Index);
        
        FFFx = FFFx(FFFx_Match_Index);
        
        FF_Profile = FF_Profile(FFx_Match_Index);
        
        FFF_Profile = FFF_Profile(FFFx_Match_Index);
        
        %The first step in Vassiliev's method is to normalize the FF
        %profile to 110% at the center
        
        FF_Profile_Normalized = FF_Profile*1.1;
        
        %The next step is to normalize the FFF beam at 100% on the FF beam
  
        %locating the index of the central axis
        zero_index = find(FFx == 0);
       
        NominalValues = [0:10:100].';
        
        PDD100 = find(NominalValues == 100);
        
        PDD80 = find(NominalValues == 80);
        
        PDD50 = find(NominalValues == 50);
        
        PDD20 = find(NominalValues == 20);
        
        %interp1 doesn't like non-integer values, so I out they go. 
        FF_Profile_Normalized(isnan(FF_Profile_Normalized))=0;
        
        %unique is used to make interp1 happy; it doesn't like duplicates.
        [CoarseFF,FF_index,CoarseFF_index] = unique(FF_Profile_Normalized(1:zero_index));
        
        
        FF_PDD_Nominal = interp1(FF_Profile_Normalized(FF_index),FFx(FF_index),NominalValues);
        
        %Taking out the non-zero values of the FFF beam
        FFF_Profile(isnan(FFF_Profile))=0;
        
        [CoarseFFF, FFF_index, CoarseFFF_index] = unique(FFF_Profile(1:zero_index));
        
        %Normalizing the FFF profile
        RenormPoint = interp1(FFFx(FFF_index),CoarseFFF,FF_PDD_Nominal(PDD100));
        
        FFF_Profile_Normalized = (FFF_Profile/RenormPoint)*100;
        
        [CoarseFFF, FFF_index, CoarseFFF_index] = unique(FFF_Profile_Normalized(1:zero_index));
        
        %And now for finding the penumbral width!
        FFF_PDD_Nominal = interp1(FFF_Profile_Normalized(FFF_index),FFFx(FFF_index),NominalValues);
        
        
        %Finding the penumbral width of the FFF beam
         
        FFF_Penumbral_Width = abs(FFF_PDD_Nominal(PDD80)-FFF_PDD_Nominal(PDD20));
        
        %Finding the penumbral width of the FF beam
       
        FF_Penumbral_Width = abs(FF_PDD_Nominal(PDD80)-FF_PDD_Nominal(PDD20));
        
        %Making a table of the penumbral widths so that they're readable.
        
        WidthTable = table;
        
        WidthTable.FF_Penumbra_cm = FF_Penumbral_Width;
        
        WidthTable.FFF_Penumbra_cm = FFF_Penumbral_Width;
        
        %Checking dose at different distances from the field edge for the
        %FFF beam
        
        OutofField_FFF = zeros(5,1);
        
        OutofField_FFF(1) = FFF_PDD_Nominal(PDD50) - 0.2;
        
        OutofField_FFF(2) = FFF_PDD_Nominal(PDD50) - 0.5;
        
        OutofField_FFF(3) = FFF_PDD_Nominal(PDD50) - 1;
        
        OutofField_FFF(4) = FFF_PDD_Nominal(PDD50) - 3;
        
        OutofField_FFF(5) = FFF_PDD_Nominal(PDD50) - 5;
        
        OutofField_FFF_Dose = interp1(FFFx(FFF_index),CoarseFFF,OutofField_FFF);
        
        %Checking dose at different distances from the field edge for the
        %FF beam
        
        OutofField_FF = zeros(5,1);
        
        OutofField_FF(1) = FF_PDD_Nominal(PDD50) - 0.2;
        
        OutofField_FF(2) = FF_PDD_Nominal(PDD50) - 0.5;
        
        OutofField_FF(3) = FF_PDD_Nominal(PDD50) - 1;
        
        OutofField_FF(4) = FF_PDD_Nominal(PDD50) - 3;
        
        OutofField_FF(5) = FF_PDD_Nominal(PDD50) - 5;
        
        OutofField_FF_Dose = interp1(FFx(FF_index),CoarseFF,OutofField_FF);
        
        %Building a table to compare the out-of-field doses
        
        FieldEdgeDistances = [2;5;10;30;50];
        
        DoseTable = table;
        
        DoseTable.DistanceFromFieldEdge_mm = FieldEdgeDistances;
        
        DoseTable.FF_Relative_Dose = OutofField_FF_Dose;
        
        DoseTable.FFF_Relative_Dose = OutofField_FFF_Dose;
        
        %Plottin' time!
        plot_title2 = char(FFdata.textdata.(char(field_names(field_name_index(depth_data))))(8,field_size));
        
        plot_title1 = char(field_names(field_name_index(depth_data)));
        
        relative_max_dose = FFF_Profile_Normalized(zero_index);
        %% 
        
        disp(char(field_names(field_name_index(depth_data))))
        disp(plot_title2)
        disp(DoseTable)
        disp(WidthTable)
        
        figure
        plot(FFx,FF_Profile_Normalized,FFFx,FFF_Profile_Normalized,FFx(zero_index),relative_max_dose,'.','MarkerSize',20,'MarkerEdgeColor','r','MarkerFaceColor','r')
        ylabel('Relative Dose (%)')
        xlabel('Distance from Central Axis (cm)')
        legend('4 MV FF Beam', '6 MV FFF Beam')
        text(FFx(zero_index),relative_max_dose+4,num2str(relative_max_dose))
        title(strcat(plot_title1,{', '},plot_title2))
        
       
      
        field_size = field_size + 1;
    end
    
    depth_data = depth_data + 1;
end