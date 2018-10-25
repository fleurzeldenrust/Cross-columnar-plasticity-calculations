% Opens all the .mat files in folder and exports the input,membrane potential and spikeindices in seperate csv files
% Does not check structs for availability of data
% 30-04-2018, Linda Wouters

% SET DIRECTORY TO CORRECT FOLDER!!!
% Importing files based on https://nl.mathworks.com/matlabcentral/answers/78626-how-do-i-import-multiple-files-from-a-folder-in-matlab#answer_88362
testfiledir = 'C:\Users\lab\Documents\Python\exported_files';
matfiles = dir(fullfile(testfiledir, '*.mat'));
min_spikes = 40
unconverted = {}

for i = 1 : length(matfiles)
   disp(i)
   disp(matfiles(i).name)

   % Load struct
   load(matfiles(i).name);

   % Set index to first inh cell with >40 spikes if it exists,
   % else, set it to first exc cell with >40 spikes
   index = 0;
   for j = 1 : length(Data)
       disp(length(Data{j}.spikeindices))
       if Data{j}.input_generation_settings.tau == 50 && length(Data{j}.spikeindices) > min_spikes
           index = j;
           break
       elseif index == 0 && length(Data{j}.spikeindices) > min_spikes
           index = j;
       end
   end
   
   % If an index was found, convert it, else append filename to
   % 'unconverted'
   if index > 0
       disp(length(Data{j}.spikeindices))
       % Set outputfile base-name
       base_name = erase(matfiles(i).name, ".mat");

       % Add exc or inh to filename
       if Data{index}.input_generation_settings.tau == 250
           disp("exc")
           base_name = strcat("exc_", int2str(index), "_", base_name)
           disp(base_name)
       elseif Data{index}.input_generation_settings.tau == 50
           disp("inh")
           base_name = strcat("inh_", int2str(index), "_", base_name)
           disp(base_name)
       end

       % Export input, membrane potential and spikeindices
       csvwrite(strcat("15_", base_name, "_input.csv"), Data{index}.input_current)
       % membrane_potential
       csvwrite(strcat("15_", base_name, "_membrane.csv"), Data{index}.membrane_potential)
       % % spikeindices
       csvwrite(strcat("15_", base_name, "_spikeindices.csv"), Data{index}.spikeindices)
   else
       % Code from
       % https://stackoverflow.com/questions/2288899/appending-string-to-matlab-array,
       % answered Feb 18 '10 at 14:06, Amro
       unconverted{end+1} = matfiles(i).name;
   end
end

% Save the files that are not converted
if length(unconverted) > 0
    disp(unconverted)
    save('unconverted_files', 'unconverted')
end
