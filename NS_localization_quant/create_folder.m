if not(isfolder([selectedimage_folder]))
    mkdir([selectedimage_folder])
end

if not(isfolder([selectedimage_folder filesep 'cell']))
    mkdir([selectedimage_folder filesep 'cell'])
end

if not(isfolder([selectedimage_folder filesep 'nucleus']))
    mkdir([selectedimage_folder filesep 'nucleus'])
end

if not(isfolder([selectedimage_folder filesep 'untrans']))
    mkdir([selectedimage_folder filesep 'untrans'])
end

if not(isfolder([selectedimage_folder filesep 'untrans' filesep 'spec']))
    mkdir([selectedimage_folder filesep 'untrans' filesep 'spec'])
end

if not(isfolder([selectedimage_folder filesep 'trans']))
    mkdir([selectedimage_folder filesep 'trans'])
end

if not(isfolder([selectedimage_folder filesep 'trans' filesep 'spec']))
    mkdir([selectedimage_folder filesep 'trans' filesep 'spec'])
end

% if not(isfolder([selectedimage_folder filesep 'trans' filesep 'RNA']))
%     mkdir([selectedimage_folder filesep 'trans' filesep 'RNA'])
% end

% if not(isfolder([selectedimage_folder filesep 'Analysis']))
%     mkdir([selectedimage_folder filesep 'Analysis'])
% end
