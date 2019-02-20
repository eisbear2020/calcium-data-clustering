classdef chirp < handle
    % includes all attributes and function for chirp stimulus
    %
    % TO DOS: 
    %      - include reading time_step, time_freq, time_amp from log file
    %
    %
    
    properties
        t  = []%duration of stimulus
        nframes= []
        freq = []%Hz
        red_intensity= []
        stim_start_times = []
        stim_end_times= []
        time_step = 3;  % Duration of bright step (3s)
        time_freq = 8; % Duration of first sine (8s)
        time_amp = 8 ;  % Duration of second sine (8s)
    end
    
    methods
        function obj = chirp()
            %CHIRP Construct an instance of this class
            %   Detailed explanation goes here
        end
             
            
        
        
        function readLogFile(obj,fname)
            
            % read the file with stimuli parameters
            if ~ isfile(fname)
                display('File not found')
            end
        
            f=fopen(fname);
            tline=fgetl(f);

            if ~contains(tline,'Full_field_chirp')
                disp('The file does not contain the full field chirp stimulus');
                fclose(f);
                return;
            end

            stim_start_times_in = [];
            stim_end_times_in = [];
            formatIn='yyyy/MM/dd HH:mm:ss:SSS';

            while ~feof(f)
                tline = fgetl(f);
                if contains(tline,'Start time: ')
                    tstr=tline(length('Start time: ')+1:end);
                    start = datetime(tstr,'InputFormat',formatIn);
                    stim_start_times_in = [stim_start_times_in,start];
                elseif contains(tline,'End time: ')
                    tstr=tline(length('End time: ')+1:end);
                    end_t = datetime(tstr,'InputFormat',formatIn);
                    stim_end_times_in = [stim_end_times_in,end_t];   
                elseif contains(tline,'Freq: ')
                    freq_str=tline(length('Freq: ')+1:end);
                    obj.freq=str2num(freq_str); %flip rate : frames per second
                elseif contains(tline,'t	 ')
                    t_str=tline(length('t	 ')+1:end);     
                    obj.t=str2num(t_str); %flip rate : frames per second
                elseif contains(tline,'nframes	 ')
                    nframes_str=tline(length('nframes	 ')+1:end);
                    obj.nframes=str2num(nframes_str);              
                elseif contains(tline,'maxRedIntensity	 ')
                    maxRedIntensity_str=tline(length('maxRedIntensity	 ')+1:end);
                    obj.red_intensity=str2num(maxRedIntensity_str);                                                         
                end
                
                obj.stim_start_times = stim_start_times_in;
                obj.stim_end_times = stim_end_times_in;
            end

            
            fclose(f);

        end
       
            
        function chirp = returnScaledChirp(obj,scaled_length)
            % returnScaledChirp re-constructs chirp scaled to length
            %  'scaled_length'. If parameter is not provided the original
            % stimulus is returned
           
            % Bright step
            step = ones(1, obj.time_step*1000);

            % Sinusoidal intensity modulation 
            % A sinusoid curve can be described as a sine function with 4 parameters: 
            % y(time) = Amplitude*sin(2*pi * d[Frequency(Time)] + Phase)

            % With increasing frequency:
            A1 = 0.5;
            t1 = 0:0.001:obj.time_freq;  % length of the wave.
            f1 = 0.5*t1;     % function of t, as it changes over time.
            ph1 = 0.5;
            sine1 = A1.*sin(2*pi*f1.*t1) + ph1;
            % With increasing amplitude
            t2 = 0:0.001:obj.time_amp;
            A2 = t2/(2*obj.time_amp);  % function of t, as it changes over time.
            f2 = -2;
            ph2 = 0.5;

            sine2 = -A2.*sin(2*pi*f2.*t2) + ph2;

            % Separation between stimuli

            delay1 = zeros(1, 2000);  % 2 s black
            delay2 = zeros(1, 3000);  % 3 s black

            delay3 = ones(1, 2000) * 0.5;  % 2 s grey
            delay4 = ones(1, 2000) * 0.5;  % 2 s grey
            delay5 = ones(1, 2000) * 0.5;  % 2 s grey

            delay6 = zeros(1, 2000);  % 2 s black

            % Full field 'chirp' signal
            chirp = [delay1 step delay2 delay3 sine1 delay4 sine2 delay5 delay6];
            nr_frames = round(obj.t/obj.freq);
            chirp_rec = zeros(1,nr_frames);

            for i=1:nr_frames
                chirp_rec(i) = chirp(round((i)*1000/60));
            end
            
            if ~exist('scaled_length','var')
            % third parameter does not exist, so default it to something
                % make last element zero for nicer visualization
                chirp_rec(end) = 0;
                chirp = chirp_rec;
            else  
                chirp = resample(chirp,scaled_length,length(chirp));
            end
        end
    end
end

