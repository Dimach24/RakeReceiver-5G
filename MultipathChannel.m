classdef MultipathChannel<handle
    methods(Static)
        function multipath_params = generateRaysParams(rays_count, ...
                min_rays_path_diff, max_rays_path_diff, ...
                samples_rate, ...
                mean_attenuation,...
                do_phase_distortion)
            % Generates uniformly distributed offsets, Rayleigh distributed gains, uniformly distributed phase shifts (if needed)
            arguments(Input)
                % totally rays (main ray included)
                rays_count (1, 1) {mustBeInteger, mustBePositive};
                % minimal distance in meters between main and reflected ray
                min_rays_path_diff (1, 1) {mustBePositive};
                % maximal distance in meters between main and reflected ray
                max_rays_path_diff (1, 1) {mustBePositive};
                % samples per 1 second
                samples_rate (1, 1) {mustBeInteger, mustBePositive};
                % mean value for Rayleigh distributed amplitude
                % attenuation (gain) factors for reflected rays
                mean_attenuation (1, 1) {mustBePositive};
                % does reflection provoke a uniform phase shift for reflected rays
                do_phase_distortion (1, 1) {mustBeNonempty} = true;
            end
            arguments(Output)
                % structure with fields {offsets, gains} which are rays_count-1 length
                multipath_params (1,1)
            end
            % mean(Rayleigh)=sqrt(pi/2)*std(Gauss), where Rayleigh=sqrt(GaussX^2+GaussY^2)
            % n_std = std(Gauss); std of normally distributed GaussX and GaussY
            n_std=mean_attenuation*sqrt(2/pi);
            multipath_params.gains=randn(1, rays_count-1)*n_std+1j*randn(1, rays_count-1)*n_std;
            if ~ do_phase_distortion
                multipath_params.gains=abs(multipath_params.gains);
            end
            max_time_delay=max_rays_path_diff/3e8; % 3e8 = speed of light
            max_samples_offset=floor(max_time_delay*samples_rate)
            min_time_delay=min_rays_path_diff/3e8;
            min_samples_offset=floor(min_time_delay*samples_rate)
            multipath_params.offsets=randi([min_samples_offset+1 max_samples_offset], 1, rays_count-1);
        end
        function signal=generateMultipathSignal(samples, multipath_params)
            rays=length(multipath_params.gains);
            N=length(samples);
            copies=zeros(rays, N);
            for i=1:rays
                copies(i, multipath_params.offsets(i):end)=...
                    samples(1:N-multipath_params.offsets(i)+1)*multipath_params.gains(i);
            end
            signal=sum(copies);
        end
    end
end