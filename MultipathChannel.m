classdef MultipathChannel
    methods(Static)
        function multipath_params = generateRaysParams(rays_count,max_rays_path_diff,sample_rate,mean_attenuation,do_phase_distortion)
            multipath_params.gains
            multipath_params.offsets
        end
        function signal=generateMultipathSignal(samples,multipath_params)
        end
    end
end