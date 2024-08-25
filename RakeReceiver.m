classdef RakeReceiver
    methods(Static)
        function processed_signal=receiveRaysByCorrPeaks(samples,peaks_shifts,peaks_values)
            % process signal with rake-receive algorithm
            arguments
                samples % signal samples array
                peaks_shifts % shifts to the reflected rays
                peaks_values % values of correlation function at the peaks
            end
            peaks_shifts=peaks_shifts(peaks_shifts>=0);
            fingers=length(peaks_shifts);
            rays=zeros(fingers,length(samples));
            coefs=abs(peaks_values).^2;
            coefs=coefs/sum(coefs);
            phase_correction=peaks_values./abs(peaks_values);
            coefs=coefs./phase_correction;
            for i=1:fingers;
                rays(i,1:length(samples)-peaks_shifts(i)+1)=...
                    samples(peaks_shifts(i):end)*coefs(i);
            end
            processed_signal=sum(rays);
        end
    end
end