function [signal,copies]= multipath_gen(samples,rays_count,offsets,gains)
    copies=zeros(rays_count,length(samples));
    for i=1:rays_count-1
        % p=rand()*2*pi
        copies(i,:)=gains(i)*[zeros(1,offsets(i)) samples(1:end-offsets(i))];%*exp(1j*p);
    end
    signal=sum(copies)+samples;
end