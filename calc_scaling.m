%Returns the best fit relaxation time scaling
%freqs are frequencies
%ref is reference harmonics
%tar is target
function relaxtime=calc_scaling(freqs,ref,tar)
scales=linspace(0.25,1.75,150);
errs=zeros(length(scales),1);

for j=1:length(scales)
    scaledfreqs=scales(j)*freqs;
    indices=scaledfreqs>=min(freqs) & scaledfreqs<=max(freqs); %ignore outer points, spline less reliable there;
    %indices=1:5; %Just use all
    %indices=1:4; %skip last one
    indices=1;
    errs(j)=norm(tar(indices)-spline(freqs,ref,scaledfreqs(indices)),'fro')^2/numel(tar(indices));
end
[val,k]=min(errs(:));
['Best RMS:' num2str(sqrt(val))]
relaxtime=scales(k)
end