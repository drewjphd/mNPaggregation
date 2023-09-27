%Returns the best fit relaxation time scaling
%freqs are frequencies
%ref is reference harmonics
%tar is target
function scalings=calc_scaling_v2(reffreqs,ref,tarfreqs,tar)
if length(reffreqs)~=length(ref) | length(tarfreqs)~=length(tar)
    'error'
    return;
end
scalings=zeros(length(tar),1);

for j=1:length(tar)
    scalings(j)=spline(ref,reffreqs,tar(j))/tarfreqs(j);
end
end
