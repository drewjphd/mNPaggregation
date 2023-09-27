%Input y(t) should be one period of the signal to be Fourier-transformed
function h=get_harmonics(t,y)
h=zeros(7,1); %Calculate first seven Fourier sine coefficients
dt=t(2)-t(1);
T=t(end)-t(1); %period
for j=1:length(h)
    h(j)=sum(y.*sin(2*pi*t*j/T)*dt*2/T); %Integration (Riemann sum; accurate for small dt)
end

%figure; plot(t,y);
%figure; bar(1:7,h(1:7)); %Show results
end