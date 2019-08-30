M=importdata('alchohal_dataset.csv');
%hilbert(M.data(1:255,2));
for k= 1:30
    phase = hilbert(M.data(28*k+1:28*k+25,2))
end
% phase1= atan2(imag(phase), real(phase));
% plot(1:255,phase1)

phase1 = atan2(imag(hilbert(M.voltage.V1(1,15615k+1:15615k+256))), real(hilbert(M.voltage.V1(1,15615k+1:15615k+256))));