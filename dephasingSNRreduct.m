N=100;
spins = [ones(1,N); zeros(1,N); zeros(1,N)];
phi = (1:N)/N*2*pi*0.9;
for ii=1:N
dephasedSpins(:,ii) = zrot(phi(ii))*spins(:,ii);
end
complexSpins = dephasedSpins(1,:)+i*dephasedSpins(2,:);
abs(sum(complexSpins))/N