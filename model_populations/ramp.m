function s=ramp(x,t,N)

Lt=length(t);
r=ones(1,Lt);
r(1:N)=1/t(N)*t(1:N);
r((Lt-N+1):Lt)=1+(t((Lt-N+1):Lt)-t(Lt-N+1))/(t(Lt-N+1)-t(Lt));

s=r.*x;