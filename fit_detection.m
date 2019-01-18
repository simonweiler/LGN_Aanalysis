 t=1000:10000
 y=test;
al=(y(t) > y(t-1)) & (y(t) > y(t+1))

dt=500;
m=0.8
al2=(y(t) - y(t-dt) > m) & (y(t) - y(t+dt) > m)