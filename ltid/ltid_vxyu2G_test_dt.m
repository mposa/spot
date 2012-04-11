function ltid_vxyu2G_test_dt(n,k,m,N)

if nargin<1, n=5; end
if nargin<2, k=3; end
if nargin<3, m=2; end
if nargin<4, N=100*n; end

G0=ltid_rand(n,k,m);
[A,B,C,D]=ssdata(G0);
U=randn(m,N);
X=randn(2*n,N);
T=1;
V=(A*X+B*U-X)/T+0.01*randn(2*n,N);
Y=C*X+D*U+0.1*randn(k,N);
[G,obj]=ltid_vxyu2G(V,X,Y,U,T,1);
fprintf(' error: obj=%f, HInf=%f out of %f\n', ...
    obj,norm(G-G0,Inf),norm(G0,Inf))
w=linspace(0,pi,200);
for i=1:k,
    for j=1:m,
        g0=squeeze(freqresp(G0(i,j),w));
        g=squeeze(freqresp(G(i,j),w));
        close(gcf)
        subplot(2,1,1);plot(w,real(g0),'.',w,real(g)); grid
        subplot(2,1,2);plot(w,imag(g0),'.',w,imag(g)); grid
        pause
    end
end