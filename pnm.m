function result=pnm(n,m)
result=(-1)^((m+n-1)/2)*(2*m+1)*(2*n+1)*factd(m-1)*factd(n)/((m-n)*(m+n+1)*factd(m)*factd(n-1));