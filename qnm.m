function result=qnm(n,m)
result=(-1)^((m+n)/2)*(2*m+1)*(2*n+1)*factd(m)*factd(n)/(m*(m+1)*n*(n+1)*factd(m-1)*factd(n-1));