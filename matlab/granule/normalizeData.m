function dnorm = normalizeData(d)

dnorm = d - min(d);
dnorm = dnorm./max(dnorm);