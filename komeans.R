fun.okc.2= function (data = data, nclust = nclust, lnorm = lnorm, tolerance = tolerance) 
{
    M = nrow(data)
    N = ncol(data)
    K = nclust
    niterations = 50
#    datanorm = apply(data, 2, fun.normalize)
    datanorm = scale(data)
    S = matrix(sample(c(0, 1), M * K, replace = TRUE), M, K)
    S = cbind(S, rep(1, M))
    W = matrix(runif(N * K), K, N)
    W = rbind(W, rep(0, N))
    sse = rep(0, niterations)
    oprevse = exp(70)
    opercentse = 1
    i = 1
    while ((i <= niterations) & (opercentse > tolerance)) {
        for (k in 1:K) {
            sminusk = S[, -k]
            wminusk = W[-k, ]
            s = as.matrix(S[, k])
            w = t(as.matrix(W[k, ]))
            dstar = datanorm - sminusk %*% wminusk
            prevse = exp(70)
            percentse = 1
            l = 1
            while ((l <= niterations) & (percentse > tolerance)) {
                for (m in 1:N) {
                  if (lnorm == 2) {
                    w[1, m] = mean(dstar[s == 1, m], na.rm = TRUE)
                  }
                  if (lnorm == 1) {
                    w[1, m] = median(dstar[s == 1, m], na.rm = TRUE)
                  }
                }
                for (m in 1:M) {
                  if (lnorm == 2) {
                    ss1 = sum((dstar[m, ] - w[1, ])^2, na.rm = TRUE)
                    ss0 = sum((dstar[m, ])^2, na.rm = TRUE)
                  }
                  if (lnorm == 1) {
                    ss1 = sum(abs(dstar[m, ] - w[1, ]), na.rm = TRUE)
                    ss0 = sum(abs(dstar[m, ]), na.rm = TRUE)
                  }
                  if (ss1 <= ss0) {
                    s[m, 1] = 1
                  }
                  if (ss1 > ss0) {
                    s[m, 1] = 0
                  }
                }
                if (sum(s) == 0) {
                  s[sample(1:length(s), 2)] = 1
                }
                if (lnorm == 2) {
                  se = sum((dstar - s %*% w)^2, na.rm = TRUE)
                }
                if (lnorm == 1) {
                  se = sum(abs(dstar - s %*% w), na.rm = TRUE)
                }
                percentse = 1 - se/prevse
                prevse = se
                l = l + 1
            }
            S[, k] = as.vector(s)
            W[k, ] = as.vector(w)
        }
        if (lnorm == 2) 
            sse[i] = sum((datanorm - S %*% W)^2, na.rm = TRUE)/sum((datanorm - 
                mean(datanorm, na.rm = TRUE))^2, na.rm = TRUE)
        if (lnorm == 1) 
            sse[i] = sum(abs(datanorm - S %*% W), na.rm = TRUE)/sum(abs(datanorm - 
                median(datanorm, na.rm = TRUE)), na.rm = TRUE)
        if (lnorm == 2) {
            ose = sum((datanorm - S %*% W)^2, na.rm = TRUE)
        }
        if (lnorm == 1) {
            ose = sum(abs(datanorm - S %*% W), na.rm = TRUE)
        }
        opercentse = (oprevse - ose)/oprevse
        oprevse = ose
        i = i + 1
    }
    if (lnorm == 2) 
        vaf = cor(as.vector(datanorm), as.vector(S %*% W), use = "complete.obs")^2
    if (lnorm == 1) 
        vaf = 1 - sse[i - 1]
     rrr = list(Data = data, Normalized.Data = datanorm, Tolerance = tolerance, 
        Groups = S[, 1:K], Centroids = round(W[1:K, ], 2), SSE.Percent = sse[1:i - 
            1], VAF = vaf)


    return(rrr)
}

komeans=function (data = data, nclust = nclust, lnorm = lnorm, nloops = nloops, tolerance = tolerance, seed = seed) 
{
    prevsse = 100
    set.seed(seed)
    for (i in 1:nloops) {
        z = fun.okc.2(data = data, nclust = nclust, lnorm = lnorm, 
            tolerance = tolerance)
        if (z$SSE.Percent[length(z$SSE.Percent[z$SSE.Percent >  0])] < prevsse) {
            prevsse = z$SSE.Percent[length(z$SSE.Percent[z$SSE.Percent >  0])]
            ind = i
            z.old = z
        }
    }
    return(list(data = z.old$Data, Normalized.Data = z.old$Normalized.Data, 
        Group = z.old$Group %*% as.matrix(2^(0:(nclust-1)) ), Centroids = z.old$Centroids, Tolerance = z.old$Tolerance, 
        SSE.Pecent = z.old$SSE.Percent, VAF = z.old$VAF, iteration = ind, 
        seed = seed))
}
