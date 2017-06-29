segmentbis <- function (Y_qc, Yhat, optK, K, sampname_qc, ref_qc, chr, lmax, mode){
    finalcall <- matrix(ncol = 9)
    lmax <- lmax - 1
    for (sampno in 1:ncol(Y_qc)) {
        message("Segmenting sample ", sampno, ": ", sampname_qc[sampno], 
            ".")
        y <- Y_qc[, sampno]
        yhat <- Yhat[[which(K == optK)]][, sampno]
        num <- length(y)
        y <- c(y, rep(0, lmax))
        yhat <- c(yhat, rep(0, lmax))
        i <- rep(1:num, rep(lmax, num))
        j <- rep(1:lmax, num) + i-1 # removed 1 to include segments of size 1 exon
        yact <- rep(0, length(i))
        lambda <- rep(0, length(i))
        for (k in 1:num) {
            yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax-1)]) # keep first element in order to get the # reads mapped to the sole exon k
            lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax-1)]) # keep first element in order to get the # reads mapped to the sole exon k
        }
        i <- i[j <= num]
        yact <- yact[j <= num]
        lambda <- lambda[j <= num]
        j <- j[j <= num]
        yact[lambda < 20] <- 20
        lambda[lambda < 20] <- 20
        if (mode == "integer") {
            chat <- round(2 * (yact/lambda))
        }
        else if (mode == "fraction") {
            chat <- 2 * (yact/lambda)
        }
        lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001)*yact
        #chat[chat > 5] <- 5 #comment to obtain the real duplication number instead of a truncated one
        if (length(lratio) > 0) { # include lratio <=0
            if (length(lratio) >= 2) { # include lratio <=0
                finalmat <- (cbind(i, j, yact, lambda, chat, lratio))#[lratio > 0, ] #keep log ratios <0 (they probably include likelihood estimation errors but are usually small)
                finalmat <- finalmat[order(-finalmat[, 6]), ]
                s <- 1
                while (s <= (nrow(finalmat))) {
                  rowstart <- finalmat[s, 1]
                  rowend <- finalmat[s, 2]
                  rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= rowstart)
                  rowsel[s] <- FALSE
                  finalmat <- finalmat[!rowsel, ]
                  if (is.vector(finalmat)) {
                    finalmat <- t(as.matrix(finalmat))
                  }
                  s <- s + 1
                }
            }
            if (length(lratio) == 1) {
                finalmat <- (cbind(i, j, yact, lambda, chat, lratio))#[lratio > 0, ] #keep log ratios <0
                finalmat <- t(as.matrix(finalmat))
            }
            finalmat <- round(finalmat, digits = 3)
            loglikeij <- cumsum(finalmat[, 6])
            mBIC <- rep(NA, length(loglikeij))
            for (s in 1:nrow(finalmat)) {
                tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
                P <- length(tau) - 2
                mbic <- loglikeij[s]
                mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
                mbic <- mbic + (0.5 - P) * log(num)
                mBIC[s] <- mbic
            }
            mBIC <- round(mBIC, digits = 3)
            if (mBIC[1] > 0) {
                finalmat <- cbind(rep(sampname_qc[sampno], nrow(finalmat)), rep(chr, nrow(finalmat)), finalmat)
                finalmat <- cbind(finalmat, mBIC) #[1:which.max(mBIC), ]) #also print CNVs with low BIC (e.g., non-significant)
                finalmat[(which.max(mBIC)+1):length(mBIC),7] <- 2 #NEWLINE: put all CNVs after the max mBIC at 2
                finalcall <- rbind(finalcall, finalmat)
            }
        }
    }
    finalcall <- finalcall[-1, ]
    if (is.vector(finalcall)) {
        finalcall <- t(as.matrix(finalcall))
    }
    st <- start(ref_qc)[as.numeric(finalcall[, 3])]
    ed <- end(ref_qc)[as.numeric(finalcall[, 4])]
    cnvtype <- rep(NA, length(st))
    cnvtype[finalcall[, 7] < 2] <- "del"
    cnvtype[finalcall[, 7] > 2] <- "dup"
    if (nrow(finalcall) == 1) {
        finalcall <- t(as.matrix(c(finalcall[, 1:2], cnvtype, st, ed, (ed - st + 1)/1000, finalcall[, 3:9])))
    }
    else {
        finalcall <- cbind(finalcall[, 1:2], cnvtype, st, ed, (ed - st + 1)/1000, finalcall[, 3:9])
    }
    colnames(finalcall) <- c("sample_name", "chr", "cnv", "st_bp", "ed_bp", "length_kb", "st_exon", "ed_exon", "raw_cov", "norm_cov", "copy_no", "lratio", "mBIC")
    rownames(finalcall) <- rep("", nrow(finalcall))
    lratio = as.numeric(finalcall[, "lratio"])
    pvalue = signif(pchisq(q = lratio * 2, df = 1, lower.tail = FALSE), 3)
    finalcall = cbind(finalcall, pvalue)
    finalcall
}
