library(Matrix)

GeCo = function(Y, VS, fixed_cov = NA) {

	N <- dim(Y)[1]; G <- dim(Y)[2]; C <- dim(Y)[3]

  	y = vectorized(Y)
  	if (is.na(fixed_cov)[1]) {
    	X = kronecker(matrix(1,ncol = 1, nrow = N),diag(G*C))
    } else {
      	X = cbind(kronecker(kronecker(as.matrix(fixed_cov),diag(G)), matrix(1,nrow = C, ncol = 1)), 
        	kronecker(matrix(1,nrow=N,ncol=1), diag(G*C)))
    }
  	qrX = qr(X); y = qr.resid(qrX, y)
  	Y = aperm(array(y, dim = c(C,G,N)), perm = c(3, 2, 1))
  	 
	S <- apply(Y, c(1,2), sum)

	Valpha  <- matrix(NA_real_, G,    G)
	SEalpha <- matrix(NA_real_, G,    G)
	Vgc     <- matrix(NA_real_, G*C, G*C)
	SEgc    <- matrix(NA_real_, G*C, G*C)

	pairs <- which(upper.tri(matrix(TRUE, G, G), diag=TRUE), arr.ind=TRUE)

	for(k in seq_len(nrow(pairs))) {
  		j  <- pairs[k,1]
  		jj <- pairs[k,2]

  		# per‐i vector a
  		a   <- (S[,j] * S[,jj] - rowSums(Y[,j,] * Y[,jj,])) / (C*(C-1))
  		mean_a <- mean(a)
  		se_a   <- sqrt(var(a)/N)

  		Valpha[j,jj]   <- Valpha[jj,j]   <- mean_a
  		SEalpha[j,jj]  <- SEalpha[jj,j]  <- se_a

  		E  <- Y[,j,] * Y[,jj,] - VS[,,j,jj]    # N×C
  		Gm <- sweep(E, 1, a, `-`)               # subtract a[i] from each row

  		m_g <- colMeans(Gm)                    # length-C
  		s_g <- sqrt(apply(Gm, 2, var)/N)       # length-C

  		rows_j  <- (0:(C-1))*G + j
  		cols_jj <- (0:(C-1))*G + jj

  		Vgc   [cbind(rows_j, cols_jj)] <- m_g
  		Vgc   [cbind(cols_jj, rows_j)] <- m_g
  		SEgc  [cbind(rows_j, cols_jj)] <- s_g
  		SEgc  [cbind(cols_jj, rows_j)] <- s_g
	}

	list('V_alpha' = Valpha, 'V_GC' = Vgc, 'SE_alpha' = SEalpha, 'SE_GC' = SEgc)


}



