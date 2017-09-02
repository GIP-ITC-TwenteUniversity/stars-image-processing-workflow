#include <Rcpp.h>
using namespace Rcpp;

double p_plus_c(NumericMatrix, int);
double p_minus_c(NumericMatrix, int);


// [[Rcpp::export]]
NumericMatrix glcmC(NumericMatrix A, int ql, int dx, int dy, int bad_val){
	int n_row = A.nrow(), n_col=A.ncol();
	
	NumericMatrix GM(ql,ql);
	int count=0;
	int a, b;
		
	for(int i_row=0; i_row < n_row; i_row++){
		for(int i_col=0; i_col < n_col; i_col++){

			a = A(i_row,i_col);

			if(a==bad_val) continue;
			
			if(i_col+dx<n_col && i_col+dx>=0 && i_row+dy<n_row && i_row+dy>=0){

				b = A(i_row+dy,i_col+dx);

				if(b==bad_val) continue;

				GM(a,b)++;
				GM(b,a)++;
				count ++;
			}
		}
	}
	
	/*
	for(a=0; i<ql; i++){
		for(b=0; j<ql; j++){
			GM(a,b)/=count;
		}
	}
	*/

	return GM;
}

// [[Rcpp::export]]
IntegerVector get_pan_pixels(IntegerVector pixels, int M, int N){

	int n = pixels.size();
	int S = 4;	// scale ratio between ms and pan
	int i, pn, pn_pan;
	int i_row, i_col;
	int j_row, j_col;
	int count = 0;

	IntegerVector out_pixels(n*S*S);

	for(i=0; i<n; i++){
		pn = pixels[i];	// ms pixel id
		i_row = floor(pn/M);
		i_col = pn - i_row * M;
		
		j_row = S * i_row;	// top left pan pixel
		j_col = S * i_col;	// top left pan pixel
		
		for(j_row=S * i_row; j_row < S * (i_row+1); j_row++){
			for(j_col=S * i_col; j_col < S * (i_col+1); j_col++){
				pn_pan = j_col + j_row*M*S;
				out_pixels[count] = pn_pan;
				count++;
			}
		}

	}

	return out_pixels;
}

// p_(x+y)(k) in Haralick's paper

// [[Rcpp::export]]
double p_plus_c(NumericMatrix A, int k){
	int ql = A.nrow();

	double val = 0.0;

	for(int i=0; i<ql; i++){
		// this loop is not necessary. test
		for(int j=0; j<ql; j++){
			if(i+j+2==k) val += A(i,j);
		}
	}

	return val;
}

// [[Rcpp::export]]
double p_plus_c_v2(NumericMatrix A, int k){
	int ql = A.nrow();
	double val = 0.0;

	if(k<ql){
		for(int i=0; i<=k; i++){
			val += A(i,k-i);
		}
	} else{
		for(int i=k-ql+1; i<ql; i++){
			val += A(i,k-i);
		}
	}

	return val;
}

// p_(x-y)(k) in Haralick's paper

// [[Rcpp::export]]
double p_minus_c(NumericMatrix A, int k){
	int ql = A.nrow();

	double val = 0.0;

	for(int i=0; i<ql; i++){
		// this loop is not necessary. test
		for(int j=0; j<ql; j++){
			if(abs(i-j)==k) val += A(i,j);
		}
	}

	return val;

}

// [[Rcpp::export]]
double p_minus_c_v2(NumericMatrix A, int k){
	int ql = A.nrow();

	double val = 0.0;
	
	int j;

	for(int i=0; i<ql; i++){
		j = i-k;
		if(j>=0) val += A(i,j);
		
		if(k>0){
			j = i+k;
			if(j<ql) val += A(i,j);
		}
	}

	return val;
}

// [[Rcpp::export]]
NumericVector p_plus_arr_c(NumericMatrix A){
	int ql = A.nrow();

	int n = 2*ql - 1;
	NumericVector out(n);
	
	for(int i=0; i<n; i++){
		out[i] = p_plus_c(A, i+2);
	}

	return out;
}

// [[Rcpp::export]]
NumericVector p_plus_arr_c_v2(NumericMatrix A){
	int ql = A.nrow();

	int n = 2*ql-1;
	NumericVector out(n);
	
	for(int k=0; k<n; k++){
		out[k] = p_plus_c_v2(A,k);
	}

	return out;
}

// [[Rcpp::export]]
NumericVector p_minus_arr_c(NumericMatrix A){
	int ql = A.nrow();

	NumericVector out(ql);
	
	for(int i=0; i<ql; i++){
		out[i] = p_minus_c(A, i);
	}

	return out;
}

// [[Rcpp::export]]
NumericVector p_minus_arr_c_v2(NumericMatrix A){
	int ql = A.nrow();

	NumericVector out(ql);
	
	for(int i=0; i<ql; i++){
		out[i] = p_minus_c_v2(A, i);
	}

	return out;
}

// [[Rcpp::export]]
List f12(NumericMatrix A, NumericVector px, NumericVector py, double HXY){

	// f12 in Haralick et al. (1973)
	// Information measures of correlation (1)

	int ql = A.nrow();
	double HX=0.0, HY=0.0, tmp1=0.0, tmp2=0.0;
	double HXY1, HXY2;
	double p, pxy;
	
	List out;

	for(int i=0; i<ql; i++){
		p = px[i];
		if(p>0){
			HX -= p*log2(p);
		}

		p = py[i];
		if(p>0){
			HY -= p*log2(p);
		}
	}

	for(int i=0; i<ql; i++){
		for(int j=0; j<ql; j++){
			pxy = px[i]*py[j];
			if(pxy>0.0){
				tmp1 -= A(j,i)*log2(pxy);
				tmp2 -= pxy*log2(pxy);
			}
			
		}
	}

	HXY1 = tmp1;
	HXY2 = tmp2;

	if(HY>HX) HX = HY;
	out["HXY1"] = HXY1;
	out["HXY2"] = HXY2;
	out["HX"] = HX;
	out["HY"] = HY;
	out["f12"] = (HXY - HXY1)/HX;
	out["f13"] = sqrt(1.0-exp(-2.0*(HXY2-HXY)));

	return out;
}


// [[Rcpp::export]]
NumericMatrix Qfun(NumericMatrix A, NumericVector px, NumericVector py){

	// f12 in Haralick et al. (1973)
	// Information measures of correlation (1)

	int ql = A.nrow();
	double pxy;

	double tmp = 0.0;
	NumericMatrix Q(ql,ql);

	
	for(int i=0; i<ql; i++){
		for(int j=0; j<ql; j++){
			tmp = 0.0;
			for(int k=0; k<ql; k++){
				pxy = px[i]*py[k];
				if(pxy>0.0) tmp += A(i,k)*A(j,k)/pxy;
			}
			Q(i,j) += tmp;
		}
	}

	return Q;
}

/*
haralick_textures <- function(A){

	ql <- dim(A)[1]
	if(dim(A)[2]!=ql) stop("Not square GCLM")

	rownames(A) <- colnames(A) <- 0:(ql-1)
	ind <- which(rowSums(A)>0)
	A <- A[ind,ind]

	eps <- 1.0e-32

	Ntex <- 14
	out <- array(NA,Ntex)
	names(out) <- tex_names

	px <- colSums(A)
	py <- rowSums(A)

	# f1 in Haralick et al. (1973)
	# Angular Second Moment
	# Note: some sources apply the square root to get energy from ASM
	out[1] <- sum(A^2)

	# f2 in Haralick et al. (1973)
	# formula is from elsewhere, simplified computationally
	out[2] <- sum(A*((row(A)-col(A))^2))

	# f3 in Haralick et al. (1973)
	mx <- sum(px*(0:(ql-1)))
	my <- sum(py*(0:(ql-1)))

	sx <- sqrt(sum(px*(((0:(ql-1))-mx)^2)))
	sy <- sqrt(sum(py*(((0:(ql-1))-my)^2)))

	if(sx*sy>eps){
		val <- sum((row(A)-1 - my)*(col(A)-1 - mx)*A)
		out[3] <- val / (sx*sy)
	}

	
	# f4 in Haralick et al. (1973)
	# Sum of Squares: variance
	# Note: m is not defined in the source. Use mx
	# Note: the summation is from 0 to ql-1, not from 1 to ql
	out[4] <- sum(A*(((col(A)-1)-mx)^2))

	# f5 in Haralick et al. (1973)
	# Inverse Difference Moment
	out[5] <- sum(A/(1+((row(A)-col(A))^2)))

	# f6 in Haralick et al. (1973)
	# Sum Average
	val <- 0
	for(i in 2:(2*ql)) val <- val + i * p_plus(A, i)
	out[6] <- val

	iarr <- 2:(2*ql)
	p_plus_arr <- array(0,length(iarr))
	for(a in 1:length(iarr)){
		p_plus_arr[a] <- p_plus(A, iarr[a])
	}

	karr <- 0:(ql-1)
	p_minus_arr <- array(0,length(karr))
	for(a in 1:length(karr)){
		p_minus_arr[a] <- p_minus(A, karr[a])
	}


	# f7 in Haralick et al. (1973)
	# Sum Variance
	mu <- out[6]
#	val <- 0
#	for(i in 2:(2*ql)) val <- val + ((i-mu)^2) * p_plus(A, i)
#	out[7] <- val
	out[7] <- sum(((iarr-mu)^2) * p_plus_arr)

	# f8 in Haralick et al. (1973)
	# Sum Entropy
	# Remark: different log bases are found in literature
	# Follow base 2 here
	val <- p_plus_arr * log2(p_plus_arr)
	out[8] <- -sum(val[p_plus_arr>0])

	# f9 in Haralick et al. (1973)
	# Remark: different log bases are found in literature
	# Follow base 2 here
	val <- as.vector(A)
	val <- val[val>0]
	val <- - val*log2(val)
	out[9] <- sum(val)

	# f10 in Haralick et al. (1973)
	# Difference variance
	out[10] <- var(p_minus_arr)

	# f11 in Haralick et al. (1973)
	# Difference entropy
	val <- p_minus_arr
	val <- val[val>0]
	val <- val*log2(val)
	out[11] <- -sum(val)

	# f12 in Haralick et al. (1973)
	# Information measures of correlation (1)
	HXY  <- out[9]

	px1 <- px[px>0]
	py1 <- py[py>0]

	HX <- -sum(px1*log2(px1))
	HY <- -sum(py1*log2(py1))

	tmp1 <- 0
	tmp2 <- 0
	for(i in 1:ql){
		for(j in 1:ql){
			pxy <- px[i]*py[j]
			if(pxy>0){
				tmp1 <- tmp1 - A[j,i]*log2(pxy)
				tmp2 <- tmp2 - pxy*log2(pxy)
			}
		}
	}
	
	HXY1 <- tmp1
	HXY2 <- tmp2

	out[12] <- (HXY - HXY1)/max(c(HX,HY))

	# f13 in Haralick et al. (1973)
	out[13] <- sqrt(1-exp(-2.0*(HXY2-HXY)))

	# f14 in Haralick et al. (1973)
	# Maximal correlation coefficient
	Q <- array(0,c(ql,ql))

	for(k in 1:ql){
		tmp <- array(0,c(ql,ql))
		for(i in 1:ql){
			for(j in 1:ql){
				pxy <- px[i]*py[k]
				if(pxy>0) tmp[i,j] <- tmp[i,j] + (A[i,k]*A[j,k])/pxy
			}
		}
		Q[,] <- Q[,] + tmp[,]
	}
	val <- eigen(Q,only.values=TRUE)$values
	out[14] <- val[2]

	return(out)
}


*/
