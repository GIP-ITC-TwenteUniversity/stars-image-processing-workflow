#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix neighbourhood_systemC(int, int);
NumericMatrix neighbourhood_system_MRF(int, int);

// [[Rcpp::export]]
NumericMatrix neighbourhood_system_MRF(int M, int N){
   
	int npix = M*N;
	int i, j, pn;
	
	NumericMatrix out(npix,17);
   
	for(i = 1; i < M-1; i++){
		for(j = 1; j < N-1; j++){
			pn = j*M + i;
			out(pn,0) = 8;
			out(pn,1) = pn - M - 1;
			out(pn,2) = 3;
			out(pn,3) = pn - M;
			out(pn,4) = 1;
			out(pn,5) = pn - M + 1;
			out(pn,6) = 2;
			out(pn,7) = pn - 1;
			out(pn,8) = 0;
			out(pn,9) = pn + 1;
			out(pn,10) = 0;
			out(pn,11) = pn + M - 1;
			out(pn,12) = 2;
			out(pn,13) = pn + M;
			out(pn,14) = 1;
			out(pn,15) = pn + M + 1;
			out(pn,16) = 3;
		}
		
		j = 0;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - 1;
		out(pn,2) = 0;
		out(pn,3) = pn + 1;
		out(pn,4) = 0;
		out(pn,5) = pn + M - 1;
		out(pn,6) = 2;
		out(pn,7) = pn + M;
		out(pn,8) = 1;
		out(pn,9) = pn + M + 1;
		out(pn,10) = 3;
		
		j = N-1;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = 3;
		out(pn,3) = pn - M;
		out(pn,4) = 1;
		out(pn,5) = pn - M + 1;
		out(pn,6) = 2;
		out(pn,7) = pn - 1;
		out(pn,8) = 0;
		out(pn,9) = pn + 1;
		out(pn,10) = 0;
	}
	
	i=0;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M;
		out(pn,2) = 1;
		out(pn,3) = pn - M + 1;
		out(pn,4) = 2;
		out(pn,5) = pn + 1;
		out(pn,6) = 0;
		out(pn,7) = pn + M;
		out(pn,8) = 1;
		out(pn,9) = pn + M + 1;
		out(pn,10) = 3;
	}

	j = 0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn + 1;
	out(pn,2) = 0;
	out(pn,3) = pn + M;
	out(pn,4) = 1;
	out(pn,5) = pn + M + 1;
	out(pn,6) = 3;
	
	
	j = N - 1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M;
	out(pn,2) = 1;
	out(pn,3) = pn - M + 1;
	out(pn,4) = 2;
	out(pn,5) = pn + 1;
	out(pn,6) = 0;

	i = M - 1;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = 3;
		out(pn,3) = pn - M;
		out(pn,4) = 1;
		out(pn,5) = pn - 1;
		out(pn,6) = 0;
		out(pn,7) = pn + M - 1;
		out(pn,8) = 2;
		out(pn,9) = pn + M;
		out(pn,10) = 1;
	}
	
	j=0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - 1;
	out(pn,2) = 0;
	out(pn,3) = pn + M - 1;
	out(pn,4) = 2;
	out(pn,5) = pn + M;
	out(pn,6) = 1;

	j=N-1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M - 1;
	out(pn,2) = 3;
	out(pn,3) = pn - M;
	out(pn,4) = 1;
	out(pn,5) = pn - 1;
	out(pn,6) = 0;

	return out;
}

// [[Rcpp::export]]
NumericVector Uprior_v3(NumericVector f1, NumericVector f2, int M, int N, NumericMatrix nsyst, NumericVector beta){
   
	int npix = M*N;
	int i, j, k, pn;
	double sum;
	NumericVector output(npix);
	//NumericMatrix nsyst(npix,17);
	//nsyst = neighbourhood_system_MRF(M,N);
   
	for(i = 0; i < M; i++){
		for(j = 0; j < N; j++){
			pn = j*M + i;
			sum = 0.0;

			for(k=1; k<2*nsyst(pn,0); k+=2){
			
				if(f1[pn]!=f2[nsyst(pn,k)]) sum += beta[nsyst(pn,k+1)];
			}

			output[pn] += sum;
		}
	}

	return output;
}

// [[Rcpp::export]]
NumericVector Uprior_v2(NumericVector f1, NumericVector f2, int M, int N, NumericVector beta){
   
	int npix = M*N;
	int i, j, pn;
	double sum;
	NumericVector output(npix);
   
	for(i = 1; i < M-1; i++){
		for(j = 1; j < N-1; j++){
			pn = j*M + i;
			sum= 0.0;

			if(f1[pn]!=f2[pn-1]) sum += beta[0];
			if(f1[pn]!=f2[pn+1]) sum += beta[0];

			if(f1[pn]!=f2[pn-M]) sum += beta[1];
			if(f1[pn]!=f2[pn+M]) sum += beta[1];

			if(f1[pn]!=f2[pn-M+1]) sum += beta[2];
			if(f1[pn]!=f2[pn+M-1]) sum += beta[2];

			if(f1[pn]!=f2[pn-M-1]) sum += beta[3];
			if(f1[pn]!=f2[pn+M+1]) sum += beta[3];

			output[pn] += sum;
		}

		j = 0;
		pn = j*M + i;
		sum = 0.0;
			
		if(f1[pn]!=f2[pn-1]) sum += beta[0];
		if(f1[pn]!=f2[pn+1]) sum += beta[0];

		if(f1[pn]!=f2[pn+M]) sum += beta[1];

		if(f1[pn]!=f2[pn+M-1]) sum += beta[2];

		if(f1[pn]!=f2[pn+M+1]) sum += beta[3];

		output[pn] += sum;

		j=N-1;
		pn = j*M + i;
		sum= 0.0;

		if(f1[pn]!=f2[pn-1]) sum += beta[0];
		if(f1[pn]!=f2[pn+1]) sum += beta[0];

		if(f1[pn]!=f2[pn-M]) sum += beta[1];

		if(f1[pn]!=f2[pn-M+1]) sum += beta[2];

		if(f1[pn]!=f2[pn-M-1]) sum += beta[3];

		output[pn] += sum;
   }
   
   i=0;
   for(j = 1; j < N-1; j++){
		pn = j*M + i;
		sum= 0.0;

		if(f1[pn]!=f2[pn+1]) sum += beta[0];

		if(f1[pn]!=f2[pn-M]) sum += beta[1];
		if(f1[pn]!=f2[pn+M]) sum += beta[1];

		if(f1[pn]!=f2[pn-M+1]) sum += beta[2];

		if(f1[pn]!=f2[pn+M+1]) sum += beta[3];

		output[pn] += sum;
   }
   
	j=0;
	pn = j*M + i;
	sum= 0.0;

	if(f1[pn]!=f2[pn+1]) sum += beta[0];

	if(f1[pn]!=f2[pn+M]) sum += beta[1];

	if(f1[pn]!=f2[pn+M+1]) sum += beta[3];

	output[pn] += sum;

	j=N-1;
	pn = j*M + i;
	sum= 0.0;

	if(f1[pn]!=f2[pn+1]) sum += beta[0];

	if(f1[pn]!=f2[pn-M]) sum += beta[1];

	if(f1[pn]!=f2[pn-M+1]) sum += beta[2];

	output[pn] += sum;

	i=M-1;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		sum= 0.0;

		if(f1[pn]!=f2[pn-1]) sum += beta[0];

		if(f1[pn]!=f2[pn-M]) sum += beta[1];
		if(f1[pn]!=f2[pn+M]) sum += beta[1];

		if(f1[pn]!=f2[pn+M-1]) sum += beta[2];

		if(f1[pn]!=f2[pn-M-1]) sum += beta[3];

		output[pn] += sum;
	}

	j=0;

	pn = j*M + i;
	sum= 0.0;

	if(f1[pn]!=f2[pn-1]) sum += beta[0];

	if(f1[pn]!=f2[pn+M]) sum += beta[1];

	if(f1[pn]!=f2[pn+M-1]) sum += beta[2];

	output[pn] += sum;
	
	
	j=N-1;
	
	pn = j*M + i;
	sum= 0.0;

	if(f1[pn]!=f2[pn-1]) sum += beta[0];

	if(f1[pn]!=f2[pn-M]) sum += beta[1];

	if(f1[pn]!=f2[pn-M-1]) sum += beta[3];

	output[pn] += sum;

	return output;
}

// [[Rcpp::export]]
NumericMatrix neighbourhood_systemC(int M, int N){
   
	int npix = M*N;
	int i, j, pn;
	
	NumericMatrix out(npix,9);
   
	for(i = 1; i < M-1; i++){
		for(j = 1; j < N-1; j++){
			pn = j*M + i;
			out(pn,0) = 8;
			out(pn,1) = pn - M - 1;
			out(pn,2) = pn - M;
			out(pn,3) = pn - M + 1;
			out(pn,4) = pn - 1;
			out(pn,5) = pn + 1;
			out(pn,6) = pn + M - 1;
			out(pn,7) = pn + M;
			out(pn,8) = pn + M + 1;
		}
		
		j = 0;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - 1;
		out(pn,2) = pn + 1;
		out(pn,3) = pn + M - 1;
		out(pn,4) = pn + M;
		out(pn,5) = pn + M + 1;
		
		j = N-1;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = pn - M;
		out(pn,3) = pn - M + 1;
		out(pn,4) = pn - 1;
		out(pn,5) = pn + 1;
	}
	
	i=0;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M;
		out(pn,2) = pn - M + 1;
		out(pn,3) = pn + 1;
		out(pn,4) = pn + M;
		out(pn,5) = pn + M + 1;
	}

	j = 0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn + 1;
	out(pn,2) = pn + M;
	out(pn,3) = pn + M + 1;
	
	
	j = N - 1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M;
	out(pn,2) = pn - M + 1;
	out(pn,3) = pn + 1;

	i = M - 1;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = pn - M;
		out(pn,3) = pn - 1;
		out(pn,4) = pn + M - 1;
		out(pn,5) = pn + M;
	}
	
	j=0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - 1;
	out(pn,2) = pn + M - 1;
	out(pn,3) = pn + M;

	j=N-1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M - 1;
	out(pn,2) = pn - M;
	out(pn,3) = pn - 1;

	return out;
}

// [[Rcpp::export]]
NumericMatrix neighbourhood_system4C(int M, int N){
   
	int npix = M*N;
	int i, j, pn;
	
	NumericMatrix out(npix,5);
   
	for(i = 1; i < M-1; i++){
		for(j = 1; j < N-1; j++){
			pn = j*M + i;
			out(pn,0) = 4;
			out(pn,1) = pn - M;
			out(pn,2) = pn - 1;
			out(pn,3) = pn + 1;
			out(pn,4) = pn + M;
		}
		
		j = 0;
		pn = j*M + i;
		out(pn,0) = 3;
		out(pn,1) = pn - 1;
		out(pn,2) = pn + 1;
		out(pn,3) = pn + M;
		
		j = N-1;
		pn = j*M + i;
		out(pn,0) = 3;
		out(pn,1) = pn - M;
		out(pn,2) = pn - 1;
		out(pn,3) = pn + 1;
	}
	
	i=0;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 3;
		out(pn,1) = pn - M;
		out(pn,2) = pn + 1;
		out(pn,3) = pn + M;
	}

	j = 0;
	pn = j*M + i;
	out(pn,0) = 2;
	out(pn,1) = pn + 1;
	out(pn,2) = pn + M;
	
	
	j = N - 1;
	pn = j*M + i;
	out(pn,0) = 2;
	out(pn,1) = pn - M;
	out(pn,2) = pn + 1;

	i = M - 1;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 3;
		out(pn,1) = pn - M;
		out(pn,2) = pn - 1;
		out(pn,3) = pn + M;
	}
	
	j=0;
	pn = j*M + i;
	out(pn,0) = 2;
	out(pn,1) = pn - 1;
	out(pn,2) = pn + M;

	j=N-1;
	pn = j*M + i;
	out(pn,0) = 2;
	out(pn,1) = pn - M;
	out(pn,2) = pn - 1;

	return out;
}

// [[Rcpp::export]]
NumericVector Grow_region_seedC(NumericVector f, int M, int N, int seed){
   
	int npix = M*N;
	int i, j, pn, pn2, cl, grow;
	int count;
	NumericVector marked(npix);
	NumericVector seg(npix);
	NumericMatrix narr(npix,9);
   
	narr = neighbourhood_systemC(M,N);
   
	pn = seed-1; // in C++ array counter starts with zero
	cl = f[pn];

	seg[count] = pn;
	marked[pn] = 1;

	count = 1;

	grow = 1;
	
	while(grow==1){
		for(i=0; i<count; i++){
			grow=0;
			pn = seg[i];
			for(j=1; j<narr(pn,0)+1;j++){
				pn2 = narr(pn,j);
				if(f[pn2]==cl && marked[pn2]==0){
					marked[pn2] = 1;
					seg[count]=pn2;
					count++;
					grow=1;
				}
			}
		}
	}

	NumericVector out(count);

	for(i=0; i<count; i++){
		out[i] = seg[i]+1; // In R array counter starts with 1
	}
	
	return out;
}

// [[Rcpp::export]]
int SegmentC(NumericVector f, int M, int N){
   
	int npix = M*N;
	int i, pn, ns;
	int pix_count=0, seg_count=0;
	NumericVector marked(npix);
	NumericVector seg(npix),out(5*npix);
	NumericMatrix narr(npix,9);
   
	narr = neighbourhood_systemC(M,N);

	//for(pn=0; pn<npix; pn++){
		pn=0;
		if(marked[pn]==0){
			seg = Grow_region_seedC(f, M, N, pn);
			seg_count++;
			ns = seg[0];
			out[pix_count] = ns;
			pix_count++;
			
			for(i=0; i<seg.size(); i++){
				out[pix_count] = seg[i];
				marked[seg[i]] = 1;
				pix_count++;
			}
		}
	//}
	
	/*NumericVector out2(pix_count);
	
	for(i=0; i<pix_count; i++){
		out2[i] = out[i]+1; // In R array counter starts with 1
	}
	*/
	
	return pix_count;
}


/*

	#u_obj1 <- array(0,npix)
	#u_obj2 <- u_obj1
	
	#upd_path <- array(0,npix)
	
	#for(pn in 1:npix){
		#if(upd_path[pn]==1) next
		#seg <- Grow_region_seed(f,M,N,pn)

		#cl <- f[pn]

		#if((cl==3)|(cl==4)|(cl==5)){
		#	tmp <- area_pdf(length(seg),meanA[cl],sdA[cl])
		#}

		#u_obj1[seg] <- tmp

		#upd_path[seg] <- 1
	#}

	#upd_path <- array(0,npix)
	
	#for(pn in 1:npix){
		#if(upd_path[pn]==1) next
		#f1 <- f
		#f1[pn] <- f_new[pn]
		
		#seg <- Grow_region_seed(f1,M,N,pn)

		#cl <- f[pn]

		#tmp <- area_pdf(length(seg),meanA[cl],sdA[cl])

		#u_obj2[seg] <- tmp

		#upd_path[seg] <- 1
	#}
	
	
	#MRF$obj <- u_obj1
	#image(MRF,attr="obj")
	
	#u1 <- u1 + u_obj1
	#u2 <- u2 + u_obj2


*/

