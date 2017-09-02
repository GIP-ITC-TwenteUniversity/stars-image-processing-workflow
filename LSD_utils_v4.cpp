#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix neighbourhood_systemC(int, int);
NumericVector Grow_region_seed(NumericVector, NumericVector, int, int, int, double, double);

// [[Rcpp::export]]
NumericMatrix convolve_2d(NumericMatrix sample, NumericMatrix kernel){
   
   int x_s = sample.nrow(), x_k = kernel.nrow();
   int y_s = sample.ncol(), y_k = kernel.ncol();

   NumericMatrix output(x_s + x_k - 1, y_s + y_k - 1);
   
   for (int row = 0; row < x_s; row++) {
     for (int col = 0; col < y_s; col++) {
       for (int i = 0; i < x_k; i++) {
         for (int j = 0; j < y_k; j++) {
           output(row + i, col + j) += sample(row, col) * kernel(i, j);
         }
       }
     }
   }
   
   return output;
}

// [[Rcpp::export]]
NumericVector convolve_2d_v3(NumericVector A, int n_col, int n_row, NumericMatrix K, int bad_val){
   
	int k_row = K.nrow();
	int k_col = K.ncol();
	
	int hw_row = (k_row-1)/2;
	int hw_col = (k_col-1)/2;
	int nr,nc;

	NumericVector output(n_row*n_col);

	for(int row = hw_row; row < n_row-hw_row; row++){
		for(int col = hw_col; col < n_col-hw_col; col++){
			if(A[row*n_col + col]==bad_val){
				output[row*n_col + col] = bad_val;
			}else{
				for(int drow = -hw_row; drow < hw_row; drow++){
					for(int dcol = -hw_col; dcol < hw_col; dcol++){
						nr = row+drow;
						nc = col+dcol;
						if(A[nr*n_col + nc]==bad_val){
							output[row*n_col + col] = bad_val;
							// exit the loop
							drow = hw_row;
							dcol = hw_col;
							continue;
						}
						//if(nr>=0 && nr<n_row && nc>=0 && nc <n_col){
							output[row*n_col + col] += A[nr*n_col + nc] * K(drow+hw_row, dcol+hw_col);
						//}
					}
				}
			}
		}
   }
   
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
NumericVector Grow_region_seed(NumericVector f_magn, NumericVector f_angle, int M, int N, int seed, double tau, double magn_min){
   
	int npix = M*N;
	int i, j, k, pn, pn2, grow;
	int count=0;
	double alpha, theta, sx, sy;
	NumericVector marked(npix);
	NumericVector seg(npix);
	NumericMatrix narr(npix,9);
   
	narr = neighbourhood_systemC(M,N);
	
	pn = seed-1; // in C++ array counter starts with zero

	seg[count] = pn;
	marked[pn] = 1;

	count++;

	grow = 1;
	
	while(grow==1){
		for(i=0; i<count; i++){
			grow=0;
			pn = seg[i];

			// mean angle of the current segment
			sx = 0.0;
			sy = 0.0;
			for(k=0;k<count;k++){
				sx+=cos(f_angle[seg[k]]);
				sy+=sin(f_angle[seg[k]]);
			}
			theta = atan2(sy,sx);
			
			for(j=1; j<narr(pn,0)+1;j++){
				pn2 = narr(pn,j);
				
				
				alpha = abs(theta - f_angle[pn2]);
			
				alpha += (alpha>PI) ? -2.0*PI : (alpha<-PI) ? (2.0*PI) : 0.0;	
			
				//if(abs(alpha)<=tau && marked[pn2]==0){
				if(abs(alpha)<=tau && f_magn[pn2]>0.0){
					marked[pn2] = 1;
					f_magn[pn2] = 0.0;
					seg[count]=pn2;
					count++;
					grow=1;
				}
			}
		}
	}

	NumericVector out(count+1);

	out[0] = count;
	for(i=0; i<count; i++){
		out[i+1] = seg[i];//+1; // In R array counter starts with 1
	}
	
	return out;
}

// [[Rcpp::export]]
NumericMatrix segm_map(NumericVector f_magn, NumericVector f_angle, int M, int N, double tau, double magn_min, double min_area){
   
	int npix = M*N;
	int i,j,k, nj, pn, pn2, grow;
	int n_left = 0;
	int count = 0;
	double alpha, theta, sx, sy;
	double magn_max = 0.0;
	NumericVector magn(npix);

	NumericMatrix output(npix,2);

	int seg_length, seg_count = 0;
	NumericVector seg(npix), seg_id(npix);
	NumericVector marked(npix);

	NumericMatrix narr(npix,9);
	
	narr = neighbourhood_systemC(M,N);
	
	// mark the low gradient magnitude pixels as inappropriate for region growing
	for(i=0; i<npix; i++){
		magn[i] = f_magn[i];
		if(magn[i]<=magn_min){
			magn[i] = 0;
		}else{
			n_left++;
			if(magn[i]>magn_max) magn_max = magn[i];
		}
	}
	
	if(magn_max<magn_min){
		return 0;
	}
	
	while(n_left>0){
		//seg_count++;
		count = 0;
		
		magn_max = 0.0;
		
		for(i=0; i<npix; i++){
			if(magn[i]>magn_max){
				magn_max = magn[i];
				pn=i;
			}
		}

		seg[count] = pn;
		marked[pn] = 1;
		magn[pn] = 0.0;

		count++;

		grow = 1;

		while(grow==1){
			// mean angle of the current segment
			sx = 0.0;
			sy = 0.0;
			for(k=0;k<count;k++){
				sx+=cos(f_angle[seg[k]]);
				sy+=sin(f_angle[seg[k]]);
			}

			theta = atan2(sy,sx);

			for(i=0; i<count; i++){
				grow=0;
				pn = seg[i];

				for(j=1; j<narr(pn,0)+1;j++){
					pn2 = narr(pn,j);
					if(magn[pn2]>0.0){
						alpha = abs(theta - f_angle[pn2]);
						alpha += (alpha>PI) ? -2.0*PI : (alpha<-PI) ? (2.0*PI) : 0.0;	

						//if(abs(alpha)<=tau && marked[pn2]==0){
						if(abs(alpha)<=tau){
							marked[pn2] = 1;
							magn[pn2] = 0.0;
							seg[count]=pn2;
							count++;
							grow=1;
						}
					}
				}
				
			}
			
		}

		if(count>min_area){
			seg_count++;				
			for(k = 0; k < count; k++){
				nj = seg[k];
				output(nj,0) = seg_count;  // segment id
				output(nj,1) = count;	  // segment size
			}
		}

		// update n_left
		n_left = 0;
		for(i=0; i<npix; i++) if(magn[i]>=magn_min) n_left++;
	}

	return output;
}

// [[Rcpp::export]]
List gaussian_kernel(double sigma){

	List out;
	
	if(sigma==0){
		return out;
	}

	int	kern_size = ceil(6.0*sigma);

	if(kern_size<3){
		kern_size = 3;
	}

	if(kern_size%2 ==0) kern_size++;

	NumericMatrix W(kern_size,kern_size);
	NumericMatrix Wx(kern_size,kern_size);
	NumericMatrix Wy(kern_size,kern_size);

	int nx = kern_size;
	int ny = kern_size;

	int mid_row = (ny-1)/2;
	int mid_col = (nx-1)/2;

	double gain = 0.0;
	double dx, dy;

	for(int i_row=0; i_row<kern_size; i_row++){
		dy = i_row - mid_row;
		for(int i_col=0; i_col<kern_size; i_col++){
			dx = i_col - mid_col;
			W(i_row,i_col) = exp(-0.5*((dx*dx+dy*dy)/(sigma*sigma)));
			Wx(i_row,i_col) = -W(i_row,i_col) * (dx/(sigma*sigma));
			Wy(i_row,i_col) = -W(i_row,i_col) * (dy/(sigma*sigma));
			gain += W(i_row,i_col);
		}
	}

	for(int i_row=0; i_row<kern_size; i_row++){
		for(int i_col=0; i_col<kern_size; i_col++){
			W(i_row,i_col) 	 /= gain;
			Wx(i_row,i_col) /= gain;
			Wy(i_row,i_col) /= gain;
		}
	}

	// NumericVector Gxx = convolve_2d_v2(A, n_col, n_row, Wxx);
	// NumericVector Gxy = convolve_2d_v2(A, n_col, n_row, Wxy);
	// NumericVector Gyy = convolve_2d_v2(A, n_col, n_row, Wyy);

	// NumericVector out(n_row*n_col);

	// for(int pn=0; pn<n_row*n_col; pn++){
		// out(pn) = (t*t)*(Gxx[pn]*Gyy[pn] - Gxy[pn]*Gxy[pn]);
	// }
	
	out["W"] = W;
	out["Wx"] = Wx;
	out["Wy"] = Wy;

	return out;
}
