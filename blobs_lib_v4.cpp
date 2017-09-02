#include <Rcpp.h>
using namespace Rcpp;

NumericVector det_HessianC(NumericVector, int, int, double);

NumericMatrix convolve_2d(NumericMatrix, NumericMatrix);
NumericVector convolve_2d_v2(NumericVector, int, int, NumericMatrix);

NumericVector parabol_interpC(NumericVector x, NumericVector y){
	NumericVector a(3);
	a[2] = (2.0*(y[1]-y[0])*(x[2]+x[0])-(y[2]-y[0])*(x[1]+x[0]))/(2.0*(2.0*y[1]-y[0]-y[2]));
	a[1] = (y[1]-y[0])/((x[1]-x[0])*(x[0]+x[1]-2*a[2]));
	a[0] = y[0] - a[1] * pow((x[0]-a[2]),2);
	return a;
}

// [[Rcpp::export]]
NumericMatrix convolve_2d_Eigen(NumericMatrix sample, NumericMatrix kernel){
   
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
NumericMatrix convolve_2d(NumericMatrix A, NumericMatrix K){
   
	int n_row = A.nrow(), k_row = K.nrow();
	int n_col = A.ncol(), k_col = K.ncol();
	
	int hw_row = (k_row-1)/2 ,hw_col = (k_col-1)/2;
	int nr,nc;

	NumericMatrix output(n_row, n_col);

	for(int row = 0; row < n_row; row++){
		for(int col = 0; col < n_col; col++){
			for(int drow = -hw_row; drow < hw_row; drow++){
				for(int dcol = -hw_col; dcol < hw_col; dcol++){
					nr = row+drow;
					nc = col+dcol;
					if(nr>=0 && nr<n_row && nc>=0 && nc <n_col){
						output(row, col) += A(nr, nc) * K(drow+hw_row, dcol+hw_col);
					}
				}
			}
		}
   }
   
   return output;
}


// [[Rcpp::export]]
NumericVector convolve_2d_v2(NumericVector A, int n_col, int n_row, NumericMatrix K){
   
	int k_row = K.nrow();
	int k_col = K.ncol();
	
	int hw_row = (k_row-1)/2;
	int hw_col = (k_col-1)/2;
	int nr,nc;

	NumericVector output(n_row*n_col);

	for(int row = hw_row; row < n_row-hw_row; row++){
		for(int col = hw_col; col < n_col-hw_col; col++){
			for(int drow = -hw_row; drow < hw_row; drow++){
				for(int dcol = -hw_col; dcol < hw_col; dcol++){
					nr = row+drow;
					nc = col+dcol;
					//if(nr>=0 && nr<n_row && nc>=0 && nc <n_col){
						output[row*n_col + col] += A[nr*n_col + nc] * K(drow+hw_row, dcol+hw_col);
					//}
				}
			}
		}
   }
   
   return output;
}

// [[Rcpp::export]]
NumericMatrix local_max_3d(NumericVector A, int n_col, int n_row, int n_scale) {
  int n = A.size();
  NumericMatrix out(n,3);
  int count=0;

  for(int k = 1; k < n_scale-1; k++) {
  	for(int j = 1; j < n_row-1; j++) {
  		for(int i = 1; i < n_col-1; i++) {
  			if(A[k*n_col*n_row+j*n_col+i]>A[(k-1)*n_col*n_row+j*n_col+i] && A[k*n_col*n_row+j*n_col+i]>A[(k+1)*n_col*n_row+j*n_col+i]){
  				if(A[k*n_col*n_row+j*n_col+i]>A[k*n_col*n_row+(j-1)*n_col+i] && A[k*n_col*n_row+j*n_col+i]>A[k*n_col*n_row+(j+1)*n_col+i]){
  					if(A[k*n_col*n_row+j*n_col+i]>A[k*n_col*n_row+j*n_col+i-1] && A[k*n_col*n_row+j*n_col+i]>A[k*n_col*n_row+j*n_col+i+1]){
  						out(count,0) = i+1;
  						out(count,1) = j+1;
  						out(count,2) = k+1;
  						count++;
  					}
  				}
  			}
  		}	
      }
   }
   
   NumericMatrix out2(count,3);
   
   for(int i=0; i <count; i++){
   		out2(i,0) = out(i,0);
   		out2(i,1) = out(i,1);
   		out2(i,2) = out(i,2);
	}      
  return out2;
}

// [[Rcpp::export]]
NumericMatrix interpolate_blob(NumericVector DetHes, int nx, int ny, int nz, NumericMatrix Pts, NumericVector ps, double xTL, double yTL, NumericVector tarr){
	int npts = Pts.nrow();
	NumericMatrix out(npts,4);
	int i,j,k,pn;
	NumericVector x(3),y(3),a(3);
	double x0,y0,radius,sum;
	int pix_count;

	int i_row, i_col;
	int min_row,max_row,min_col,max_col;
	double weight;
	double sum_w;
	double dist;


  for(int pt=0; pt<npts; pt++){
  	i = Pts(pt,0)-1;
	j = Pts(pt,1)-1;
	k = Pts(pt,2)-1;
	
	pn = j*nx + i;

	// Scale
	x[0] = tarr[k-1];
	x[1] = tarr[k];
	x[2] = tarr[k+1];

	y[0] = DetHes[(k-1)*nx*ny+j*nx+i];
	y[1] = DetHes[k*nx*ny+j*nx+i];
	y[2] = DetHes[(k+1)*nx*ny+j*nx+i];
	
	a = parabol_interpC(x,y);
	
	radius = sqrt(2.0*a[2])*0.5*(ps[0]+ps[1]);
	out(pt,2) = 2.0*radius;
	out(pt,3) = a[0];

	// Position x
	x[0] = (double)-1.0;
	x[1] = (double)0.0;
	x[2] = (double)1.0;

	y[0] = DetHes[k*nx*ny+j*nx+i-1];
	y[1] = DetHes[k*nx*ny+j*nx+i];
	y[2] = DetHes[k*nx*ny+j*nx+i+1];
	
	a = parabol_interpC(x,y);
	
	x0 = a[2]*ps[0] + xTL+i*ps[0];

	out(pt,0) = x0;
	
	// Position y
	x[0] = (double)-1.0;
	x[1] = (double)0.0;
	x[2] = (double)1.0;

	y[0] = DetHes[k*nx*ny+(j-1)*nx+i];
	y[1] = DetHes[k*nx*ny+j*nx+i];
	y[2] = DetHes[k*nx*ny+(j+1)*nx+i];
	
	a = parabol_interpC(x,y);
	
	y0 = -a[2]*ps[1]+(yTL-j*ps[1]);

	out(pt,1) = y0;

	//ndvi
	/*
	sum = 0.0;
	sum_w = 0.0;
	
	pix_count = 0;

	// This is slow. Refine.
	min_row = floor((x0-xTL-radius)/ps[0]);
	max_row = ceil((x0-xTL+radius)/ps[0]);
	min_col = floor((yTL-y0-radius)/ps[1]);
	max_col = ceil((yTL-y0+radius)/ps[1]);
	
	if(min_row < 0) min_row = 0;
	if(min_col < 0) min_col = 0;
	if(max_row > nx) max_row = nx;
	if(max_col > ny) max_col = ny;
	
	// Add weights from the Gaussian!
	for(i_row=min_row; i_row<max_row; i_row++){
		for(i_col=min_col; i_col<max_col; i_col++){
		
			pn = i_col*nx + i_row;
			
			dist = sqrt(pow(xTL+i*ps[0]-x0,2.0)+pow(yTL-j*ps[1]-y0,2.0));
		
			if(dist <= radius){
				//weight = exp(-0.5*pow(dist/radius,2.0));
				//sum += weight*ndvi[i];
				sum += ndvi[i];
				//sum_w += weight;
				pix_count++;
			}
		}
	}
	out(pt,4) = sum / (double)pix_count;
	//out(pt,4) = sum / sum_w;
	*/
  }
	return out;
}

// [[Rcpp::export]]
NumericVector measure_blob_ndvi(List Blobs, int nx, int ny, NumericVector ndvi, NumericVector ps, double xTL, double yTL){
	
	int i,j,k,pn;

	double x0,y0,radius,sum;
	int pix_count;

	long int i_row, i_col;
	long int min_row,max_row,min_col,max_col;
	long int di;
	long int dj;

	double weight;
	double sum_w;
	double dist;
	
	NumericVector xarr = Blobs["x"];
	NumericVector yarr = Blobs["y"];
	NumericVector darr = Blobs["d"];

	int nblobs = xarr.size();
	NumericVector out(nblobs);


	for(int blob_id=0; blob_id<nblobs; blob_id++){
		x0 = xarr[blob_id];
		y0 = yarr[blob_id];
		radius = 0.5 * darr[blob_id];
		
		i = round((x0-xTL)/ps[0]);
		j = round((yTL-y0)/ps[1]);
		pn = i+ j*nx;
		
		// Measure ndvi at the central point:
		//out[blob_id] = ndvi[pn];

		// Measure ndvi at the entire circle
		di = ceil(radius/ps[0]);
		dj = ceil(radius/ps[1]);

		min_col = i - di;
		max_col = i+di;
		min_row = j-dj;
		max_row = j+dj;
		
		if(min_row < 0) min_row = 0;
		if(min_col < 0) min_col = 0;
		if(max_row > nx) max_row = ny;
		if(max_col > ny) max_col = nx;
		
		//ndvi
		sum = 0.0;
		sum_w = 0.0;
		pix_count = 0;

		for(i_row=min_row; i_row<max_row; i_row++){
			for(i_col=min_col; i_col<max_col; i_col++){
			
				pn = i_row*nx + i_col;
				
				dist = sqrt(pow(xTL+i_col*ps[0]-x0,2.0)+pow(yTL-i_row*ps[1]-y0,2.0));
			
				if(dist <= radius){
					//weight = exp(-0.5*pow(dist/radius,2.0));
					weight = 1.0;
					sum += weight*ndvi[pn];
					//sum += ndvi[pn];
					sum_w += weight;
					//pix_count++;
				}
			}
		}

		//out[blob_id] = sum / (double)pix_count;
		out[blob_id] = sum / sum_w;

	}
	
	return out;
}


// [[Rcpp::export]]
NumericVector measure_blob_ndvi_dens(List Blobs, int nx, int ny, NumericVector ndvi, NumericVector ps, double xTL, double yTL, double ndvi_thresh){
	
	int i,j,k,pn;

	double x0,y0,radius,sum;
	int pix_count, veg_pix_count;

	long int i_row, i_col;
	long int min_row,max_row,min_col,max_col;
	long int di;
	long int dj;

	double weight;
	double sum_w;
	double dist;
	
	NumericVector xarr = Blobs["x"];
	NumericVector yarr = Blobs["y"];
	NumericVector darr = Blobs["d"];

	int nblobs = xarr.size();
	NumericVector out(nblobs);


	for(int blob_id=0; blob_id<nblobs; blob_id++){
		x0 = xarr[blob_id];
		y0 = yarr[blob_id];
		radius = 0.5 * darr[blob_id];
		
		i = round((x0-xTL)/ps[0]);
		j = round((yTL-y0)/ps[1]);
		pn = i+ j*nx;
		
		// Measure ndvi at the central point:
		//out[blob_id] = ndvi[pn];

		// Measure ndvi at the entire circle
		di = ceil(radius/ps[0]);
		dj = ceil(radius/ps[1]);

		min_col = i - di;
		max_col = i+di;
		min_row = j-dj;
		max_row = j+dj;
		
		if(min_row < 0) min_row = 0;
		if(min_col < 0) min_col = 0;
		if(max_row > nx) max_row = ny;
		if(max_col > ny) max_col = nx;
		
		//ndvi
		sum = 0.0;
		sum_w = 0.0;
		pix_count = 0;
		veg_pix_count = 0;

		for(i_row=min_row; i_row<max_row; i_row++){
			for(i_col=min_col; i_col<max_col; i_col++){
			
				pn = i_row*nx + i_col;
				
				dist = sqrt(pow(xTL+i_col*ps[0]-x0,2.0)+pow(yTL-i_row*ps[1]-y0,2.0));
			
				if(dist < radius){
					if(ndvi[pn]>=ndvi_thresh) veg_pix_count++;
					pix_count++;
				}
			}
		}

		//out[blob_id] = sum / (double)pix_count;
		out[blob_id] = (double)veg_pix_count / (double)pix_count;
	}

	return out;
}



// [[Rcpp::export]]
NumericVector measure_blob_prob(List Blobs, int nx, int ny, NumericMatrix data, NumericVector ps, double xTL, double yTL, NumericVector mu, NumericMatrix Cov, double prob_thresh){

	int Nb = data.ncol();
	
	int i,j,k,l,pn;

	double x0,y0,radius,sum, sum2;
	NumericVector DN(Nb), res(Nb);
	int pix_count, veg_pix_count;

	long int i_row, i_col;
	long int min_row,max_row,min_col,max_col;
	long int di;
	long int dj;

	double weight;
	double sum_w;
	double dist;

	NumericVector xarr = Blobs["x"];
	NumericVector yarr = Blobs["y"];
	NumericVector darr = Blobs["d"];

	int nblobs = xarr.size();
	NumericVector out(nblobs);

	for(int blob_id=0; blob_id<nblobs; blob_id++){
		x0 = xarr[blob_id];
		y0 = yarr[blob_id];
		radius = 0.5 * darr[blob_id];
		
		i = round((x0-xTL)/ps[0]);
		j = round((yTL-y0)/ps[1]);
		pn = i+ j*nx;
		
		// Measure ndvi at the central point:
		//out[blob_id] = ndvi[pn];

		// Measure ndvi at the entire circle
		di = ceil(radius/ps[0]);
		dj = ceil(radius/ps[1]);

		min_col = i - di;
		max_col = i+di;
		min_row = j-dj;
		max_row = j+dj;
		
		if(min_row < 0) min_row = 0;
		if(min_col < 0) min_col = 0;
		if(max_row > nx) max_row = ny;
		if(max_col > ny) max_col = nx;
		
		sum = 0.0;
		pix_count = 0;

		for(i_row=min_row; i_row<max_row; i_row++){
			for(i_col=min_col; i_col<max_col; i_col++){
			
				pn = i_row*nx + i_col;
				
				dist = sqrt(pow(xTL+i_col*ps[0]-x0,2.0)+pow(yTL-i_row*ps[1]-y0,2.0));
			
				if(dist < radius){
					for(k=0;k<Nb;k++){
						DN[k] = data(pn,k);
						res[k] = DN[k]-mu[k];
					}
				
					sum2 = 0.0;
				
					for(k=0;k<Nb;k++){
						for(l=0;l<Nb;l++){
							sum2+=res[k]*Cov(k,l)*res[l];
						}
					}

					//sum +=  exp(-0.5*sum2);
					if(sum2 >= prob_thresh) sum += 1.0;
					//sum +=  sum2;

					pix_count++;
				}
			}
		}

		out[blob_id] = sum / (double)pix_count;
	}

	return out;
}




// [[Rcpp::export]]
List measure_blob_ndvi_profile(List Blob, int nx, int ny, NumericVector ndvi, NumericVector ps, double xTL, double yTL){
	
	int i,j,k,pn;

	double x0,y0,radius,sum;
	int pix_count, veg_pix_count;

	long int i_row, i_col;
	long int min_row,max_row,min_col,max_col;
	long int di;
	long int dj;

	double weight;
	double sum_w;
	double dist;
	
	x0 = Blob["x"];
	y0 = Blob["y"];
	radius = Blob["d"];
	radius *=0.5;

	List out;
	
	int npix;
	npix = ceil(M_PI*pow(radius,2.0)/(ps[0]*ps[1]));
	
	NumericVector r(npix);
	NumericVector ndvi_r(npix);

	i = round((x0-xTL)/ps[0]);
	j = round((yTL-y0)/ps[1]);
	pn = i+ j*nx;
	
	// Measure ndvi at the central point:
	//out[blob_id] = ndvi[pn];

	// Measure ndvi at the entire circle
	di = ceil(radius/ps[0]);
	dj = ceil(radius/ps[1]);

	min_col = i - di;
	max_col = i+di;
	min_row = j-dj;
	max_row = j+dj;
	
	if(min_row < 0) min_row = 0;
	if(min_col < 0) min_col = 0;
	if(max_row > nx) max_row = ny;
	if(max_col > ny) max_col = nx;
	
	//ndvi
	pix_count = 0;

	for(i_row=min_row; i_row<max_row; i_row++){
		for(i_col=min_col; i_col<max_col; i_col++){
		
			pn = i_row*nx + i_col;
			
			dist = sqrt(pow(xTL+i_col*ps[0]-x0,2.0)+pow(yTL-i_row*ps[1]-y0,2.0));
		
			if(dist < radius){
				r[pix_count] = dist;
				ndvi_r[pix_count] = ndvi[pn];
				pix_count++;
			}
		}
	}

		
	out["r"] = r;
	out["ndvi_r"] = ndvi_r;

	return out;
}




// [[Rcpp::export]]
NumericMatrix remove_duplicates(NumericMatrix A){
   
   int n_row = A.nrow(), n_col = A.ncol();
   int count=0,nmatch;
   double x,y,r,xa,ya,ra,dist;

   NumericMatrix B(n_row, n_col);
   
   for (int i = 0; i < n_row; i++){
		x = A(i,0);
		y = A(i,1);
		r = A(i,2);

		nmatch=0;

		for(int j = 0; j < i; j++){
			xa=A(j,0);
			ya=A(j,1);
			ra=A(j,2);
			
			dist = pow(x-xa,2.0)+pow(y-ya,2.0)+pow(r-ra,2.0);
			
			if(dist<=3.0){
				nmatch++;
			}
		}

		if(nmatch==0){
			for(int k = 0; k<n_col; k++){
				B(count,k)=A(i,k);
			}
			count++;
		}
   }

   NumericMatrix out(count, n_col);

	for(int i = 0; i < count; i++){
		for(int k = 0; k < n_col; k++){
				out(i,k)=B(i,k);
			}
	}
	
   return out;
}


// Computes determinant of Hessian matrix (or Laplacian) for a fixed value of scale t

// [[Rcpp::export]]
NumericVector det_HessianC(NumericVector A, int n_col, int n_row, double t){

	double sigma = sqrt(t);

	if(t==0){
		return A;
	}

	int	kern_size = ceil(6.0*sigma);

	if(kern_size<3){
		kern_size = 3;
	}

	if(kern_size%2 ==0) kern_size++;

	NumericMatrix W(kern_size,kern_size);
	NumericMatrix Wxx(kern_size,kern_size);
	NumericMatrix Wxy(kern_size,kern_size);
	NumericMatrix Wyy(kern_size,kern_size);

	int nx = kern_size;
	int ny = kern_size;

	int i0 = (nx-1)/2;
	int j0 = (ny-1)/2;

	double gain = 0.0;
	double dx, dy;

	for(int i=0; i<kern_size; i++){
		dx = i - i0;
		for(int j=0; j<kern_size; j++){
			dy = j - j0;
			W(i,j) = exp(-0.5*((dx*dx+dy*dy)/(sigma*sigma)));
			Wxx(i,j) = W(i,j) * ((dx*dx)/t - 1.0) / t;
			Wxy(i,j) = W(i,j) * dx*dy/(t*t);
			Wyy(i,j) = W(i,j) * ((dy*dy)/t - 1.0) / t;
			gain += W(i,j);
		}
	}

	for(int i=0; i<kern_size; i++){
		for(int j=0; j<kern_size; j++){
			W(i,j) 	 /= gain;
			Wxx(i,j) /= gain;
			Wxy(i,j) /= gain;
			Wyy(i,j) /= gain;
		}
	}

	NumericVector Gxx = convolve_2d_v2(A, n_col, n_row, Wxx);
	NumericVector Gxy = convolve_2d_v2(A, n_col, n_row, Wxy);
	NumericVector Gyy = convolve_2d_v2(A, n_col, n_row, Wyy);

	NumericVector out(n_row*n_col);

	for(int pn=0; pn<n_row*n_col; pn++){
		out(pn) = (t*t)*(Gxx[pn]*Gyy[pn] - Gxy[pn]*Gxy[pn]);
	}

	return out;
}


// [[Rcpp::export]]
NumericVector DOH_spaceC(NumericVector A, int n_col, int n_row, NumericVector tarr){

	int ns = tarr.size();

	int i_row, i_col, k, offset;

	NumericVector detHes(n_row * n_col * ns);

	for(k=0; k<ns; k++){
		offset = k*n_row*n_col;
		NumericVector DH = det_HessianC(A, n_col, n_row, tarr[k]);
		for(i_row=0; i_row<n_row; i_row++){
			for(i_col=0; i_col<n_col; i_col++){
				detHes[offset+i_row*n_col+i_col] = DH(i_row*n_col+i_col);
			}
		}
	}
	return detHes;
}

