#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector tree_mask_from_blobs(NumericMatrix blobs, int n_col, int n_row, double xTL, double yTL, NumericVector ps, double buffer){

	int nblobs = blobs.nrow(), npix = n_col*n_row;

	double pix_size = 0.5*(ps[0]+ps[1]);

	int i, j, di, dj, pn, i_row, i_col;
	int min_row, max_row, min_col, max_col;

	double x0, y0, r0, d;

	IntegerVector out(npix);
	
	for(pn=0; pn<npix; pn++){
		out[pn] = 0;
	}

	for(int id=0; id<nblobs; id++){
		x0 = blobs(id,0);
		y0 = blobs(id,1);
		r0 = blobs(id,2);

		i = round((x0-xTL)/ps[0]);
		j = round((yTL-y0)/ps[1]);

		di = ceil(r0/ps[0]);
		dj = ceil(r0/ps[1]);

		min_col = i - (di + buffer + 1);
		max_col = i + (di + buffer + 1);
		min_row = j - (dj + buffer + 1);
		max_row = j + (dj + buffer + 1);

		if(min_row < 0) min_row = 0;
		if(min_col < 0) min_col = 0;
		if(max_row > n_row) max_row = n_row;
		if(max_col > n_col) max_col = n_col;

		for(i_row=min_row; i_row<max_row; i_row++){
			for(i_col=min_col; i_col<max_col; i_col++){
			
				pn = i_row*n_col + i_col;
				if(out[pn]==1) continue;
				
				d = sqrt(pow(xTL+i_col*ps[0]-x0,2.0)+pow(yTL-i_row*ps[1]-y0,2.0));
			
				if(d <= r0+buffer*pix_size){
					out[pn] = 1;
				}
			}
		}

	}

	return out;
}

