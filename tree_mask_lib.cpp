#include <Rcpp.h>
using namespace Rcpp;

List uv_tangentC(double,double, double, double);
double project_pollock_quantitativeC(double, double, double, double, double);
NumericVector project_pollock_quantitative_matrixC(NumericVector, NumericVector, NumericVector, double, NumericVector);
NumericMatrix project_pollockC(double, double, double, double, double, double, double, double, double);
IntegerVector tree_mask_from_blobs(NumericMatrix, int, int, double, double, NumericVector, double);

// [[Rcpp::export]]
List uv_tangentC(double n, double cx, double theta, double af){
	
	double der0, vmax, v, u, min_diff;
	double dudv;
	double ut, vt;
	List out;

	if(af<1.0e-06){
		vt = 1.0;
		ut = 0.0;
	
		out["ut"] = ut;
		out["vt"] = vt;
		return out;
	}
	
	der0 = -1.0/(af*tan(theta));

	vmax = sqrt(1.0-pow(cx,2));
	
	if(vmax<1.e-16){
		out["ut"] = 0.0;
		out["vt"] = 0.0;

		return out;
	}

	min_diff = 1.0e+32;  // huge number

	vt = -1.0;	//default value, used if no better match is found

	long int nsteps = floor(vmax/0.001);
	NumericVector der(nsteps);
	long int i;
	double diff;

	for(i=0; i<nsteps; i++){

		v = i*0.001;
		
		u = (1.0 - pow(cx*cx + v*v,n/2.0));
		u = pow(u,1.0/n);
		dudv = -v * pow(1.0-pow(u,n),1.0-2.0/n) * pow(u,1.0-n);
		
		der[i] = dudv;

		diff = dudv-der0;
		if(diff<0) diff=-diff;
		
		if(diff<min_diff){
			min_diff = diff;
			vt = v;
		}
	}
	ut = pow(1.0-pow(cx*cx+vt*vt,n/2.0),1.0/n);

	out["ut"] = ut;
	out["vt"] = vt;
	
	//out["der0"] = der0;
	//out["der"]  = der;
	//out["vmax"] = vmax;
	//out["min_diff"]=min_diff;

	return out;
}


// [[Rcpp::export]]
double project_pollock_quantitativeC(double h1, double h2, double R, double theta, double n){

	double af = h2/R;
	double z_tan, r_tan;
	double top_model;
	//double Dt_model;

	double cx = 0.0;
	List tmp = uv_tangentC(n, cx, theta, af);

	z_tan = tmp["ut"];
	z_tan *= h2;
	r_tan = tmp["vt"];
	r_tan *= R;

	top_model = 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta));
	//Dt_model <- r_tan+R+z_tan*tan(theta)

	return top_model;
}

// [[Rcpp::export]]
NumericMatrix transform_xyC(NumericMatrix xy, NumericVector xy0, double alpha){
	
	long int npix=xy.nrow();
	long int i;
	double x0,y0;
	
	double x, y, phi, r;
	
	NumericMatrix out(npix,2);

	x0 = xy0[0];
	y0 = xy0[1];
	
	for(i=0; i<npix; i++){
		x = xy(i,0);
		y = xy(i,1);

		phi = atan2(y,x);
		r   = sqrt(x*x+y*y);

		phi += alpha;

		out(i,0) = x0 + r*cos(phi); 
		out(i,1) = y0 + r*sin(phi);
	}

	return out;
}

// [[Rcpp::export]]
NumericMatrix project_pollockC(double h1, double h2, double R, double x0, double y0, double theta, double alpha, double n, double pix){

	double af = h2/R;
	double z_tan, r_tan, top_model;
	double x, y, cx = 0.0;

	int n_shad = ceil((2.0*R)/pix) + 1;
	double step = (R*2.0)/((double) (n_shad-1));
	List tmp;

	//NumericVector xarr(n_shad), yarr(n_shad);
	NumericMatrix xyarr(n_shad,2);
	NumericVector xy0(2);

	xy0[0] = x0;
	xy0[1] = y0;

	for(int i=0; i<n_shad; i++){
		xyarr(i,0) = -R + ((double) i) * step;
		xyarr(i,1) = h1*tan(theta) - sqrt(-xyarr(i,0)*xyarr(i,0)+R*R);
	}

	NumericMatrix xy_shad_lower = transform_xyC(xyarr,xy0,alpha-3.0*M_PI/2.0);
	
	for(int i=0; i<n_shad; i++){
		x = +R - ((double)i) * step;
		cx = x/R;
		tmp = uv_tangentC(n,cx,theta,af);
		z_tan = tmp["ut"];
		z_tan *= h2;
		r_tan = tmp["vt"];
		r_tan *= sqrt(R*R-x*x);

		xyarr(i,0) = x;
		xyarr(i,1) = r_tan + (h1+z_tan)*tan(theta);
	}

	NumericMatrix xy_shad_upper = transform_xyC(xyarr,xy0,alpha-3.0*M_PI/2.0);

	NumericMatrix xy_shad(xy_shad_lower.nrow()+xy_shad_upper.nrow(),2);

	for(int i=0; i<n_shad; i++){
		xy_shad(i,0)=xy_shad_lower(i,0);
		xy_shad(i,1)=xy_shad_lower(i,1);
	}
	
	for(int i=0; i<n_shad; i++){
		xy_shad(n_shad+i,0)=xy_shad_upper(i,0);
		xy_shad(n_shad+i,1)=xy_shad_upper(i,1);
	}

	return xy_shad;
}

// [[Rcpp::export]]
NumericVector project_pollock_quantitative_matrixC(NumericVector h1, NumericVector h2, NumericVector R, double theta, NumericVector n){

	long int n_row = h1.size();
	long int i;
	
	NumericVector top_model(n_row);
	
	for(i=0; i<n_row; i++){
		top_model[i] = project_pollock_quantitativeC(h1[i], h2[i], R[i], theta, n[i]);
	}

	return top_model;
}

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


// [[Rcpp::export]]
NumericMatrix shadow_contour(double h1, double h2, double R, double x0, double y0, double theta, double alpha, double n, int n_sectors){

	double af = h2/R;
	double z_tan, r_tan, top_model;
	double x, y, cx = 0.0;

	int n_shad = n_sectors;
	double step = (R*2.0)/((double) n_shad-1);
	List tmp;

	//NumericVector xarr(n_shad), yarr(n_shad);
	NumericMatrix xyarr(n_shad,2);
	NumericVector xy0(2);
	xy0[0] = x0;
	xy0[1] = y0;

	double phi;
	double d_phi = M_PI/((double)n_sectors);
	
	for(int i=0; i<n_sectors; i++){
		phi = M_PI + d_phi*((double) i);
		x = R*cos(phi);
		y = R*sin(phi);
		xyarr(i,0) = x;
		xyarr(i,1) = y + h1*tan(theta);
	}

	NumericMatrix xy_shad_lower = transform_xyC(xyarr,xy0,alpha-3.0*M_PI/2.0);
	
	for(int i=0; i<n_sectors; i++){
		phi = d_phi*((double) i);
		x = R*cos(phi);
		y = R*sin(phi);

		cx = x/R;
		tmp = uv_tangentC(n,cx,theta,af);
		z_tan = tmp["ut"];
		z_tan *= h2;
		r_tan = tmp["vt"];
		r_tan *= y;

		xyarr(i,0) = x;
		xyarr(i,1) = r_tan + (h1+z_tan)*tan(theta);
	}

	NumericMatrix xy_shad_upper = transform_xyC(xyarr,xy0,alpha-3.0*M_PI/2.0);

	NumericMatrix xy_shad(xy_shad_lower.nrow()+xy_shad_upper.nrow()+1,2);

	for(int i=0; i<n_shad; i++){
		xy_shad(i,0)=xy_shad_lower(i,0);
		xy_shad(i,1)=xy_shad_lower(i,1);
	}
	
	for(int i=0; i<n_shad; i++){
		xy_shad(n_shad+i,0)=xy_shad_upper(i,0);
		xy_shad(n_shad+i,1)=xy_shad_upper(i,1);
	}
	xy_shad(2*n_shad,0) = xy_shad(0,0);
	xy_shad(2*n_shad,1) = xy_shad(0,1);
	
	return xy_shad;
}
