#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix geo_transformC(NumericMatrix, double, double);
List uv_tangentC(double,double, double, double);
double project_pollock_quantitativeC(double, double, double, double, double);
NumericVector project_pollock_quantitative_matrixC(NumericVector, NumericVector, NumericVector, double, NumericVector);
NumericMatrix transform_xyC(NumericMatrix, NumericVector, double);
NumericMatrix project_pollockC(double, double, double, double, double, double, double, double, double);

NumericMatrix geo_transformC(NumericMatrix A, double dx, double dy){
	long int n = A.nrow(), ncol=A.ncol();
	NumericMatrix out(n,ncol);

	for(long int i=0; i<n; i++){
		out(i,0) = A(i,0) + dx;
		out(i,1) = A(i,1) + dy;

		for(long int j=2; j<ncol; j++){
			out(i,j) = A(i,j);
			out(i,j) = A(i,j);
		}
	}

	return out;
}

// [[Rcpp::export]]
NumericMatrix match_candC(NumericMatrix A, NumericMatrix B, double dx, double dy, double r_search){
	long int na = A.nrow();
	long int nb = B.nrow();

	NumericMatrix out(na*nb,2);
	long int i, j, count;
	double xa, ya, ra;
	double xb, yb, rb;
	double rmin, rmax, ratio, dist;
	
	B = geo_transformC(B, dx, dy);

	count=0;
	
	for(i=0; i<na; i++){
		xa = A(i,0);
		ya = A(i,1);
		ra = 0.5*A(i,2);

		for(j=0; j<nb ; j++){
			xb = B(j,0);
			yb = B(j,1);
			rb = 0.5*B(j,2);

			dist = sqrt(pow(xb - xa,2) + pow(yb - ya,2));

			if(dist<=r_search){
				if(ra>rb){
					rmin = rb;
					rmax = ra;
				} else{
					rmin = ra;
					rmax = rb;
				}

				ratio = rmin/rmax;
				
				if(ratio>=0.9){
					out(count,0)= i+1;
					out(count,1)= j+1;
					count++;
				}
			}
		}
	}

	NumericMatrix out_short(count,2);

	for(long int k=0; k<count; k++){
		out_short(k,0) = out(k,0);
		out_short(k,1) = out(k,1);
	}      

	return out_short;
}

// [[Rcpp::export]]
NumericMatrix optimizer_geomC(NumericMatrix A, NumericMatrix B, NumericVector sarr, NumericMatrix matched_pairs, double req_precision, long int Nsamp){

	long int i,j,l,iter,s,n_cons;
	long int count=0;
	double res,sum;
	double xa,ya,xb,yb,dx,dy;
	NumericMatrix samp(Nsamp,5), Btmp(B.nrow(),B.ncol());
	NumericVector residuals(matched_pairs.nrow());

	for(iter=0;iter<Nsamp;iter++){
	//iter=0;
		
		s = sarr[iter]-1;

		i = matched_pairs(s,0)-1;
		j = matched_pairs(s,1)-1;

		// Estimate model parameters
		xa = A(i,0);
		ya = A(i,1);
		
		xb = B(j,0);
		yb = B(j,1);

		dx = xa-xb;
		dy = ya-yb;

		Btmp = geo_transformC(B,dx,dy);
		
		n_cons=0;
		sum = 0.0;
			
		for(l=0;l<matched_pairs.nrow();l++){
			//if(l!=s){
				i = matched_pairs(l,0)-1;
				j = matched_pairs(l,1)-1;

				xa = A(i,0);
				ya = A(i,1);
				
				xb = Btmp(j,0);
				yb = Btmp(j,1);

				res = sqrt(pow(xa-xb,2)+pow(ya-yb,2));
				residuals[l]=res;
				
				if(res<req_precision){
					
					sum+=res;
					n_cons++;
				}
			//}
		}

		if(n_cons>0){
			samp(count,0)=s+1;
			samp(count,1)=n_cons;
			samp(count,2)=dx;
			samp(count,3)=dy;
			samp(count,4)=sum/(double)n_cons;

			count++;
		}
	}

	return samp;
}

// [[Rcpp::export]]
NumericMatrix select_matchC(NumericMatrix A, NumericMatrix B, NumericMatrix matched_pairs, double dx, double dy, double req_precision){

	long int i,j,l;
	long int count=0;
	double res;
	double xa,ya,xb,yb;
	NumericVector val(matched_pairs.nrow());

	B = geo_transformC(B,dx,dy);

	for(l=0; l<matched_pairs.nrow(); l++){
		i = matched_pairs(l,0)-1;
		j = matched_pairs(l,1)-1;

		xa = A(i,0);
		ya = A(i,1);

		xb = B(j,0);
		yb = B(j,1);

		res = sqrt(pow(xa-xb,2)+pow(ya-yb,2));

		if(res<req_precision){
			val[count] = l;
			count++;
		}
	
	}

	NumericMatrix select_match(count,2);
	
	for(l=0; l<count; l++){
		select_match(l,0) = matched_pairs(val[l],0);
		select_match(l,1) = matched_pairs(val[l],1);
	}

	return select_match;
}

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
	
	List tmp;

	//List tmp = uv_tangentC(n, cx, theta, af);
	//z_tan = tmp["ut"];
	//z_tan *= h2;
	//r_tan = tmp["vt"];
	//r_tan *= R;

	//top_model = 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta));
	//Dt_model <- r_tan+R+z_tan*tan(theta)

	int n_shad = ceil((2.0*R)/pix) + 1;
	double step = (R*2.0)/((double) n_shad-1);

	//NumericVector xarr(n_shad), yarr(n_shad);
	NumericMatrix xyarr(n_shad,2);
	NumericVector xy0(2);
	
	xy0[0] = x0;
	xy0[1] = y0;

	for(int i=0; i<n_shad; i++){
		//xarr[i] = -R + i * step;
		//yarr[i] = h1*tan(theta) - sqrt(-xarr[i]*xarr[i]+R*R);
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

		y = r_tan + (h1+z_tan)*tan(theta);
		xyarr(i,0) = x;
		xyarr(i,1) = y;
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
NumericMatrix project_shad_mask(double h1, double h2, double R, double x0, double y0, double theta, double alpha, double n, double pix, double buffer){

	double af = h2/R;
	double z_tan, r_tan, top_model;
	double x, y, cx = 0.0;
	
	List tmp;

	int n_shad = ceil((2.0*R)/pix) + 1;
	double step = (R*2.0)/((double) n_shad-1);

	//NumericVector xarr(n_shad), yarr(n_shad);
	NumericMatrix xyarr(n_shad,2);
	NumericVector xy0(2);
	
	xy0[0] = x0;
	xy0[1] = y0;

	for(int i=0; i<n_shad; i++){
		//xarr[i] = -R + i * step;
		//yarr[i] = h1*tan(theta) - sqrt(-xarr[i]*xarr[i]+R*R);
		xyarr(i,0) = -R + ((double) i) * step;
		xyarr(i,1) = h1*tan(theta) - sqrt(-xyarr(i,0)*xyarr(i,0)+R*R) - buffer;
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

		y = r_tan + (h1+z_tan)*tan(theta);
		xyarr(i,0) = x;
		xyarr(i,1) = y + buffer;
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
List eval_shadC(NumericVector params, double x0, double y0, double R, double theta_sat, double alpha_sat, double theta_sun, double alpha_sun, NumericMatrix xy_shad_pix, NumericVector ps_pan){

	const double h0 = 30.0;

	int n_shad_pix = xy_shad_pix.nrow();

	double h  = params[0]*h0;
	double hf = params[1];
	double n  = params[2];

	double pix = 0.5*(ps_pan[0]+ps_pan[1]);
	
	double h1, h2, top_shift;
	double x1, y1, x2, y2, r, phi;

	NumericVector x(n_shad_pix),y(n_shad_pix);

	h1 = h*hf;
	h2 = h - h1;

	top_shift = project_pollock_quantitativeC(h1,h2,R,theta_sat,n);
	x1 = x0 + top_shift * cos(alpha_sat);
	y1 = y0 + top_shift * sin(alpha_sat);


	for(int i=0; i<n_shad_pix; i++){
		x[i] = xy_shad_pix(i,0) - x1;
		y[i] = xy_shad_pix(i,1) - y1;

		r = sqrt(x[i]*x[i]+y[i]*y[i]);
		phi = atan2(y[i],x[i]);
		phi -= alpha_sun-3.0*M_PI/2.0;
		x[i] = r*cos(phi);
		y[i] = r*sin(phi);
	}

	NumericMatrix shad_pol = project_pollockC(h1,h2,R,0,0,theta_sun,3*M_PI/2,n,0.25);

	// Get pixels inside the projected shadow region
	int n_out = 0, n_in =0;
	int n_all = n_shad_pix;

	for(int i=0; i<n_shad_pix; i++){
		if((x[i]< -R) || (x[i] > R)) n_out ++;
	}

	int n_shad = ceil((2.0*R)/pix) + 1;
	double step = (R*2.0)/((double) n_shad-1);

	double shad_area = 0.0;

	for(int i=0; i<n_shad; i++){
		x1 = -R + ((double) i) * step;
		x2 = x1 + step;
		
		y1 = (shad_pol(i,1) + shad_pol(i+1,1))/2.0;
		y2 = (shad_pol(i+n_shad,1)+shad_pol(i+1+n_shad,1))/2.0;

		if(x[i]>=x1 && x[i]<x2){
			if(y[i]>=y1 && y[i] <=y2){
				n_in++;
			}else{
				n_out++;
			}
		}

		shad_area += (x2-x1)*(y2-y1);
	}

	shad_area /= (ps_pan[0]*ps_pan[1]);
	
	double val = 2.0 * n_in - n_out - 2.0 * (shad_area - n_in);
	
	// normalize by area of the shadow mask
	val /= n_shad_pix;

	List out;
	out["shad_pol"]=shad_pol;
	out["n_shad_pix"]=n_shad_pix;
	out["ts"]=top_shift;
	out["x1"]=x1;
	out["y1"]=y1;
	
	out["x"]=x;
	out["y"]=y;

	out["val"] = val;
	out["n_in"] = n_in;
	out["n_out"] = n_out;
	out["shad_area"] = shad_area;

	return out;
	
	
//	return val;

}


// [[Rcpp::export]]
NumericMatrix project_shadowC(double h1, double h2, double R, double x0, double y0, double theta, double alpha, double n, double pix, double buffer){

	double af = h2/R;
	double z_tan, r_tan, top_model;
	double x, y, cx = 0.0;
	
	List tmp;

	//List tmp = uv_tangentC(n, cx, theta, af);
	//z_tan = tmp["ut"];
	//z_tan *= h2;
	//r_tan = tmp["vt"];
	//r_tan *= R;

	//top_model = 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta));
	//Dt_model <- r_tan+R+z_tan*tan(theta)

	int n_shad = ceil((2.0*R)/pix) + 1;
	double step = (R*2.0)/((double) n_shad-1);

	//NumericVector xarr(n_shad), yarr(n_shad);
	NumericMatrix xyarr(n_shad,2);
	NumericVector xy0(2);
	
	xy0[0] = x0;
	xy0[1] = y0;

	for(int i=0; i<n_shad; i++){
		//xarr[i] = -R + i * step;
		//yarr[i] = h1*tan(theta) - sqrt(-xarr[i]*xarr[i]+R*R);
		xyarr(i,0) = -R + ((double) i) * step;
		xyarr(i,1) = h1*tan(theta) - sqrt(-xyarr(i,0)*xyarr(i,0)+R*R) - buffer;
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

		y = r_tan + (h1+z_tan)*tan(theta) + buffer;
		xyarr(i,0) = x;
		xyarr(i,1) = y;
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
