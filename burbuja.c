/*

	INSTALACION DE FFTW3
		./configure --prefix="/home/$USER/.fftw3" --enable-threads --enable-openmp --enable-shared=yes
		make
		sudo make install


	gcc -I$HOME/.fftw3/include burbuja.c -L$HOME/.fftw3/lib -fopenmp -lfftw3_omp -lfftw3 -lm -o bubbles
	gcc -I/mnt/MD1200B/hernandez/crosales/mpi-examples/include fare.c -L/mnt/MD1200B/hernandez/crosales/mpi-examples/lib -fopenmp -lfftw3_omp -lfftw3 -lm -o fare

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#include <fftw3.h>
#include <sys/time.h>
#include "1Dprint.h"
#include <sys/stat.h>
#include <sys/types.h>

int Qv_u(double [],int ,double ,double );
int Qv_w(double [],int ,double ,double );
int Fdqvdz(double [],int ,double ,double );

void MAX_Real(double *,double *,double *,double *,double *,double *,double *,
int ,int ,double ,double ,double ,double ,double );
void ABS_Real(double *,double *,double *,double *,double *,double *,double *,
int ,int ,double ,double ,double ,double ,double );
void MIN_Real(double *,double *,double *,double *,double *,double *,double *,
int ,int ,double ,double ,double ,double ,double );
void Fluct(double *,double *,double *,double *,double *,int ,int ,double ,
double ,double );
void Thermo(double *,double *,double *,double *,double *,int ,double ,double );
void MAX_Cplx(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,int ,int );
void MAX_Flux(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
int ,int );

void leerarchivo(char *,double *,int ,int ,int );
double max_vecs(double [], int );
double max_vecp(double *,int );
double max_vecp2(double *,int ,int );
double min_vecp(double *,int );
double abs_vecp(double *,int );
double max_complexvar(fftw_complex *,int );
double max_qt_level(double *,int ,int );
double min_qt_level(double *,int ,int );
double max_varhat_level(fftw_complex *,int ,int );
int varprima(double *,int ,int );

void Imprime1D(double *,int ,int ,int ,double );
void Imprime2D(double *,int ,int ,int ,int ,int ,double );
void Imprime3D(double *,int ,int ,int ,int ,int ,int ,double );
void Impresion1D(double *,double *,int ,int ,double ,double );
void Impresion2D(double *,double *,double *,double *,double *,
double *,int ,int ,int ,int ,double ,double ,double );
void Impresion3D(double *,double *,double *,double *,double *,
double *,int ,int ,int ,int ,int ,double ,double ,double );
int mat(int ,int ,int ,int ,int );
void AbrirArchivos(int );
double timeval_diff(struct timeval *a,struct timeval *b)
{
	return (double)(a->tv_sec + (double)a->tv_usec/1000000) -
							 (double)(b->tv_sec + (double)b->tv_usec/1000000);
}
int upad_real(double *,double *,int ,int ,int ,int ,int ,int );
int upad_cplx0(fftw_complex *,fftw_complex *,int ,int ,int ,int ,int ,int );
int upad_cplx(fftw_complex *,fftw_complex *,int ,int ,int ,int ,int ,int ,int );
int pad_real(double *,double *,int ,int ,int ,int ,int ,int );
int pad_cplx(fftw_complex *,fftw_complex *,int ,int ,int ,int ,int ,int );
void f_print(double *,int ,int );
double promedio(double *,int ,int );
double minimum(double ,double ,double );
double min(double ,double );
double max(double ,double );
double promedioa(double *,int );
double suma(double *,double *,double ,int ,int ,int );

void moisture(double *,double *,double *,double [],double [],int ,int );
void NoProm(double *,int ,int );

const int NUM_X = 128;
const int NUM_Y = 128;
const int NUM_Z = 101;
const int NUM_Y2 = 65;

// Linear Part
int LinFlux(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,double *,int ,int ,
double ,double ,double ,double );

int Fourier(int ,double *,double *,double *,fftw_plan ,double *,double *,int [],
int [],fftw_complex *);
int dzbar(int ,int ,double *, double *,double ,double []);

int NonFlux(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
double *,double *,double *,double *,double *,double *,fftw_complex *,
double *,double *,double *,double , double ,double ,double ,int [],int [],
fftw_complex *,fftw_complex *,fftw_plan );

int FluxH1(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
int ,int ,double ,double ,double );

int FluxH23(fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
fftw_complex *,fftw_complex *,fftw_complex *,fftw_complex *,
int ,int ,double ,double ,double ,double );

int Thomas(fftw_complex *,fftw_complex *,fftw_complex *,int ,int ,double ,
double ,double ,double ,double *);


int main(int argc, char *argv[])
{
	/*
	El código se divide en las siguientes partes:
	- Definir y dar valor a los parámetros
	- Crear las variables a utilizar
	- Crear los perfiles promedio de saturación
	- Generar las condiciones iniciales de las variables (espacio físico)
	- Generar las condiciones iniciales de las variables (espacio de fourier)
		- Correr el Runge-Kutta
		- Actualizar después de cada Runge-Kutta
		- Al término del 3er RK, actualizar las variables termodinámicas
		- Modificar los perfiles promedio
	
	Agosto 30
	Editar datos sobre burbujas y perfiles segun Majda (2002)
	*/

	int m,nx_phy,ny_phy,nx_fou,ny_fou,nx_upd,ny_upd;
	int savedata2D,savedata3D,fileindex,bt,pi_[5],bardata;
	double movies1D,movies2D,movies3D;
	double Lx,Ly,Lz,dx,dy,dz,drx,dry,zk;
	double Ti,Tfinal,dt,dt0,dtinv,nu0,nu1,CFL,mindt,max_theta,N_Th;
	int ix,iy,iz,ic,i,d_upd,d_phy,d_fou;
	_Complex IM=1.0*I;
	struct timeval t_ini,t_fin;
	double secs;

	// Grid space sizes
	m = NUM_Z;
	nx_phy = NUM_X;
	ny_phy = NUM_Y;
	nx_fou = NUM_X;
	ny_fou = NUM_Y2;
	double Th_bar[m];
	double Th_bar2[m-2];

	// Grid space padding sizes
	nx_upd = 2*floor(nx_phy/3)+1;
	ny_upd = ny_phy/3 + 1;
	d_upd = nx_upd*ny_upd;
	d_phy = nx_phy*ny_phy;
	d_fou = nx_fou*ny_fou;

	// Domain size
	Lx = 20.0*4.8154;	// Non-dimensionless
	Ly = 20.0*4.8154;	// Non-dimensionless
	Lz = 1.0;			// Non-dimensionless

	dx = Lx/(1.0*nx_phy); dy = Ly/(1.0*ny_phy);	dz = Lz/(1.0*(m-1));

	printf("0.5*Lz = %e\n",Lz*0.5);
	printf("dz = %e\n",dz);
	printf("coefficient = %e\n",0.5*Lz/dz);
	//return 0;
	double T0,PrT,Ra,Rvd,Q0,PrQ,qs,Vs,LcpT0,qvs0,qv0,mp,parTe,parqt,parqr,Gamma;
	T0 = 1.0;
	PrT = 7.0;
	Ra = 2.0;
	Rvd = 0.622;
	Q0 = 1.0;
	PrQ = 1.0;
	qs = 1.0;
	Vs = 1.0;
	LcpT0 = 25.0/3.0;
	qvs0 = 2.0;
	qv0 = 1.8;
	mp = 1.0;
	Gamma = Ra/PrT;
	parTe = 1.0;
	parqt = -(LcpT0 - Rvd)*qs;
	parqr = (LcpT0 - Rvd - 1.0)*qs;

	double alpha1,alpha2,alpha3;
	alpha1 = 29.0/96.0;
	alpha2 = -3.0/40.0;
	alpha3 = 1.0/6.0;
	
	double beta1,beta2,beta3;
	beta1 = 37.0/60.0;
	beta2 = 5.0/24.0;
	beta3 = 1.0/6.0;
	
	double gamma1,gamma2,gamma3;
	gamma1 = 8.0/15.0;
	gamma2 = 5.0/12.0;
	gamma3 = 3.0/4.0;
	
	double eta1,eta2;
	eta1 = -17.0/60.0;
	eta2 = -5.0/12.0;

	double Re,PeT,PeQ;
	Re = 1.0;
	PeT = PrT;
	PeQ = PrQ;

	// Final time
	Tfinal = 20.0; // In hours
	
	// Frames
	movies1D = 2000;
	movies2D = 200; // In frames per hour
	movies3D = 2; 	// In frames per hour

	// Wave numbers are multipli of those
	drx = 2.0*M_PI/Lx;
	dry = 2.0*M_PI/Ly;

	pi_[0] = nx_upd;
	pi_[1] = ny_upd;
	pi_[2] = d_upd;
	pi_[3] = d_phy;
	pi_[4] = d_fou;

	int file = 1;
	char carpet[100];
	sprintf(carpet,"Run_%d",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Run_%d/Data1D",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Run_%d/Data2D",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Run_%d/Data3D",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Run_%d/Videos",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Run_%d/Parameters",file);
	mkdir(carpet,0777);
	//AbrirArchivos(file);

	//--------------------------------------------------------------------------
	// FFTW OPENMP initialize
	int nthreads;
	fftw_init_threads();
	nthreads = omp_get_max_threads();
 	fftw_plan_with_nthreads(nthreads);
  	printf("Using %u OpenMP threads.\n",nthreads);

	// rin and cout creation
	double *rin;
	rin = fftw_alloc_real(d_phy);
	fftw_complex *cout;
	cout=fftw_alloc_complex(d_fou);

	// FFTW forward and backward plans creation
	fftw_plan planf,planb;
	planf = fftw_plan_dft_r2c_2d(nx_phy,ny_phy,rin,cout,FFTW_PATIENT);
	planb = fftw_plan_dft_c2r_2d(nx_phy,ny_phy,cout,rin,FFTW_PATIENT);
	printf("passing fftw openmp directrices\n");

	//--------------------------------------------------------------------------
	// Variables definitions
	// grid 1
	double *u,*v,*o;
	u = fftw_alloc_real(d_phy*(m+1));
	v = fftw_alloc_real(d_phy*(m+1));
	o = fftw_alloc_real(d_phy*(m+1));

	// grid 2
	double *w,*Th,*Te,*qv,*qr,*qt;
	 w = fftw_alloc_real(d_phy*m);
	Th = fftw_alloc_real(d_phy*m);
	Te = fftw_alloc_real(d_phy*m);
	qv = fftw_alloc_real(d_phy*m);
	qr = fftw_alloc_real(d_phy*m);
	qt = fftw_alloc_real(d_phy*m);

	double array_aux[3];
	// L1,L2 norms
	double *dt_l1,*dt_l2;
  	dt_l1 = fftw_alloc_real(d_phy*m);
  	dt_l2 = fftw_alloc_real(d_phy*m);

	// Vanishing viscosity variables
	double *enu1,*enu2,*enu3;
	enu1 = fftw_alloc_real(d_upd);
	enu2 = fftw_alloc_real(d_upd);
	enu3 = fftw_alloc_real(d_upd);

	double *enu1_h,*enu2_h,*enu3_h;
	enu1_h = fftw_alloc_real(d_upd);
	enu2_h = fftw_alloc_real(d_upd);
	enu3_h = fftw_alloc_real(d_upd);

	double *enu1_w,*enu2_w,*enu3_w;
	enu1_w = fftw_alloc_real(d_upd);
	enu2_w = fftw_alloc_real(d_upd);
	enu3_w = fftw_alloc_real(d_upd);

	// Complex variables
	int dim_u = d_upd*(m+1);
	int dim_w = d_upd*m;

	// original variables complex
	fftw_complex *uh,*vh,*oh;
	uh = fftw_alloc_complex(dim_u);
	vh = fftw_alloc_complex(dim_u);
	oh = fftw_alloc_complex(dim_u);
	
	fftw_complex *wh,*qrh,*boh;//,*Teh,*qth;
	wh = fftw_alloc_complex(dim_w);
	qrh = fftw_alloc_complex(dim_w);
	boh = fftw_alloc_complex(dim_w);
	
	// psi variables
	fftw_complex *V10,*V11,*V12,*L1,*N10,*N11,*N12;
	V10 = fftw_alloc_complex(dim_u);
	V11 = fftw_alloc_complex(dim_u);
	V12 = fftw_alloc_complex(dim_u);
	 L1 = fftw_alloc_complex(dim_u);
	N10 = fftw_alloc_complex(dim_u);
	N11 = fftw_alloc_complex(dim_u);
	N12 = fftw_alloc_complex(dim_u);	
	
	// phi variables
	fftw_complex *V20,*V21,*V22,*L2,*N20,*N21,*N22;
	V20 = fftw_alloc_complex(dim_w);
	V21 = fftw_alloc_complex(dim_w);
	V22 = fftw_alloc_complex(dim_w);
	 L2 = fftw_alloc_complex(dim_w);
	N20 = fftw_alloc_complex(dim_w);
	N21 = fftw_alloc_complex(dim_w);
	N22 = fftw_alloc_complex(dim_w);

	// Theta_e variables
	fftw_complex *V30,*V31,*V32,*L3,*N30,*N31,*N32;
	V30 = fftw_alloc_complex(dim_w);
	V31 = fftw_alloc_complex(dim_w);
	V32 = fftw_alloc_complex(dim_w);
	 L3 = fftw_alloc_complex(dim_w);
	N30 = fftw_alloc_complex(dim_w);
	N31 = fftw_alloc_complex(dim_w);
	N32 = fftw_alloc_complex(dim_w);

	// Q_t variables
	fftw_complex *V40,*V41,*V42,*L4,*N40,*N41,*N42;
	V40 = fftw_alloc_complex(dim_w);
	V41 = fftw_alloc_complex(dim_w);
	V42 = fftw_alloc_complex(dim_w);
	 L4 = fftw_alloc_complex(dim_w);
	N40 = fftw_alloc_complex(dim_w);
	N41 = fftw_alloc_complex(dim_w);
	N42 = fftw_alloc_complex(dim_w);

	fftw_complex *H1,*H2,*H3,*H4;
	H1 = fftw_alloc_complex(dim_u);
	H2 = fftw_alloc_complex(dim_w);
	H3 = fftw_alloc_complex(dim_w);
	H4 = fftw_alloc_complex(dim_w);

	printf("passing variables creation\n");
	//--------------------------------------------------------------------------
	// padding definitions
	int vp_[6];
	int piece = floor(nx_phy/3);
	vp_[0] = piece + 1;
	vp_[1] = 0;
	vp_[2] = 0;
	vp_[3] = piece;
	vp_[4] = nx_phy - piece;
	vp_[5] = piece + 1;
	printf("passing unpadding definitions\n");
	for(i=0;i<6;i++) printf("vp_[%d] = %d\n",i,vp_[i]);
	//return 0;

	//--------------------------------------------------------------------------
	// fourier derivatives definition
	// kx,ky,kk unpadding
	double *kxup,*kyup,*kkup;
	kyup = fftw_alloc_real(d_fou);
	kxup = fftw_alloc_real(d_fou);
	kkup = fftw_alloc_real(d_fou);

	// kx,ky,kk padding
	double *kx,*ky,*kk,*kkinv;
	ky = fftw_alloc_real(d_upd);
	kx = fftw_alloc_real(d_upd);
	kk = fftw_alloc_real(d_upd);
	kkinv = fftw_alloc_real(d_upd);

	for(ix=0;ix<nx_fou;++ix){
		for(iy=0;iy<ny_fou;++iy){
			i = mat(nx_fou,ny_fou,ix,iy,0);
			kyup[i] = 1.0*iy;
		}
	}

	int sep;
	sep = nx_fou/2 + 1;
	for(ix=0;ix<sep;++ix){
		for(iy=0;iy<ny_fou;++iy){
			i = mat(nx_fou,ny_fou,ix,iy,0);
			kxup[i] = 1.0*ix;
		}
	}for(ix=sep;ix<nx_fou;++ix){
		for(iy=0;iy<ny_fou;++iy){
			i = mat(nx_fou,ny_fou,ix,iy,0);
			kxup[i] = 1.0*(ix-nx_fou);
		}
	}

	for(i=0;i<d_fou;++i){
		kxup[i] = drx*kxup[i];
		kyup[i] = dry*kyup[i];
		kkup[i] = pow(kxup[i],2.0)+pow(kyup[i],2.0);
	}

	upad_real(kx,kxup,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],0);
	upad_real(ky,kyup,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],0);
	upad_real(kk,kkup,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],0);

	upad_real(kx,kxup,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],0);
	upad_real(ky,kyup,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],0);
	upad_real(kk,kkup,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],0);

	kkinv[0] = 0.0;
	for(i=1;i<d_upd;++i) kkinv[i] = 1.0/kk[i];
	printf("passing fourier derivatives\n");

	free(kkup);
	free(kyup);
	free(kxup);

	//--------------------------------------------------------------------------
	// Thermodynamic parameters
	double qvs[m],qt_bg[m];
	Qv_w(qvs,m,dz,qvs0);
	Qv_w(qt_bg,m,dz,qv0);
	//for(iz=0;iz<m;iz++) printf("qvs = %e, qv_bg = %e at %d\n",qvs[iz],qv_bg[iz],iz);
	//--------------------------------------------------------------------------
	// Initial data
	for(i=0;i<d_phy*m;i++){
		 w[i] = 0.0;
		qv[i] = 0.0;
		qr[i] = 0.0;
		qt[i] = 0.0;
		Th[i] = 0.0;
		Te[i] = 0.0;
	}for(i=0;i<d_phy*(m+1);++i){
		u[i] = 0.0;
		v[i] = 0.0;
		o[i] = 0.0;
	}
	//printf("stop the calls1\n");

	//--------------------------------------------------------------------------
	double dphy;
	dphy = 1.0*d_phy;
	double *r_num;
	r_num = fftw_alloc_real(d_phy*m);
	srand(time(NULL));
	for(i=0;i<d_phy*m;++i) r_num[i] = rand()/((double) RAND_MAX + 1.0) - 0.5;
	// numeros aleatorios no normalizados, con media no zero para cada Z
	//printf("stop the calls2\n");

	//--------------------------------------------------------------------------
	// temperature perturbation (for Th fluctuations as in FARE)
	double aux_Th,aux_qt;
	for(iz=0;iz<m;++iz){
		zk = iz*dz;
		for(ic=0;ic<d_phy;++ic){
			i = ic + iz*d_phy;
			Th[i] = Th[i] + 2.0*mp*r_num[i];
			if(zk < 0.5*Lz) qt[i] = qt[i] + 0.1*2.0*r_num[i];
			//Th[i] = Th[i] + Th_bg[iz]; // Perturbation + steady state
		}//NoProm(Th,d_phy,iz);
		aux_Th = max_qt_level(Th,d_phy,iz);
		aux_qt = max_qt_level(qt,d_phy,iz);
		//printf("max Th = %e, max qt = %e, at %d\n",aux_Th,aux_qt,iz);
	}free(r_num);
	//printf("stop the calls4\n");

	//--------------------------------------------------------------------------

	// Initialization of complex variables
	// Each variable needs to initiate with mean 0
	for(i=0;i<d_phy;++i) rin[i] = 0.0;
	for(i=0;i<d_fou;++i) cout[i] = 0.0 + 0.0*IM;

	for(i=0;i<d_upd;++i){
		for(iz=0;iz<m;++iz){
			wh[i+d_upd*iz] = 0.0 + 0.0*IM;
			qrh[i+d_upd*iz] = 0.0 + 0.0*IM;
			V20[i+d_upd*iz] = 0.0 + 0.0*IM;
			V30[i+d_upd*iz] = 0.0 + 0.0*IM;
			V40[i+d_upd*iz] = 0.0 + 0.0*IM;			
		}for(iz=0;iz<m+1;++iz){
			uh[i+d_upd*iz] = 0.0 + 0.0*IM;
			vh[i+d_upd*iz] = 0.0 + 0.0*IM;
			oh[i+d_upd*iz] = 0.0 + 0.0*IM;
			V10[i+d_upd*iz] = 0.0 + 0.0*IM;
		}
	}
	//printf("stop the calls zero\n");

	// moist bubbles (for qt fluctuation as in FARE)
	double qt_tot,Th_tot,Te_tot;
	for(iz=0;iz<m+1;++iz){
		// u to u0hat
		for(i=0;i<d_phy;++i) rin[i] = u[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) uh[i] = 0.0 + 0.0*IM;
		upad_cplx(uh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(uh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);

		// v to v0hat
		for(i=0;i<d_phy;++i) rin[i] = v[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) vh[i] = 0.0 + 0.0*IM;
		upad_cplx(vh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(vh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);

		// vorticity and streamfunction
		for(i=0;i<d_upd;i++){
			oh[i+d_upd*iz] = IM*kx[i]*vh[i+d_upd*iz]
							- IM*ky[i]*uh[i+d_upd*iz];
			V10[i+d_upd*iz] = -kkinv[i]*oh[i+d_upd*iz];
		}oh[iz*m] = 0.0 + 0.0*IM;
		V10[iz*m] = 0.0 + 0.0*IM;
		
		// o0hat to o
		for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
		pad_cplx(cout,oh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
		pad_cplx(cout,oh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
		fftw_execute(planb);
		for(i=0;i<d_phy;++i) o[i+d_phy*iz] = rin[i];
	}for(iz=0;iz<m;++iz){
		// w to w0hat
		for(i=0;i<d_phy;++i) rin[i] = w[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) wh[i] = 0.0 + 0.0*IM;
		upad_cplx(wh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(wh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
		wh[iz*m] = 0.0 + 0.0*IM;
		for(i=0;i<d_upd;i++) V20[i+d_upd*iz] = -kkinv[i]*wh[i+d_upd*iz];
		V20[iz*m] = 0.0 + 0.0*IM;
		for(i=0;i<d_upd;++i) wh[i+d_upd*iz] = -kk[i]*V20[i+d_upd*iz];
		for(i=0;i<d_fou;++i) cout[i] = 0.0 + 0.0*IM;
		pad_cplx(cout,wh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
		pad_cplx(cout,wh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
		fftw_execute(planb);
		for(i=0;i<d_phy;++i) w[i+d_phy*iz] = rin[i];
		
		NoProm(qt,d_phy,iz);
		NoProm(Th,d_phy,iz);
		for(ic=0;ic<d_phy;ic++){
			i = ic + d_phy*iz;
			qt_tot = qt[i] + qt_bg[iz]; // perturbation + steady state
			Th_tot = Th[i];// + Th_bg[iz]; // perturbation + steady state
			if(qt_tot - qvs[iz] > 0.0)
				qr[i] = qt_tot - qvs[iz]; // rain water
			else qr[i] = 0.0;
			qv[i] = qt_tot - qr[i];
			if(qv[i] < 0.0) qv[i] = 0.0;
			Te[i] = Th_tot + LcpT0*Q0*qv[i];
			qt[i] = qv[i] + qr[i];
		}

		// Te to Tehat
		for(i=0;i<d_phy;++i) rin[i] = Te[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) V30[i] = 0.0 + 0.0*IM;
		upad_cplx(V30,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(V30,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
		V30[iz*m] = 0.0 + 0.0*IM;
		for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
		pad_cplx(cout,V30,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
		pad_cplx(cout,V30,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
		fftw_execute(planb);
		for(i=0;i<d_phy;++i) Te[i+d_phy*iz] = rin[i];

		// qt to qthat
		for(i=0;i<d_phy;++i) rin[i] = qt[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) V40[i] = 0.0 + 0.0*IM;
		upad_cplx(V40,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(V40,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
		V40[iz*m] = 0.0 + 0.0*IM;
		
		// qthat to qt
		for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
		pad_cplx(cout,V40,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
		pad_cplx(cout,V40,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
		fftw_execute(planb);
		for(i=0;i<d_phy;++i) qt[i+d_phy*iz] = rin[i];

		for(ic=0;ic<d_phy;ic++){
			i = ic + d_phy*iz;
			qt_tot = qt[i] + qt_bg[iz]; // perturbation + steady state
			if(qt_tot - qvs[iz] > 0.0)
				qr[i] = qt_tot - qvs[iz]; // rain water
			else qr[i] = 0.0;
			Th[i] = Te[i] - LcpT0*Q0*(qt[i]-qr[i]);
		}

		// qr to qrhat
		for(i=0;i<d_phy;++i) rin[i] = qr[i+d_phy*iz];
		fftw_execute(planf);
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) qrh[i] = 0.0 + 0.0*IM;
		upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
		upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
		qrh[iz*m] = 0.0 + 0.0*IM;
		
		// qrhat to qr
		for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
		pad_cplx(cout,qrh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
		pad_cplx(cout,qrh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
		fftw_execute(planb);
		for(i=0;i<d_phy;++i) qr[i+d_phy*iz] = rin[i];
		
		for(i=d_upd*iz;i<d_upd*(iz+1);++i) boh[i] = parTe*V30[i] + parqt*V40[i] - parqr*qrh[i];
		//aux_Th = max_qt_level(Te,d_phy,iz);
		//aux_qt = max_qt_level(qt,d_phy,iz);
		//printf("max Te = %e, max qt = %e, ",aux_Th,aux_qt);
		//aux_qt = max_qt_level(qr,d_phy,iz);
		//printf("max qr = %e at %d\n",aux_qt,iz);
	}


	//----------------------------------------------------------------------
	//							INITIAL CONDITIONS
	//----------------------------------------------------------------------
	printf("\n	         PRINTING IN 1D\n");
	//Impresion1D(qt_bg,Te_bg,m+1,file,Ths,qs);

	printf("\n	         PRINTING IN 2D\n");
	Impresion2D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,1.0,1.0,1.0);

	printf("\n	         PRINTING IN 3D\n");
	//Impresion3D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,fileindex,Us,Ths,qs);

	// Saving data
	   bardata = 1;
	 fileindex = 1;
	savedata2D = 1;

	bt = 0;
	Ti = 0.0;
	array_aux[0] = 1.0;//pow(M_PI/(10.0*9.0),-1.0);

	// Time iteration start
	Ti = 0.0;
	dt0 = 0.144;// 2.5*pow(10.0,-3.0);

	CFL = 0.01;
	mindt = 1.0;
	max_theta = 0.0;
	dt = dt0;
	//return 0;
	double ui,vi,wi;
	double wtime = 0.0;
	int it_large = 0;
	int cont = 0;
	
	printf("stop the calls before\n");
	//return 0;
	do{
		printf("\n\n\n\n-------------------------------------------------\n\n");
		printf("	      At step %d  \n\n\n",bt);
		//printf("   Elapsed time = %4.11f seg\n",Ti*Ts*60.0*60.0);
		//printf("   Resting time = %4.11f seg\n",(Tfinal-Ti)*Ts*60.0*60.0);

		//----------------------------------------------------------------------
		//							TIME STEP
		//----------------------------------------------------------------------
		wtime = omp_get_wtime();
		for(iz=0;iz<m;++iz) Th_bar[iz] = promedio(Th,d_phy,iz);
		for(iz=0;iz<m-3;++iz) Th_bar2[iz] = 0.5*fabs(Th_bar[iz+3] - Th_bar[iz+1])/dz;
		N_Th = sqrt(81.0 + 8.1*max_vecs(Th_bar2,m-3));
		dt = dt0;

		#pragma omp parallel for private(iz,ic,ui,vi,wi) shared(dt_l1,dt_l2)
		for(i=0;i<d_phy;++i){
			for(iz=0;iz<m;++iz){
				ic = i + d_phy*iz;
				ui = 0.5*(u[ic] + u[ic + d_phy])/dx;
				vi = 0.5*(v[ic] + v[ic + d_phy])/dy;
				wi = w[ic]/dz;
				dt_l1[ic] = fabs(ui) + fabs(vi) + fabs(wi);
				dt_l2[ic] = sqrt(pow(ui,2.0) + pow(vi,2.0) + pow(wi,2.0));
			}
		}array_aux[1] = max_vecp(dt_l1,d_phy*m);
		//ImprimirPar(file,16,array_aux[1]);
		array_aux[2] = max_vecp(dt_l2,d_phy*m);
		//ImprimirPar(file,17,array_aux[2]);
		dtinv = max_vecs(array_aux,3);
		//printf("             V1 = %.15e\n",dt*Ts*60.0*60.0);
		//printf("             V2 = %.15e\n",CFL/dtinv*Ts*60.0*60.0);
		//printf("             V3 = %.15e\n",CFL*M_PI/(10.0*N_Th)*Ts*60.0*60.0);

		//dt = CFL/dtinv;
		dt = minimum(dt,CFL/dtinv,M_PI/(10.0*N_Th));
		bt = bt + 1;
		Ti = Ti + dt;		//ImprimirPar(file,3,Ti*Ts);

		printf("           N_th = %.15e\n",N_Th);
		printf("             dt = %.15e\n",dt);//*Ts*60.0*60.0);
		//ImprimirPar(file,0,dt*Ts);

		//----------------------------------------------------------------------
		// 						FIRST STEP OF TIME-EVOLUTION
		//----------------------------------------------------------------------
		//for(iz=0;iz<m;iz++) printf("max qr at %d = %e\n",iz,max_qt_level(qr,d_phy,iz));
		LinFlux(L1,L2,L3,L4,V10,V20,V30,V40,kk,d_upd,m,Re,PeT,PeQ,dz);

		NonFlux(N10,N20,N30,N40,u,v,o,w,Te,qt,V20,kx,ky,kk,
		Gamma,Vs/dz,PeT,PeQ,pi_,vp_,qrh,boh,planf);
		//printf("nonlinflux ready\n");

		FluxH1(H1,H2,H3,H4,V10,V20,V30,V40,L1,L2,L3,L4,N10,N20,N30,N40,
		d_upd,m,alpha1,gamma1,dt);
		double aux_H1,aux_H2;
		for(iz=0;iz<m;iz++){
			aux_H1 = max_varhat_level(H1,d_upd,iz);
			aux_H2 = max_varhat_level(H2,d_upd,iz);
			printf("max H1 = %e, max H2 = %e\n",aux_H1,aux_H2);
		}
		//printf("flux ready\n");

		// actualizar phi
		Thomas(H1,H2,V21,d_upd,m,beta1,dz,dt,Re,kk);
		//printf("thomas ready\n");
		for(i=0;i<d_upd;++i){
			for(iz=1;iz<m;iz++){
				ic = i + d_upd*iz;
				V11[ic] = (1.0/(1.0 - beta1*dt*(-kk[i])/Re))*(H1[ic] 
									+ (beta1*dt/dz)*(V21[ic] - V21[ic-d_upd]));
			}
		}

		// actualizar Te, qt
		for(i=0;i<d_upd;++i){
			for(iz=0;iz<m;++iz){
				ic = i + d_upd*iz;
				V31[ic] = (1.0/(1.0 - beta1*dt*(-kk[i])/PeT))*H3[ic];
				V41[ic] = (1.0/(1.0 - beta1*dt*(-kk[i])/PeQ))*H4[ic];				
			}
		}
		
		//----------------------------------------------------------------------
		// 							C2R FOURIER TRANSFORMS
		//----------------------------------------------------------------------
		for(iz=0;iz<m;++iz){
			for(i=0;i<d_upd;++i) wh[i+d_upd*iz] = -kk[i]*V21[i+d_upd*iz];
			
			// w1hat to w
			for(i=0;i<d_fou;++i) cout[i] = 0.0 + 0.0*IM;
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) w[i+d_phy*iz] = rin[i];

			// Tehat to Te
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V31,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V31,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) Te[i+d_phy*iz] = rin[i];

			// qthat to qt
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V41,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V41,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) qt[i+d_phy*iz] = rin[i];

			// qr creation
			for(ic=0;ic<d_phy;ic++){
				i = ic + d_phy*iz;
				qt_tot = qt[i] + qt_bg[iz]; // perturbation + steady state
				if(qt_tot - qvs[iz] > 0.0)
					qr[i] = qt_tot - qvs[iz]; // rain water
				else qr[i] = 0.0;
			}
			for(i=0;i<d_phy;++i) rin[i] = qr[i+d_phy*iz];
			fftw_execute(planf);
			for(i=d_upd*iz;i<d_upd*(iz+1);++i) qrh[i] = 0.0 + 0.0*IM;
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
			qrh[iz*m] = 0.0 + 0.0*IM;

			for(i=d_upd*iz;i<d_upd*(iz+1);++i) boh[i] = parTe*V31[i] + parqt*V41[i] + parqr*qrh[i];
			
			printf("qt = %e,Te = %e, qr = %e, w = %e\n",max_qt_level(qt,d_phy,iz),max_qt_level(Te,d_phy,iz),max_qt_level(qr,d_phy,iz),max_qt_level(w,d_phy,iz));
		}for(iz=0;iz<m+1;++iz){
			for(i=0;i<d_upd;++i){
				uh[i+d_upd*iz] = -IM*ky[i]*V11[i+d_upd*iz];
				vh[i+d_upd*iz] =  IM*kx[i]*V11[i+d_upd*iz];
				oh[i+d_upd*iz] = -kk[i]*V11[i+d_upd*iz];
			}
			
			// u1hat to u
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) u[i+d_phy*iz] = rin[i];

			// v1hat to v
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) v[i+d_phy*iz] = rin[i];

			// o1hat to o
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) o[i+d_phy*iz] = rin[i];

			printf("u = %e,v = %e, o = %e\n",max_qt_level(u,d_phy,iz),max_qt_level(v,d_phy,iz),max_qt_level(o,d_phy,iz));
		}

		//----------------------------------------------------------------------
		// 						SECOND STEP OF TIME-EVOLUTION
		//----------------------------------------------------------------------
		LinFlux(L1,L2,L3,L4,V11,V21,V31,V41,kk,d_upd,m,Re,PeT,PeQ,dz);

		NonFlux(N11,N21,N31,N41,u,v,o,w,Te,qt,V21,kx,ky,kk,
		Gamma,Vs/dz,PeT,PeQ,pi_,vp_,qrh,boh,planf);

		FluxH23(H1,H2,H3,H4,V11,V21,V31,V41,L1,L2,L3,L4,N11,N21,N31,N41,
								N10,N20,N30,N40,d_upd,m,alpha2,gamma2,eta1,dt);

		// actualizar phi
		Thomas(H1,H2,V22,d_upd,m,beta2,dz,dt,Re,kk);
		for(i=0;i<d_upd;++i){
			for(iz=1;iz<m;iz++){
				ic = i + d_upd*iz;
				V12[ic] = (1.0/(1.0 - beta2*dt*(-kk[i])/Re))*(H1[ic] 
									+ (beta2*dt/dz)*(V22[ic] - V22[ic-d_upd]));
			}
		}

		// actualizar Te, qt
		for(i=0;i<d_upd;++i){
			for(iz=0;iz<m;++iz){
				ic = i + d_upd*iz;
				V32[ic] = (1.0/(1.0 - beta2*dt*(-kk[i])/PeT))*H3[ic];
				V42[ic] = (1.0/(1.0 - beta2*dt*(-kk[i])/PeQ))*H4[ic];				
			}
		}

		//----------------------------------------------------------------------
		// 							C2R FOURIER TRANSFORMS
		//----------------------------------------------------------------------
		for(iz=0;iz<m;++iz){
			for(i=0;i<d_upd;++i) wh[i+d_upd*iz] = -kk[i]*V22[i+d_upd*iz];
			
			// w1hat to w
			for(i=0;i<d_fou;++i) cout[i] = 0.0 + 0.0*IM;
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) w[i+d_phy*iz] = rin[i];

			// Tehat to Te
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V32,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V32,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) Te[i+d_phy*iz] = rin[i];

			// qthat to qt
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V42,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V42,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) qt[i+d_phy*iz] = rin[i];

			// qr creation
			for(ic=0;ic<d_phy;ic++){
				i = ic + d_phy*iz;
				qt_tot = qt[i] + qt_bg[iz]; // perturbation + steady state
				if(qt_tot - qvs[iz] > 0.0)
					qr[i] = qt_tot - qvs[iz]; // rain water
				else qr[i] = 0.0;
			}
			for(i=0;i<d_phy;++i) rin[i] = qr[i+d_phy*iz];
			fftw_execute(planf);
			for(i=d_upd*iz;i<d_upd*(iz+1);++i) qrh[i] = 0.0 + 0.0*IM;
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
			qrh[iz*m] = 0.0 + 0.0*IM;

			for(i=d_upd*iz;i<d_upd*(iz+1);++i) boh[i] = parTe*V32[i] + parqt*V42[i] + parqr*qrh[i];
		}for(iz=0;iz<m+1;++iz){
			for(i=0;i<d_upd;++i){
				uh[i+d_upd*iz] = -IM*ky[i]*V12[i+d_upd*iz];
				vh[i+d_upd*iz] =  IM*kx[i]*V12[i+d_upd*iz];
				oh[i+d_upd*iz] = -kk[i]*V12[i+d_upd*iz];
			}
			
			// u1hat to u
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) u[i+d_phy*iz] = rin[i];

			// v1hat to v
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) v[i+d_phy*iz] = rin[i];

			// o1hat to o
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) o[i+d_phy*iz] = rin[i];
		}

		//----------------------------------------------------------------------
		// 						THIRD STEP OF TIME-EVOLUTION
		//----------------------------------------------------------------------
		LinFlux(L1,L2,L3,L4,V12,V22,V32,V42,kk,d_upd,m,Re,PeT,PeQ,dz);

		NonFlux(N12,N22,N32,N42,u,v,o,w,Te,qt,V22,kx,ky,kk,
		Gamma,Vs/dz,PeT,PeQ,pi_,vp_,qrh,boh,planf);

		FluxH23(H1,H2,H3,H4,V12,V22,V32,V42,L1,L2,L3,L4,N12,N22,N32,N42,
								N11,N21,N31,N41,d_upd,m,alpha3,gamma3,eta2,dt);

		// actualizar phi
		Thomas(H1,H2,V20,d_upd,m,beta3,dz,dt,Re,kk);
		for(i=0;i<d_upd;++i){
			for(iz=1;iz<m;iz++){
				ic = i + d_upd*iz;
				V10[ic] = (1.0/(1.0 - beta3*dt*(-kk[i])/Re))*(H1[ic]
									 + (beta3*dt/dz)*(V20[ic] - V20[ic-d_upd]));
			}
			V10[i] = -V10[i+d_upd];
			V10[i+m*d_upd] = V10[i+(m-1)*d_upd];
		}

		// actualizar Te, qt
		for(i=0;i<d_upd;++i){
			for(iz=0;iz<m;++iz){
				ic = i + d_upd*iz;
				V30[ic] = (1.0/(1.0 - beta3*dt*(-kk[i])/PeT))*H3[ic];
				V40[ic] = (1.0/(1.0 - beta3*dt*(-kk[i])/PeQ))*H4[ic];				
			}
		}


		//----------------------------------------------------------------------
		// 							C2R FOURIER TRANSFORMS
		//----------------------------------------------------------------------
		for(iz=0;iz<m;++iz){
			for(i=0;i<d_upd;++i) wh[i+d_upd*iz] = -kk[i]*V20[i+d_upd*iz];
			
			// w1hat to w
			for(i=0;i<d_fou;++i) cout[i] = 0.0 + 0.0*IM;
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,wh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) w[i+d_phy*iz] = rin[i];

			// Tehat to Te
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V30,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V30,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) Te[i+d_phy*iz] = rin[i];

			// qthat to qt
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,V40,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,V40,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) qt[i+d_phy*iz] = rin[i];

			// qr creation
			for(ic=0;ic<d_phy;ic++){
				i = ic + d_phy*iz;
				qt_tot = qt[i] + qt_bg[iz]; // perturbation + steady state
				if(qt_tot - qvs[iz] > 0.0)
					qr[i] = qt_tot - qvs[iz]; // rain water
				else qr[i] = 0.0;
				Th[i] = qt[i] - LcpT0*Q0*(qt[i]-qr[i]);
			}
			for(i=0;i<d_phy;++i) rin[i] = qr[i+d_phy*iz];
			fftw_execute(planf);
			for(i=d_upd*iz;i<d_upd*(iz+1);++i) qrh[i] = 0.0 + 0.0*IM;
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],d_phy,iz);
			upad_cplx(qrh,cout,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],d_phy,iz);
			qrh[iz*m] = 0.0 + 0.0*IM;

			for(i=d_upd*iz;i<d_upd*(iz+1);++i) boh[i] = parTe*V30[i] + parqt*V40[i] + parqr*qrh[i];
		}for(iz=0;iz<m+1;++iz){
			for(i=0;i<d_upd;++i){
				uh[i+d_upd*iz] = -IM*ky[i]*V10[i+d_upd*iz];
				vh[i+d_upd*iz] =  IM*kx[i]*V10[i+d_upd*iz];
				oh[i+d_upd*iz] = -kk[i]*V10[i+d_upd*iz];
			}
			
			// u1hat to u
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,uh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) u[i+d_phy*iz] = rin[i];

			// v1hat to v
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,vh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) v[i+d_phy*iz] = rin[i];

			// o1hat to o
			for(i=0;i<d_fou;++i) cout[i]=0.0 + 0.0*IM;
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[0],vp_[1],vp_[2],iz);
			pad_cplx(cout,oh,nx_upd,ny_upd,vp_[3],vp_[4],vp_[5],iz);
			fftw_execute(planb);
			for(i=0;i<d_phy;++i) o[i+d_phy*iz] = rin[i];
		}

		//----------------------------------------------------------------------
		// 							SAVE DATA 2D & 3D (Complete variables)
		//----------------------------------------------------------------------
		wtime = omp_get_wtime()-wtime;
		if(Ti>=(1.0*savedata2D*Tfinal/movies2D)){
			printf("\n	         PRINTING IN 2D\n");
			Impresion2D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,1.0,1.0,1.0);
			if(savedata2D%12 == 0){
				printf("\n	         PRINTING IN 3D\n");
				//Impresion3D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,fileindex,Us,Ths,qs);
				fileindex=fileindex+1;
			}savedata2D=savedata2D+1;
		}
//		bt = bt + 1;
		//----------------------------------------------------------------------
		//
		//----------------------------------------------------------------------
		ImprimirPar(file,18,wtime);
		printf("Tiempo de respuesta paralelo2 :%f \n",wtime);
		//printf("     Relacion loop--time-step :%f \n",wtime/(Ts*dt*60.0*60.0));
		//ImprimirPar(file,4,wtime/(Ts*dt*60.0*60.0));
	}while(Ti<Tfinal && bt < 2); //*** end of while **************
	printf("stop the calls after\n");


	printf("\n	         PRINTING IN 1D\n");
	//Impresion1D(qv_bg,Te_bg,m+1,file,Ths,qs);

	printf("\n	         PRINTING IN 2D at %d\n",savedata2D);
	//Impresion2D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,Us,Ths,qs);

	printf("\n	         PRINTING IN 3D at %d\n",fileindex);
	//Impresion3D(u,v,w,qv,qr,Th,nx_phy,ny_phy,m,file,fileindex,Us,Ths,qs);

	fftw_destroy_plan(planf);
	fftw_destroy_plan(planb);

	free(u);
	free(v);
	free(o);
	free(w);
	free(Th);
	free(Te);

	free(qv);
	free(qr);
	free(qt);

	free(kx);
	free(ky);
	free(kk);
	free(kkinv);

	free(enu1);
	free(enu2);
	free(enu3);

	free(enu1_h);
	free(enu2_h);
	free(enu3_h);

	free(enu1_w);
	free(enu2_w);
	free(enu3_w);

	free(rin);
	fftw_free(cout);
	free(dt_l1);
	free(dt_l2);

	fftw_free(L1);
	fftw_free(V10);
	fftw_free(V11);
	fftw_free(V12);
	fftw_free(N10);
	fftw_free(N11);
	fftw_free(N12);

	fftw_free(L2);
	fftw_free(V20);
	fftw_free(V21);
	fftw_free(V22);
	fftw_free(N20);
	fftw_free(N21);
	fftw_free(N22);

	fftw_free(L3);
	fftw_free(V30);
	fftw_free(V31);
	fftw_free(V32);
	fftw_free(N30);
	fftw_free(N31);
	fftw_free(N32);

	fftw_free(L4);
	fftw_free(V40);
	fftw_free(V41);
	fftw_free(V42);
	fftw_free(N40);
	fftw_free(N41);
	fftw_free(N42);
}

int LinFlux(fftw_complex *L1,fftw_complex *L2,fftw_complex *L3,fftw_complex *L4,
fftw_complex *V1,fftw_complex *V2,fftw_complex *V3,fftw_complex *V4,double *kk,
int d,int m,double Re,double PeT,double PeQ,double dz)
{
	int i,iz,ic;
	_Complex IM = 0.0 + 1.0*I;
	#pragma omp parallel private(iz,ic) shared(L1,L2,L3,L4)
	for(i=0;i<d;++i){
		for(iz=1;iz<m;++iz){
			ic = i + d*iz;
			L1[ic] = (pow(-kk[i],2.0)/Re)*V1[ic] + 
						((-kk[i])/dz)*(V2[ic]-V2[ic-d]);
		}for(iz=1;iz<m-1;++iz){
			ic = i + d*iz;
			L2[ic] = (pow(-kk[i],2.0)/Re)*V2[ic] - 
						(1.0/dz)*(V1[ic+d] - V1[ic]);
			L3[ic] = ((-kk[i])/PeT)*V3[ic];
			L4[ic] = ((-kk[i])/PeQ)*V4[ic];
		}
		
		// Boundaries
		L1[i] = 0.0 + 0.0*IM;
		L2[i] = 0.0 + 0.0*IM;
		L3[i] = ((-kk[i])/PeT)*V3[i];
		L4[i] = ((-kk[i])/PeQ)*V4[i];

		L1[i + d*m] = 0.0 + 0.0*IM;
		L2[i + d*(m-1)] = 0.0 + 0.0*IM;
		L3[i + d*(m-1)] = ((-kk[i])/PeT)*V3[i + d*(m-1)];
		L4[i + d*(m-1)] = ((-kk[i])/PeQ)*V4[i + d*(m-1)];
	}return 0;
}

int Fourier(int m,double *Var,double *u,double *v,fftw_plan planf1,double *kx,
double *ky,int vp_[],int pi_[],fftw_complex *N)
{
	int iz,i,ic;
	_Complex IM=0.0 + 1.0*I;
	fftw_complex *cout;
	cout=fftw_alloc_complex(pi_[4]);
	double *rin;
	rin = fftw_alloc_real(pi_[3]);
	
	double *u_prod,*v_prod;
	u_prod = fftw_alloc_real(pi_[3]);
	v_prod = fftw_alloc_real(pi_[3]);

	fftw_complex *u_prodh,*v_prodh;
	u_prodh = fftw_alloc_complex(pi_[2]);
	v_prodh = fftw_alloc_complex(pi_[2]);

	#pragma omp parallel for shared(N)
	for(i=0;i<pi_[2]*m;++i){
		N[i] = 0.0 + 0.0*IM;
	}

	for(iz=0;iz<m;++iz){
		#pragma omp parallel for shared(u_prod,v_prod)
		for(i=0;i<pi_[3];++i){
			ic = i + pi_[3]*iz;
			u_prod[i] = u[ic]*Var[ic];
			v_prod[i] = v[ic]*Var[ic];
		}
		// u_var to u_varhat
		for(i=0;i<pi_[4];++i) cout[i] = 0.0 + 0.0*IM;
		for(i=0;i<pi_[3];++i)  rin[i] = u_prod[i];
		fftw_execute(planf1);
		upad_cplx0(u_prodh,cout,pi_[0],pi_[1],vp_[0],vp_[1],vp_[2],pi_[3]);
		upad_cplx0(u_prodh,cout,pi_[0],pi_[1],vp_[3],vp_[4],vp_[5],pi_[3]);

		// v_var to v_varhat
		for(i=0;i<pi_[4];++i) cout[i] = 0.0 + 0.0*IM;
		for(i=0;i<pi_[3];++i)  rin[i] = v_prod[i];
		fftw_execute(planf1);
		upad_cplx0(v_prodh,cout,pi_[0],pi_[1],vp_[0],vp_[1],vp_[2],pi_[3]);
		upad_cplx0(v_prodh,cout,pi_[0],pi_[1],vp_[3],vp_[4],vp_[5],pi_[3]);

		#pragma omp parallel for shared(N)
		for(i=0;i<pi_[2];++i){
			N[i+pi_[2]*iz] = - IM*kx[i]*u_prodh[i] - IM*ky[i]*v_prodh[i];
		}
	}
	
	free(u_prod);
	free(v_prod);
	
	fftw_free(u_prodh);
	fftw_free(v_prodh);
	return 0;
}

int dzbar(int m,int d,double *w,double *var,double Pr,double dzvar[])
{
	int i,iz;
	double wvarbar,wvar;
	wvarbar = 0.0;
	for(i=0;i<d*m;i++) wvarbar = wvarbar + var[i]*w[i];
	for(iz=0;iz<m;iz++){
		wvar = 0.0;
		for(i=0;i<d;i++) wvar = wvar + var[i+d*iz]*w[i+d*iz];
		dzvar[iz] = -1.0 + Pr*(m*wvar - wvarbar)/(d*m);
	}return 0;
}

int NonFlux(fftw_complex *N1,fftw_complex *N2,fftw_complex *N3,fftw_complex *N4,
double *u,double *v,double *o,double *w,double *Te,double *qt,fftw_complex *phi,
double *kx,double *ky,double *kk,double Gamma,double VTdz,double PrT,double PrQ,
int pi_[],int vp_[],fftw_complex *qr,fftw_complex *bo,fftw_plan planf1)
{
	double wtime = 0.0;
	wtime = omp_get_wtime();

	int i,iz,ic;
	_Complex IM = 0.0 + 1.0*I;
	
	int nx_phy,ny_phy,nx_fou,ny_fou,m;
	nx_phy = NUM_X;
	ny_phy = NUM_Y;
	nx_fou = NUM_X;
	ny_fou = NUM_Y2;
	m = NUM_Z;

	int nx_upd,ny_upd,d_upd,d_phy,d_fou;
	nx_upd = pi_[0];
	ny_upd = pi_[1];
	d_upd = pi_[2];
	d_phy = pi_[3];
	d_fou = pi_[4];

	#pragma omp parallel for shared(N1)
	for(i=0;i<d_upd*(m+1);++i){
		N1[i] = 0.0 + 0.0*IM;
	}

	#pragma omp parallel for shared(N2,N3,N4)
	for(i=0;i<d_upd*m;++i){
		N2[i] = 0.0 + 0.0*IM;
		N3[i] = 0.0 + 0.0*IM;
		N4[i] = 0.0 + 0.0*IM;
	}

	// Jacobianos
	Fourier(m+1, o,u,v,planf1,kx,ky,vp_,pi_,N1);
	Fourier( m , w,u,v,planf1,kx,ky,vp_,pi_,N2);
	Fourier( m ,Te,u,v,planf1,kx,ky,vp_,pi_,N3);
	Fourier( m ,qt,u,v,planf1,kx,ky,vp_,pi_,N4);

	// Derivadas de estados base
	double dzTe[m],dzqt[m];
	dzbar(m,d_phy,w,Te,PrT,dzTe);
	dzbar(m,d_phy,w,qt,PrQ,dzqt);

	#pragma omp parallel for shared(N2,N3,N4) private(iz,i)
	for(ic=0;ic<d_upd;++ic){
		for(iz=0;iz<m-1;iz++){
			i = ic + d_upd*iz;
			N2[i] = N2[i] + Gamma*bo[i];
			N3[i] = N3[i] - (-kk[ic])*dzTe[iz]*phi[i];
			N4[i] = N4[i] - (-kk[ic])*dzqt[iz]*phi[i] + VTdz*(qr[i+d_upd] - qr[i]);
		}iz=m-1;
		i = ic + d_upd*iz;
		N2[i] = N2[i] + Gamma*bo[i];
		N3[i] = N3[i] - (-kk[ic])*dzTe[iz]*phi[i];
		N4[i] = N4[i] - (-kk[ic])*dzqt[iz]*phi[i];
	}

	wtime = omp_get_wtime()-wtime;
	return 0;
}

int FluxH1(fftw_complex *H1,fftw_complex *H2,fftw_complex *H3,fftw_complex *H4,
fftw_complex *V1,fftw_complex *V2,fftw_complex *V3,fftw_complex *V4,
fftw_complex *L1,fftw_complex *L2,fftw_complex *L3,fftw_complex *L4,
fftw_complex *N1,fftw_complex *N2,fftw_complex *N3,fftw_complex *N4,
int d_upd,int m,double a,double g,double dt)
{
	int i,iz,ic;
	#pragma omp parallel private(iz,ic) shared(H1,H2,H3,H4)
	for(i=0;i<d_upd;++i){
		for(iz=0;iz<m+1;++iz){
			ic = i + d_upd*iz;
			H1[ic] = V1[ic] + dt*(a*L1[ic] + g*N1[ic]);
		}for(iz=0;iz<m;++iz){
			ic = i + d_upd*iz;
			H2[ic] = V2[ic] + dt*(a*L2[ic] + g*N2[ic]);
			H3[ic] = V3[ic] + dt*(a*L3[ic] + g*N3[ic]);
			H4[ic] = V4[ic] + dt*(a*L4[ic] + g*N4[ic]);			
		}
	}return 0;
}

int FluxH23(fftw_complex *H1,fftw_complex *H2,fftw_complex *H3,fftw_complex *H4,
fftw_complex *V1,fftw_complex *V2,fftw_complex *V3,fftw_complex *V4,
fftw_complex *L1,fftw_complex *L2,fftw_complex *L3,fftw_complex *L4,
fftw_complex *N1,fftw_complex *N2,fftw_complex *N3,fftw_complex *N4,
fftw_complex *NN1,fftw_complex *NN2,fftw_complex *NN3,fftw_complex *NN4,
int d_upd,int m,double a,double g,double e,double dt)
{
	int i,iz,ic;
	#pragma omp parallel private(iz,ic) shared(H1,H2,H3,H4)
	for(i=0;i<d_upd;++i){
		for(iz=0;iz<m+1;++iz){
			ic = i + d_upd*iz;
			H1[ic] = V1[ic] + dt*(a*L1[ic] + g*N1[ic] + e*NN1[ic]);
		}for(iz=0;iz<m;++iz){
			ic = i + d_upd*iz;
			H2[ic] = V2[ic] + dt*(a*L2[ic] + g*N2[ic] + e*NN2[ic]);
			H3[ic] = V3[ic] + dt*(a*L3[ic] + g*N3[ic] + e*NN2[ic]);
			H4[ic] = V4[ic] + dt*(a*L4[ic] + g*N4[ic] + e*NN4[ic]);			
		}
	}return 0;
}

int Thomas(fftw_complex *H1,fftw_complex *H2,fftw_complex *phi,int d_upd,int m,
double b,double dz,double dt,double Re,double *kk)
{
	//printf("inside Thomas 1\n");
	int i,ic,iz;
	double b1,b2;
	b2 = (b*dt/dz);
	printf("b2 = %e\n",b2);
	double A[m-2],B[m-2],C[m-2],CT[m-2];
	fftw_complex D[m-2],DT[m-2];
	A[0] = 0.0;
	for(iz=1;iz<m-2;iz++) A[iz] = -pow(b2,2.0);
	for(iz=0;iz<m-3;iz++) C[iz] = -pow(b2,2.0);
	C[m-3] = 0.0;
	//printf("all ingredients Thomas 2\n");
	for(ic=0;ic<d_upd;ic++){
		b1 = (1.0-b*dt*(-kk[ic])/Re);
		for(iz=0;iz<m-2;iz++) B[iz] = 2.0*pow(b2,2.0) + (-kk[ic])*pow(b1,2.0);
		for(iz=0;iz<m-2;iz++) D[iz] = (-kk[ic])*b1*H2[ic+d_upd*(iz+1)]
		 + b1*(H1[ic+d_upd*(iz+2)] - H1[ic+d_upd*(iz+1)]);
		CT[0] = C[0]/B[0];
		for(iz=1;iz<m-2;iz++) CT[iz] = C[iz]/(B[iz]-C[iz-1]*A[iz]);
		DT[0] = D[0]/B[0];
		for(iz=1;iz<m-2;iz++) DT[iz] = (D[iz] - D[iz-1]*A[iz])/(B[iz]-CT[iz-1]*A[iz]);
		phi[ic+d_upd*(m-2)] = DT[m-3];
		for(iz=m-4;iz>-1;iz--) phi[ic+d_upd*(iz+1)] = DT[iz] - CT[iz]*phi[ic+d_upd*(iz+2)];
	}return 0;
}

void NoProm(double *var,int d,int iz)
{
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<d;i++) sum = sum + var[i + d*iz];
	for(i=0;i<d;i++) var[i + d*iz] = var[i + d*iz] - sum/(1.0*d);
}

double suma(double *var,double *w,double dz,int iz,int nx,int ny)
{
	int ic,level;
	double vp[nx*ny],vm[nx*ny];
	double auxp = 0.0,auxm = 0.0;
	for(ic=0;ic<nx*ny;ic++){
		level = ic + nx*ny*iz;
		// terminar de modificar
		vp[ic] = var[level]*w[level];
		vm[ic] = var[level-nx*ny]*w[level-nx*ny];
	}
	auxp = promedioa(vp,nx*ny);
	auxm = promedioa(vm,nx*ny);
	return (auxp-auxm)/dz;
}

int varprima(double *var,int d_phy,int m)
{
	int ic,level,iz;
	double suma;
	for(iz=0;iz<m;iz++){
		suma = promedio(var,d_phy,iz);
		for(ic=0;ic<d_phy;ic++){
			level = ic + d_phy*iz;
			var[level] = var[level] - suma;
		}
	}return 0;
}


int mat(int L,int M,int x,int y,int z)
{
	return y+M*(x+L*z);
}

int upad_real(double *var_new,double *var_old,int nx_upd,int ny_upd,int len_x,
int st_x1,int st_x2,int iz)
{
	int ix,iy,i,j;
	for(ix=0;ix<len_x;++ix){
		for(iy=0;iy<ny_upd;++iy){
			i = mat(NUM_X,NUM_Y2,ix+st_x1,iy,0);
			j = mat(nx_upd,ny_upd,ix+st_x2,iy,iz);
			var_new[j] = var_old[i];
		}
	}return 0;
}

int upad_cplx0(fftw_complex *var_new,fftw_complex *var_old,int nx_upd,int ny_upd,
int len_x,int st_x1,int st_x2,int dim)
{
	int ix,iy,i,j;
	for(ix=0;ix<len_x;++ix){
		for(iy=0;iy<ny_upd;++iy){
			i = mat(NUM_X,NUM_Y2,ix+st_x1,iy,0);
			j = mat(nx_upd,ny_upd,ix+st_x2,iy,0);
			var_new[j] = var_old[i]/(1.0*dim);
		}
	}return 0;
}

int pad_real(double *var_new,double *var_old,int nx_upd,int ny_upd,int len_x,
int st_x1,int st_x2,int iz)
{
	int ix,iy,i,j;
	for(ix=0;ix<len_x;++ix){
		for(iy=0;iy<ny_upd;++iy){
			i = mat(NUM_X,NUM_Y2,ix+st_x1,iy,0);
			j = mat(nx_upd,ny_upd,ix+st_x2,iy,iz);
			var_new[i] = var_old[j];
		}
	}return 0;
}

int pad_cplx(fftw_complex *var_new,fftw_complex *var_old,int nx_upd,int ny_upd,
int len_x,int st_x1,int st_x2,int iz)
{
	int ix,iy,i,j;
	for(ix=0;ix<len_x;++ix){
		for(iy=0;iy<ny_upd;++iy){
			i = mat(NUM_X,NUM_Y2,ix+st_x1,iy,0);
			j = mat(nx_upd,ny_upd,ix+st_x2,iy,iz);
			var_new[i] = var_old[j];
		}
	}return 0;
}

void f_print(double *var,int n1,int n2)
{
	int ix,iy;
	for(ix=0;ix<n1;++ix){
		for(iy=0;iy<n2;++iy) printf("%.0f ",floor(var[iy+n2*ix]));
		printf("         ix = %d \n",ix);
	}
}

double promedio(double *var,int dim,int level)
{
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;++i) sum = sum + var[i+dim*level];
	return sum/((double) dim);
}

double promedioa(double *var,int dim)
{
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;++i) sum = sum + fabs(var[i]);
	return (1.0/dim)*sum;
}


int Qv_w(double q[], int m, double dz, double qv0)
{
	double zk, k0, k1, k2, k3, k4, pz, aux1, aux2;
	int iz;

	k0 = 18.04;
	k1 = 3.27;
	k2 = 0.1;
	k3 = 0.1804;
	k4 = 3.48;

	for(iz=0;iz<m;++iz){
		zk = iz*dz;
		aux1 = 1.0 + k2*zk;
		aux2 = 1.0 - k1*log(1.0*aux1);
		pz = pow(1.0*aux2,k4);
		q[iz] = (qv0/pz)*exp(k0 - k0/(aux1*aux2));
	}
	return 0;
}

int Fdqvdz(double d[],int m, double dz, double qv)
{
	double zk, k0, k1, k2, k3, k4, c1, pz;
	double n1,n2,d1,d2;
	int iz;

	k0 = 18.04;
	k1 = 3.27;
	k2 = 0.1;
	k3 = 0.1804;
	k4 = 3.48;

	for(iz=0; iz<m; ++iz){
		zk = iz*dz;
		pz = pow(1.0 - k1*log(1.0 + k2*zk),k4);
		c1 = (qv/pz)*exp(-k0*(1.0/((1.0-k1*log(1.0+k2*zk))*(1.0+k2*zk))-1.0));
		n1 = k1*k2*k4;
		d1 = (1.0 - k1*log(1.0 + k2*zk))*(1.0 + k2*zk);
		n2 = k0*k2*(1.0 - k1 - k1*log(1.0 + k2*zk));
		d2 = pow((1.0 - k1*log(1.0 + k2*zk))*(1.0 + k2*zk),2.0);
		d[iz] = k1*k2*k4/((1.0 - k1*log(1.0 + k2*zk))*(1.0 + k2*zk));
		d[iz] = c1*(n1/d1 + n2/d2);
	}
	return 0;
}

int Qv_u(double q[], int m, double dz, double qv0)
{
	double zk, k0, k1, k2, k3, k4, pz, aux1, aux2;
	int iz;

	k0 = 18.04;
	k1 = 3.27;
	k2 = 0.1;
	k3 = 0.1804;
	k4 = 3.48;

	for(iz=1;iz<m-1;++iz){
		zk = (iz + 0.5)*dz;
		aux1 = 1.0 + k2*zk;
		aux2 = 1.0 - k1*log(1.0*aux1);
		pz = pow(1.0*aux2,k4);
		q[iz] = (qv0/pz)*exp(k0 - k0/(aux1*aux2));
	}q[0] = 2.0*q[1] - q[2];
	q[m] = 2.0*q[m-1] - q[m-2];
	return 0;
}

int PVe_M(double *qt,double *Te,double *o,int d_phy,int m,double fdz,double dTe[],double dqt[],double *PVe,double *M)
{
	int i,ic,iz;
	// Definition of M
	# pragma omp parallel for private(ic,i) shared(M)
	for(iz=0;iz<m;iz++){
		for(ic=0;ic<d_phy;ic++){
			i = ic + iz*d_phy;
			M[i] = qt[i] - (dqt[iz]/dTe[iz])*Te[i];
		}
	}

	// Definition of PVe
	# pragma omp parallel for private(ic,i) shared(PVe)
	for(iz=1;iz<m;iz++){
		for(ic=0;ic<d_phy;ic++){
			i = ic + iz*d_phy;
			PVe[i] = o[i] + fdz*(Te[i]/dTe[iz] - Te[i-d_phy]/dTe[iz-1]);
		}
	}// Boundary conditions for PVe
	# pragma omp parallel for shared(PVe)
	for(ic=0;ic<d_phy;ic++){
		PVe[ic] = -PVe[ic+d_phy];
		PVe[ic+m*d_phy] = PVe[ic+(m-1)*d_phy];
	}

	return 0;
}

int ThetaE_Qt(double *M,double *PVe,double *o,int d_phy,int m,double dzf,double dTe[],double dqt[],double *Te,double *qt)
{
	int i,ic,iz;
	# pragma omp parallel for private(iz,i) shared(Te)
	for(ic=0;ic<d_phy;ic++){
		Te[ic + (m-1)*d_phy] = 0.0;
		for(iz=m-2;iz>-1;iz--){
			i = ic + iz*d_phy;
			Te[i] = dTe[iz]*(Te[i+d_phy]/dTe[iz+1] + dzf*(o[i+d_phy] - PVe[i+d_phy]));
		}
	}

	# pragma omp parallel for private(iz,i) shared(qt)
	for(ic=0;ic<d_phy;ic++){
		for(iz=m-1;iz>-1;iz--){
			i = ic + iz*d_phy;
			qt[i] = M[i] + (dqt[iz]/dTe[iz])*Te[i];
		}
	}
	return 0;
}

int F_PVe_M(fftw_complex *qt,fftw_complex *Te,fftw_complex *o,int d_upd,int m,double fdz,double dTe[],double dqt[],fftw_complex *PVe,fftw_complex *M)
{
	int i,ic,iz;
	// Definition of M
	# pragma omp parallel for private(ic,i) shared(M)
	for(iz=0;iz<m;iz++){
		for(ic=0;ic<d_upd;ic++){
			i = ic + iz*d_upd;
			M[i] = qt[i] - (dqt[iz]/dTe[iz])*Te[i];
		}
	}

	// Definition of PVe
	# pragma omp parallel for private(ic,i) shared(PVe)
	for(iz=1;iz<m;iz++){
		for(ic=0;ic<d_upd;ic++){
			i = ic + iz*d_upd;
			PVe[i] = o[i] + fdz*(Te[i]/dTe[iz] - Te[i-d_upd]/dTe[iz-1]);
		}
	}// Boundary conditions for PVe
	# pragma omp parallel for shared(PVe)
	for(ic=0;ic<d_upd;ic++){
		PVe[ic] = -PVe[ic+d_upd];
		PVe[ic+m*d_upd] = PVe[ic+(m-1)*d_upd];
	}

	return 0;
}

int F_ThetaE_Qt(fftw_complex *M,fftw_complex *PVe,fftw_complex *o,int d_upd,int m,double dzf,double dTe[],double dqt[],fftw_complex *Te,fftw_complex *qt)
{
/*
	int i,ic,iz;
	fftw_complex IM = 0.0 + 0.0*I;
	# pragma omp parallel for private(iz,i) shared(Te)
	for(ic=0;ic<d_upd;ic++){
		Te[ic + (m-1)*d_upd] = 0.0 + 0.0*IM;
		for(iz=m-2;iz>-1;iz--){
			i = ic + iz*d_upd;
			Te[i] = Te[i+d_upd] + 0.5*dzf*(dTe[iz+1] + dTe[iz])*(o[i+d_upd] - PVe[i+d_upd]);
		}
	}

	# pragma omp parallel for private(iz,i) shared(qt)
	for(ic=0;ic<d_upd;ic++){
		for(iz=m-1;iz>-1;iz--){
			i = ic + iz*d_upd;
			qt[i] = M[i] + (dqt[iz]/dTe[iz])*Te[i];
		}
	}
*/
	int i,ic,iz;
	fftw_complex IM = 0.0 + 0.0*I;
	# pragma omp parallel for private(iz,i) shared(qt)
	for(ic=0;ic<d_upd;ic++){
		qt[ic + (m-1)*d_upd] = M[ic + (m-1)*d_upd];
		for(iz=m-2;iz>-1;iz--){
			i = ic + iz*d_upd;
			qt[i] = M[i] + dqt[iz]*((qt[i+d_upd] - M[i+d_upd])/dqt[iz+1] + dzf*(o[i+d_upd] - PVe[i+d_upd]));
		}
	}

	# pragma omp parallel for private(iz,i) shared(Te)
	for(ic=0;ic<d_upd;ic++){
		for(iz=m-1;iz>-1;iz--){
			i = ic + iz*d_upd;
			Te[i] = (dTe[iz]/dqt[iz])*(qt[i] - M[i]);
		}
	}
	return 0;
}


void leerarchivo(char *nombre,double *var,int nx,int ny,int m)
{
	int ix,iy,iz;
    FILE *archivo=NULL;
    archivo=fopen(nombre,"r");
    if(archivo==NULL){
        puts("No se pudo abrir el archivo.");
        exit(0);
    }for(iz=0;iz<m;++iz){
		for(iy=0;iy<ny;++iy){
			for(ix=0;ix<nx;++ix) fscanf(archivo,"%lf",&var[iy+ny*ix+nx*ny*iz]);
		}
	}fclose(archivo);
}

double minimum(double a,double b,double c)
{
	double aux1,aux2;
	aux1=min(a,b);
	aux2=min(aux1,c);
	return aux2;
}

double min(double a,double b){
	double aux;
	aux=a;
	if(b<aux) aux=b;
	return aux;
}

double max(double a,double b){
	double aux;
	aux=a;
	if(b>aux) aux=b;
	return aux;
}

double max_qt_level(double *qt,int dim,int iz)
{
	double max;
	int i;
	max = qt[iz*dim];
	for(i=1;i<dim;i++) if(qt[i+dim*iz]>max) max = qt[i+dim*iz];
	return max;
}

double min_qt_level(double *qt,int dim,int iz)
{
	double min;
	int i;
	min = qt[iz*dim];
	for(i=1;i<dim;i++) if(qt[i+dim*iz]<min) min = qt[i+dim*iz];
	return min;
}

double max_vecs(double a[], int dim)
{
	int i;
	double max;
	max=a[0];
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i) if(a[i]>max) max=a[i];
	//printf("max = %.15e\n",vmax);
	return max;
}


double max_vecp2(double *a,int iz,int dim)
{
	int i,max_index;
	double vmax;
	vmax = a[dim*iz];
	max_index=0;
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i){
		if(vmax < a[i+dim*iz]){
			vmax = a[i+dim*iz];
			max_index = i;
		}
	}
	//printf(" the max val is find at %d\n",max_index);
	return vmax;
}

double max_vecp(double *a,int dim)
{
	int i;
	double vmax;
	vmax = a[0];
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i){
		if(vmax < a[i]) vmax = a[i];
	}
	//printf("max = %.15e\n",vmax);
	return vmax;
}

double min_vecp(double *a,int dim)
{
	int i;
	double vmin;
	vmin = a[0];
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i){
		if(vmin > a[i]) vmin = a[i];
	}
	//printf("max = %.15e\n",vmax);
	return vmin;
}


double abs_vecp(double *a,int dim)
{
	int i;
	double vmax;
	vmax = fabs(a[0]);
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i){
		if(vmax < fabs(a[i])) vmax = fabs(a[i]);
	}
	//printf("max = %.15e\n",vmax);
	return vmax;
}

double max_complexvar(fftw_complex *var,int dim)
{
	int i,max_index;
	double vmax;
	vmax = cabs(var[0]);
	max_index = 0;
	//printf("max = %.15e\n",vmax);
	for(i=1;i<dim;++i){
		if(vmax < cabs(var[i])){
			vmax = cabs(var[i]);
			max_index = i;
		}
	}
	return vmax;
}

double max_varhat_level(fftw_complex *var,int dim,int iz)
{
	double max;
	int i;
	max = cabs(var[iz*dim]);
	for(i=1;i<dim;i++) if(cabs(var[i+dim*iz])>max) max = cabs(var[i+dim*iz]);
	return max;
}

void Fluct(double *u,double *v,double *w,double *qt,double *Te,int dim,int m,double Us,double Ths,double qs)
{
	double d_u = dim*(m+1),d_w = dim*m;
	printf("\n\n    Value of all Fluctiating-variables\n");
	printf("         Max u =  %.15e\n",Us*max_vecp(u,d_u));
	printf("         Min u = %.15e\n",Us*min_vecp(u,d_u));
	printf("\n");
	printf("         Max v =  %.15e\n",Us*max_vecp(v,d_u));
	printf("         Min v = %.15e\n",Us*min_vecp(v,d_u));
	printf("\n");
	printf("         Max w =  %.15e\n",Us*max_vecp(w,d_w));
	printf("         Min w = %.15e\n",Us*min_vecp(w,d_w));
	printf("\n");
	printf("        Max qt =  %.15e\n",qs*max_vecp(qt,d_w));
	printf("        Min qt = %.15e\n",qs*min_vecp(qt,d_w));
	printf("\n");
	printf("        Max Te =  %.15e\n",Ths*max_vecp(Te,d_w));
	printf("        Min Te = %.15e\n",Ths*min_vecp(Te,d_w));
	printf("\n");
}

void Thermo(double *qv,double *qr,double *Th,double *qt,double *Te,int dim,double Ths,double qs)
{
	printf("\n\n    Value of all Thermo-variables\n");
	printf("        Max qv = %.15e\n",qs*max_vecp(qv,dim));
	printf("        Min qv = %.15e\n\n",qs*min_vecp(qv,dim));
	printf("        Max qr = %.15e\n",qs*max_vecp(qr,dim));
	printf("        Min qr = %.15e\n\n",qs*min_vecp(qr,dim));
	printf("        Max qt = %.15e\n",qs*max_vecp(qt,dim));
	printf("        Min qt = %.15e\n\n",qs*min_vecp(qt,dim));
	printf("        Max Th = %.15e\n",Ths*max_vecp(Th,dim));
	printf("        Min Th = %.15e\n\n",Ths*min_vecp(Th,dim));
	printf("        Max Te = %.15e\n",Ths*max_vecp(Te,dim));
	printf("        Min Te = %.15e\n\n",Ths*min_vecp(Te,dim));
}

void ABS_Real(double *u,double *v,double *p,double *o,double *w,double *qt,double *Te,int d_phy,int m,double Us,double Ls,double Ts,double Ths,double qs)
{
	printf("\n\n    Abs value of all variables\n");
	printf("        Abs p = %.15e\n",Us*Ls*abs_vecp(p,d_phy*(m+1)));
	printf("        Abs u = %.15e\n",Us*abs_vecp(u,d_phy*(m+1)));
	printf("        Abs v = %.15e\n",Us*abs_vecp(v,d_phy*(m+1)));
	printf("        Abs o = %.15e\n",Ts*abs_vecp(o,d_phy*(m+1)));
	printf("        Abs w = %.15e\n",Us*abs_vecp(w,d_phy*m));
	printf("       Abs qt = %.15e\n",qs*abs_vecp(qt,d_phy*m));
	printf("       Abs Te = %.15e\n",Ths*abs_vecp(Te,d_phy*m));
}

void MAX_Real(double *u,double *v,double *p,double *o,double *w,double *qt,double *Te,int d_phy,int m,double Us,double Ls,double Ts,double Ths,double qs)
{
	printf("\n\n    Max value of all variables\n");
	printf("        Max p = %.15e\n",Us*Ls*max_vecp(p,d_phy*(m+1)));
	printf("        Max u = %.15e\n",Us*max_vecp(u,d_phy*(m+1)));
	printf("        Max v = %.15e\n",Us*max_vecp(v,d_phy*(m+1)));
	printf("        Max o = %.15e\n",Ts*max_vecp(o,d_phy*(m+1)));
	printf("        Max w = %.15e\n",Us*max_vecp(w,d_phy*m));
	printf("       Max qt = %.15e\n",qs*max_vecp(qt,d_phy*m));
	printf("       Max Te = %.15e\n",Ths*max_vecp(Te,d_phy*m));
}

void MIN_Real(double *u,double *v,double *p,double *o,double *w,double *qt,double *Te,int d_phy,int m,double Us,double Ls,double Ts,double Ths,double qs)
{
	printf("\n\n    Min value of all variables\n");
	printf("        Min p = %.15e\n",Us*Ls*min_vecp(p,d_phy*(m+1)));
	printf("        Min u = %.15e\n",Us*min_vecp(u,d_phy*(m+1)));
	printf("        Min v = %.15e\n",Us*min_vecp(v,d_phy*(m+1)));
	printf("        Min o = %.15e\n",Ts*min_vecp(o,d_phy*(m+1)));
	printf("        Min w = %.15e\n",Us*min_vecp(w,d_phy*m));
	printf("       Min qt = %.15e\n",qs*min_vecp(qt,d_phy*m));
	printf("       Min Te = %.15e\n",Ths*min_vecp(Te,d_phy*m));
}

void MAX_Cplx(fftw_complex *u,fftw_complex *v,fftw_complex *p,fftw_complex *o,fftw_complex *w,fftw_complex *qt,fftw_complex *Te,int d_upd,int m)
{
	printf("\n\n    Max value of all hat variables\n");
	printf("        Max p = %.15e\n",max_complexvar(p,d_upd*(m+1)));
	printf("        Max u = %.15e\n",max_complexvar(u,d_upd*(m+1)));
	printf("        Max v = %.15e\n",max_complexvar(v,d_upd*(m+1)));
	printf("        Max o = %.15e\n",max_complexvar(o,d_upd*(m+1)));
	printf("        Max w = %.15e\n",max_complexvar(w,d_upd*m));
	printf("       Max qt = %.15e\n",max_complexvar(qt,d_upd*m));
	printf("       Max Te = %.15e\n",max_complexvar(Te,d_upd*m));
}

void MAX_Flux(fftw_complex *fo,fftw_complex *fw,fftw_complex *fqt,fftw_complex *fTe,int d_upd,int m)
{
	printf("\n\n    Max value of the 4 hat variables\n");
	printf("        Max fo = %.15e\n",max_complexvar(fo,d_upd*(m+1)));
	printf("        Max fw = %.15e\n",max_complexvar(fw,d_upd*m));
	printf("       Max fqt = %.15e\n",max_complexvar(fqt,d_upd*m));
	printf("       Max fTe = %.15e\n",max_complexvar(fTe,d_upd*m));
}

void AbrirArchivos(int fileindex)
{
    FILE *salida=NULL;
	char eps_name[80];

	// bar
	sprintf(eps_name,"Run_%d/Data1D/Qt_bar.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/Te_bar.dat",fileindex);
	salida=fopen(eps_name,"w");

	// background
	sprintf(eps_name,"Run_%d/Data1D/Qt_bg.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/Te_bg.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/Qvs.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/dqt_bg.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/RC.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data1D/Mr.dat",fileindex);
	salida=fopen(eps_name,"w");

	// Z max
	sprintf(eps_name,"Run_%d/Data2D/QrZMax.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/QvZMax.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/ThZMax.dat",fileindex);
	salida=fopen(eps_name,"w");

	// Z avg
	sprintf(eps_name,"Run_%d/Data2D/QrZavg.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/QvZavg.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/ThZavg.dat",fileindex);
	salida=fopen(eps_name,"w");

	// y half
	sprintf(eps_name,"Run_%d/Data2D/QrYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/QvYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/ThYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/UYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/VYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
	sprintf(eps_name,"Run_%d/Data2D/WYHalf.dat",fileindex);
	salida=fopen(eps_name,"w");
}

void Impresion1D(double *qt,double *Te,int m,int fileindex,double Ths,double qs)
{
	Imprime1D(qt,1,m,fileindex,qs);
	Imprime1D(Te,2,m,fileindex,Ths);
}

void Imprime1D(double *Y,int option,int m,int fileindex,double sc)
{
	FILE *salida=NULL;
	char eps_name[80];
	int iz;
	switch(option){
		case 1:{
			sprintf(eps_name,"Run_%d/Data1D/Qt_bar.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 2:{
			sprintf(eps_name,"Run_%d/Data1D/Te_bar.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 3:{
			sprintf(eps_name,"Run_%d/Data1D/Qt_bg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 4:{
			sprintf(eps_name,"Run_%d/Data1D/Te_bg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 5:{
			sprintf(eps_name,"Run_%d/Data1D/Qvs.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 6:{
			sprintf(eps_name,"Run_%d/Data1D/dqt_bg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 7:{
			sprintf(eps_name,"Run_%d/Data1D/RC.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 8:{
			sprintf(eps_name,"Run_%d/Data1D/Mr.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		}if(salida==NULL){
		printf("No se pudo abrir el archivo %d",option);
		printf("--> %s\n",eps_name);
		exit(0);
	}for(iz=0;iz<m;++iz) fprintf(salida," %.15e\n",sc*Y[iz]);
	fclose(salida);
	salida=NULL;
}

void Impresion2D(double *u,double *v,double *w,double *qv,double *qr,double *Th,int nx,int ny,int m,int fileindex,double Us,double Ths,double qs)
{
	// Z max
	Imprime2D(qr,1,m,ny,nx,fileindex,qs);
	Imprime2D(qv,2,m,ny,nx,fileindex,qs);
	Imprime2D(Th,3,m,ny,nx,fileindex,Ths);

	// Z avg
	Imprime2D(qr,11,m,ny,nx,fileindex,qs);
	Imprime2D(qv,12,m,ny,nx,fileindex,qs);
	Imprime2D(Th,13,m,ny,nx,fileindex,Ths);

	// ny half
	Imprime2D(qr,21,m,ny,nx,fileindex,qs);
	Imprime2D(qv,22,m,ny,nx,fileindex,qs);
	Imprime2D(Th,23,m,ny,nx,fileindex,Ths);
	Imprime2D(u,24,m+1,ny,nx,fileindex,Us);
	Imprime2D(v,25,m+1,ny,nx,fileindex,Us);
	Imprime2D(w,26,m,ny,nx,fileindex,Us);
}


void Imprime2D(double *Y,int option,int m,int ny,int nx,int fileindex,double sc)
{
    FILE *salida=NULL;
	char eps_name[80];
    int ix,iy,iz;
	double aux;

	switch(option){
		// Z max
		case 1:{
			sprintf(eps_name,"Run_%d/Data2D/QrZMax.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 2:{
			sprintf(eps_name,"Run_%d/Data2D/QvZMax.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 3:{
			sprintf(eps_name,"Run_%d/Data2D/ThZMax.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		// Z avg
		case 11:{
			sprintf(eps_name,"Run_%d/Data2D/QrZavg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 12:{
			sprintf(eps_name,"Run_%d/Data2D/QvZavg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 13:{
			sprintf(eps_name,"Run_%d/Data2D/ThZavg.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		// ny half
		case 21:{
			sprintf(eps_name,"Run_%d/Data2D/QrYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
            }break;
		case 22:{
			sprintf(eps_name,"Run_%d/Data2D/QvYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 23:{
			sprintf(eps_name,"Run_%d/Data2D/ThYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 24:{
			sprintf(eps_name,"Run_%d/Data2D/UYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 25:{
			sprintf(eps_name,"Run_%d/Data2D/VYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
		case 26:{
			sprintf(eps_name,"Run_%d/Data2D/WYHalf.dat",fileindex);
			salida=fopen(eps_name,"a");
			}break;
	}if(salida==NULL){
		puts("No se pudo abrir el archivo.");
		exit(0);
	}switch(option){
		case 1:
		case 2:
		case 3:{
			double MV[m];
			aux=0.0;
			for(ix=0;ix<nx;++ix){
				for(iy=0;iy<ny;++iy){
					for(iz=0;iz<m;++iz)	MV[iz] = Y[mat(nx,ny,ix,iy,iz)];
					aux = max_vecs(MV,m);
					fprintf(salida,"%.15e\n",sc*aux);
				}
			}
		}break;
		case 11:
		case 12:
		case 13:{
			for(ix=0;ix<nx;++ix){
				for(iy=0;iy<ny;++iy){
					aux=0.0;
					for(iz=0;iz<m;++iz) aux=aux+Y[mat(nx,ny,ix,iy,iz)];
					fprintf(salida,"%.15e\n",sc*aux/m);
				}
			}
		}break;
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		case 26:{
			iy=ny/2;
			for(iz=0;iz<m;++iz){
				for(ix=0;ix<nx;++ix) fprintf(salida,"%.15e\n",sc*Y[mat(nx,ny,ix,iy-1,iz)]);
			}
		}break;
	}fclose(salida);
	salida=NULL;
}

void Impresion3D(double *u,double *v,double *w,double *qv,double *qr,double *Th,int nx,int ny,int m,int file,int fileindex,double Us,double Ths,double qs)
{
	Imprime3D(qr,1,m,ny,nx,file,fileindex,qs);
	Imprime3D(qv,2,m,ny,nx,file,fileindex,qs);
	Imprime3D(Th,3,m,ny,nx,file,fileindex,Ths);
	Imprime3D(u,4,m+1,ny,nx,file,fileindex,Us);
	Imprime3D(v,5,m+1,ny,nx,file,fileindex,Us);
	Imprime3D(w,6,m,ny,nx,file,fileindex,Us);
}

void Imprime3D(double *X,int option,int m,int ny,int nx,int file,int fileindex,double sc) //--- For printing, m=m, n=ny, r=nx
{
    FILE *salida=NULL;
	char eps_name[80];
    int ix,iy,iz;

    switch(option){
		case 1:{
			sprintf(eps_name,"Run_%d/Data3D/Qr_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
		case 2:{
			sprintf(eps_name,"Run_%d/Data3D/Qv_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
		case 3:{
			sprintf(eps_name,"Run_%d/Data3D/Th_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
		case 4:{
			sprintf(eps_name,"Run_%d/Data3D/U_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
		case 5:{
			sprintf(eps_name,"Run_%d/Data3D/V_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
		case 6:{
			sprintf(eps_name,"Run_%d/Data3D/W_%d.dat",file,fileindex);
			salida=fopen(eps_name,"w");
			}break;
	}if(salida==NULL){
		puts("No se pudo abrir el archivo.");
		exit(0);
	}for(iz=0;iz<m;++iz){
		for(iy=0;iy<ny;++iy){
			for(ix=0;ix<nx;++ix){
				if(X[mat(nx,ny,ix,iy,iz)]<0.0) fprintf(salida," %.15e\n",sc*X[mat(nx,ny,ix,iy,iz)]);
				else fprintf(salida,"  %.15e\n",sc*X[mat(nx,ny,ix,iy,iz)]);
			}
		}
	}fclose(salida);
	salida=NULL;
}

int upad_cplx(fftw_complex *var_new,fftw_complex *var_old,int nx_upd,int ny_upd,
int len_x,int st_x1,int st_x2,int dim,int iz)
{
	int ix,iy,i,j;
	for(ix=0;ix<len_x;++ix){
		for(iy=0;iy<ny_upd;++iy){
			i = mat(NUM_X,NUM_Y2,ix+st_x1,iy,0);
			j = mat(nx_upd,ny_upd,ix+st_x2,iy,iz);
			var_new[j] = var_old[i]/(1.0*dim);
		}
	}return 0;
}

