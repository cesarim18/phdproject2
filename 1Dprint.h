// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 

void ImprimirPar(int file,int elec,double value)
{
	char eps_name[80];
	FILE *salida=NULL;
	switch(elec){
		case 0:{
			sprintf(eps_name,"Run_%d/Parameters/dt.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 1:{
			sprintf(eps_name,"Run_%d/Parameters/phat.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 2:{
			sprintf(eps_name,"Run_%d/Parameters/div.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 3:{
			sprintf(eps_name,"Run_%d/Parameters/time.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 4:{
			sprintf(eps_name,"Run_%d/Parameters/relation.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 16:{
			sprintf(eps_name,"Run_%d/Parameters/norml1.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 17:{
			sprintf(eps_name,"Run_%d/Parameters/norml2.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 18:{
			sprintf(eps_name,"Run_%d/Parameters/timewhile.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 19:{
			sprintf(eps_name,"Run_%d/Parameters/max_qv_ini.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 20:{
			sprintf(eps_name,"Run_%d/Parameters/max_qr.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 21:{
			sprintf(eps_name,"Run_%d/Parameters/qvs.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 22:{
			sprintf(eps_name,"Run_%d/Parameters/max_qv.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 23:{
			sprintf(eps_name,"Run_%d/Parameters/max_qt.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 24:{
			sprintf(eps_name,"Run_%d/Parameters/min_w.dat",file);
			salida=fopen(eps_name,"a");
			}break;
		case 25:{
			sprintf(eps_name,"Run_%d/Parameters/max_w.dat",file);
			salida=fopen(eps_name,"a");
			}break;
	}if(salida==NULL){
		puts("No se pudo abrir el archivo.");
		exit(0);
	}fprintf(salida,"%.15e\n",value);
	fclose(salida);
    salida=NULL;
}
