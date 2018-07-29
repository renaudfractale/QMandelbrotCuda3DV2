#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10

#include <iostream> // prompt Output
#include <fstream> //File Output
#include <math.h> //lib mayh
#include <stdio.h> // lib stantard
#include <cuda_fp16.h> // lib CUDA
#include <windows.h>
#include <fstream>
#include <string>


#define Dim_isFix 0;
#define Dim_end 20.0f;
#define Dim_start -20.0f;
#define Dim_NbStep 4;
#define Dim_step 10.0f;

#define ITER_MAX 255;
#define ITER_MIN 1;
#define ITER_isFix 0;

#define DEV 1;
#define FILTER 0;
#define POWER 2.0f;
#define ISSHOW 1;
#define RMAX 4.0f;

#define NBPOINTS 64;

#define MODEAUTOMANU false;
//Strcture state

 struct struct_Stat_float_T {
	float Xmin = Dim_end
	float Xmax = Dim_start
	float Wmin = Dim_end
	float Wmax = Dim_start
	float Ymin = Dim_end
	float Ymax = Dim_start
	float Zmin = Dim_end
	float Zmax = Dim_start
	float Wstep = 0.0f;
	float Xstep = 0.0f;
	float Ystep = 0.0f;
	float Zstep = 0.0f;
	unsigned long NbPoint = 0;
 };


// struct sur la gestion des dimensions
 struct struct_P_float_T {
	 int isFix = Dim_isFix
	 float start = Dim_start
	 float end = Dim_end
	 int NbStep = Dim_NbStep
	 float step = Dim_step
		int coef = 1;
} ;


struct 	 struct_Iter_T  {
	int max = ITER_MAX
	int min = ITER_MIN
	int isFix = ITER_isFix
} ;

typedef struct 	struct_FileName_T {
	char root[100];
	char txt[110];
	char csv[110];
	char stat[110];
	char histo[110];
	char stl[110];
} struct_FileName_T;

// struct sur la gestion param�tres d'entr�
struct struct_P_All_T {
	struct_P_float_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	struct_Iter_T Iter;
	int dev = DEV
		int filter = FILTER
		float power = POWER
		int isShow = ISSHOW
		struct_FileName_T nameFile;
	float rMax = RMAX
		int NbPointByStep = NBPOINTS
		bool modeAM = MODEAUTOMANU
} ;

struct 	struct_P_Simulation_T {
	//Quaternions
	struct_P_float_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	struct_Iter_T Iter;
	float rMax;
	//Parametrer variable systematique
	float power;
	int max;
} ;

struct struct_Q_T {
	float x;
	float y;
	float z;
	float w;
} ;

struct struct_QD_T {
	double x;
	double y;
	double z;
	double w;
};


__device__  void CreateQ_By_float(struct_QD_T *out, double x, double y, double z, double w)
{
	out->x = x;
	out->y = y;
	out->z = z;
	out->w = w;
}

//__device__  float  Get_QNorm(struct_Q_T *Q)
//{
//	float R2 = Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w;
//	float N = norm4df(Q->x, Q->y, Q->z, Q->w);
//	float R = sqrtf(Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w);
//	return R;
//}

__device__ void Get_QPow(struct_QD_T *Q, float power)
{
	
	double A = norm4d(Q->x, Q->y, Q->z, Q->w);
	//printf("A=%f 1er\n", A);
	double theta = 0.0;
	double B = 0.0;
	double R = 0.0;
	double Temp = 0.0;
	double Temp2 = 0.0;
	if (power > 0.0f && A>0.00001)
	{
		//double coef = 1.0;
		//if (A<1.0)
		//{
		//	//printf("1.0f > %f(A) *******\n", A);
		//	coef = 1 / A;
		//	Q->x = __dmul_rn(Q->x, coef);
		//	Q->y = __dmul_rn(Q->y, coef);
		//	Q->z = __dmul_rn(Q->z, coef);
		//	Q->w = __dmul_rn(Q->w, coef);
		//	if (Q->x == NAN || Q->y == NAN || Q->z == NAN || Q->w == NAN)
		//		return;
		//	//printf("coef=%f\n", coef);
		//	A = 1.0;
		//	//printf("A=%f=%f 2eme\n", A, norm4df(Q->x, Q->y, Q->z, Q->w));

		//}
		
		
		/*printf("Q->w / A = %f/%f \n",  Q->w, A);
		printf("Q->w / A = %f \n", Q->w / A, Q->w,A);*/

		//printf("%f +++++++++\n", A);
		Temp = __ddiv_rd(Q->w, A);
		if (abs(Temp) <= 1.0)
		{
			Temp2 = acos(Temp);
			theta = __dmul_rn(Temp2, power);
		}
		else
		{
			printf("Error __ddiv_rd(Q->w, A)=__ddiv_rd(%f, %f)=%f", Q->w, A, Temp);
		}
		
		if (abs(theta) <= 0.001)
			theta = 0.0;

		//B = __dmul_rn(A,A) - __dmul_rn(Q->w,Q->w);
		B = norm3d(Q->x, Q->y, Q->z);
		//caluer debug

	/*	R = pow(A / coef, (double)(power));
		if (R == NAN )
			return;*/
		double C = log(A); //__ddiv_rn(A , coef)
		C = __dmul_rn( C,power);
		if (C > -2.661474648)
			R = exp(C);
		else
		{
			Q->w = 0.0;
			Q->x = 0.0;
			Q->y = 0.0;
			Q->z = 0.0;
			return;
			//printf("%f +++++++++\n", R);
		}
			
		//R = powf(A / coef, pow);
		Temp = sin(theta);
		if (abs(Temp) > 0.001)
		{
			Temp2 = __dmul_rn(R, Temp);
			if (abs(B) > 0.01 || abs(Q->x) > 0.01)
				Q->x = __dmul_rn(Temp2, (Q->x / B));
			else
				Q->x = 0;

			Temp2 = __dmul_rn(R, Temp);
			if (abs(B) > 0.01 || abs(Q->y) > 0.01)
				Q->y = __dmul_rn(Temp2, (Q->y / B));
			else
				Q->y = 0;

			Temp2 = __dmul_rn(R, Temp);
			if (abs(B) > 0.01 || abs(Q->z) > 0.01)
				Q->z = __dmul_rn(Temp2, (Q->z / B));
			else
				Q->z = 0;
		}
		else
		{
			Q->x = 0.0;
			Q->y = 0.0;
			Q->z = 0.0;

		}


		Temp = cos(theta);
		if (abs(Temp) > 0.001)
			Q->w = __dmul_rn(R, Temp);
		else
			Q->w = 0.0;


	}
	else
	{
		//printf("%f --------\n", A);
		Q->w = 0.0;
		Q->x = 0.0;
		Q->y = 0.0;
		Q->z = 0.0;
	}
}


__host__  void CreateQ_By_float_H(struct_Q_T *out, float x, float y, float z, float w)
{
	out->x = x;
	out->y = y;
	out->z = z;
	out->w = w;
}

__host__  float  Get_QNorm_H(struct_Q_T *Q)
{
	return sqrtf(Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w);
}

__host__ void Get_QPow_H(struct_Q_T *Q, float pow)
{
	float A = Get_QNorm_H(Q);
	float theta = 0.0f;
	float B = 0.0f;
	float R = 0.0f;
	if (pow > 0.0f && A>0.000001f)
	{
		float coef = 1.0f;
		if (A<1.0f)
		{
			//printf("%f *******\n", A);
			coef = 1 / A;
			Q->x *= coef;
			Q->y *= coef;
			Q->z *= coef;
			Q->w *= coef;

		}
		A = Get_QNorm_H(Q);
		//printf("%f +++++++++\n", A);
		theta = acosf(Q->w / A)*pow;
		B = sqrt(A*A - Q->w*Q->w);
		R = exp2f(logf(A / coef)* pow);
		Q->x = R*sinf(theta)*(Q->x / B);
		Q->y = R*sinf(theta)*(Q->y / B);
		Q->z = R*sinf(theta)*(Q->z / B);
		Q->w = R*cosf(theta);

	}
	else
	{
		//printf("%f --------\n", A);
		Q->w = 0.0f;
		Q->x = 0.0f;
		Q->y = 0.0f;
		Q->z = 0.0f;

	}
}

__host__ int  GetQIter_H(struct_P_Simulation_T *P_Simulation_DEVICE, int  *x_filter, int  *y_filter, int *z_filter, int *w_filter)
{
	//int Tempindex = 0;
	struct_Q_T Q_Current;
	float w, x, y, z;
	int iter_computed;
	//X
	x = ((float)*x_filter)*P_Simulation_DEVICE->X.step + P_Simulation_DEVICE->X.start;
	//Y
	y = ((float)*y_filter)*P_Simulation_DEVICE->Y.step + P_Simulation_DEVICE->Y.start;
	//Z
	z = ((float)*z_filter)*P_Simulation_DEVICE->Z.step + P_Simulation_DEVICE->Z.start;
	//W
	w = ((float)*w_filter)*P_Simulation_DEVICE->W.step + P_Simulation_DEVICE->W.start;



	CreateQ_By_float_H(&Q_Current, x, y, z, w);

	for (iter_computed = 0; iter_computed <= P_Simulation_DEVICE->Iter.max; iter_computed++)
	{
		Get_QPow_H(&Q_Current, P_Simulation_DEVICE->power);
		Q_Current.x += x;
		Q_Current.y += y;
		Q_Current.z += z;
		Q_Current.w += w;

		if (Get_QNorm_H(&Q_Current) > P_Simulation_DEVICE->rMax)
		{
			if (iter_computed > 0)
				iter_computed--;
			return iter_computed;
		}
	}
	if (iter_computed > 0)
		iter_computed--;
	return iter_computed;
}
__host__ bool  FilterQ_H(int *Filter, int *Nx, int *Ny, int *Nz, int *Nw, int iter,	struct_P_Simulation_T *P_Simulation)
{
	if (*Filter == 0)
		return true;
	int iter_computed = 0;

	int pasx = 1;
	if (P_Simulation->X.NbStep == 1)
		pasx = 0;

	int pasy = 1;
	if (P_Simulation->Y.NbStep == 1)
		pasy = 0;

	int pasz = 1;
	if (P_Simulation->Z.NbStep == 1)
		pasz = 0;

	int pasw = 1;
	if (P_Simulation->W.NbStep == 1)
		pasw = 0;


	for (int x_filter = *Nx - pasx; x_filter <= *Nx + pasx; x_filter++)
	{
		for (int y_filter = *Ny - pasy; y_filter <= *Ny + pasy; y_filter++)
		{
			for (int z_filter = *Nz - pasz; z_filter <= *Nz + pasz; z_filter++)
			{
				for (int w_filter = *Nw - pasw; w_filter <= *Nw + pasw; w_filter++)
				{
					iter_computed = GetQIter_H(P_Simulation, &x_filter, &y_filter, &z_filter, &w_filter);
					if (*Filter == 1)
					{
						if (iter_computed != iter)
							return true;
					}
					else //filter==2
					{
						if (iter_computed == 0)
							return true;
					}
				}
			}
		}
	}
	return false;
}


int menu(int argc, char *argv[], struct_P_All_T *Config)
{
	//Config

	strcpy(Config->nameFile.root, "O");
	if (argc == 1)
	{
		Config->W.end = -0.3375f;
		Config->W.start = -0.3375f;
		Config->W.step = 1.0f;
		Config->W.isFix = 1;
		Config->W.NbStep = 0;

		float borne = -0.75f;
		Config->X.isFix = 2;
		Config->X.start = borne;
		Config->X.end = 0.0f;
		Config->X.NbStep = 1;

		Config->Y.isFix = 2;
		Config->Y.start = borne;
		Config->Y.end = 0.0f;
		Config->Y.NbStep = 1;

		Config->Z.isFix = 2;
		Config->Z.start = borne;
		Config->Z.end = 0.0f;
		Config->Z.NbStep = 1;

		/*Config->Z.end = -0.0;
		Config->Z.start = -0.0;
		Config->Z.step = 1.0;
		Config->Z.isFix = 1;
		Config->Z.NbStep = 1;*/



		Config->filter = 2;


		Config->modeAM = false;
		//Config->Iter.isFix = 1;
		//Config->Iter.max = 2;
		Config->Iter.min = 10;
		//Config->rMax = 2;
	}




	//Arg Help
	char Str_H[] = "-h";
	char Str_Help[] = "--help";

	//Arg X
	char Str_xFix[] = "-x";
	char Str_xMax[] = "-xmax";
	char Str_xMin[] = "-xmin";
	char Str_xNbStep[] = "-xNbStep";
	//Arg W
	char Str_wFix[] = "-w";
	char Str_wMax[] = "-wmax";
	char Str_wMin[] = "-wmin";
	char Str_wNbStep[] = "-wNbStep";
	//Arg Y
	char Str_yFix[] = "-y";
	char Str_yMax[] = "-ymax";
	char Str_yMin[] = "-ymin";
	char Str_yNbStep[] = "-yNbStep";
	//Arg Z
	char Str_zFix[] = "-z";
	char Str_zMax[] = "-zmax";
	char Str_zMin[] = "-zmin";
	char Str_zNbStep[] = "-zNbStep";
	//Arg Iter
	char Str_IterFix[] = "-iter";
	char Str_IterMax[] = "-iterMax";
	char Str_IterMin[] = "-iterMin";
	//Arg Dev
	char Str_dev[] = "-device";
	//Arg Filter
	char Str_filter[] = "-filter";
	//Arg Power
	char Str_power[] = "-power";
	//Arg IsShow
	char Str_IsShow[] = "-isShow";
	//Arg Output File
	char Str_Out[] = "-o";
	//Arg Rmax
	char Str_rMax[] = "-rMax";

	//Arg Rmax
	char Str_mode[] = "-mode";
	char Str_modeAuto[] = "auto";
	char Str_modeManu[] = "manu";
	//Si il y a des Arguments
	if (argc > 1)
	{
		// dédection la commande -help pu -h
		for (int i = 1; i < argc; i++)
		{
			if (strcmp(argv[i], Str_Help) == 0 || strcmp(argv[i], Str_H) == 0)
			{
				//Affiche Help
				FILE *fileman;
				char line[1000];
				fileman = fopen("man", "r");
				if (fileman != NULL)
				{
					while (std::fgets(line, 1000, fileman))
					{
						std::cout << line;
					}
					fclose(fileman);
					return 0; //Fin du programme
				}
				else
				{
					std::cout << "Error  : man no found" << "\n";
					return -1;
				}
			}
		}
		//Verification : si nb d'arguments est paire --> erreur
		if (argc % 2 == 0)
		{
			std::cout << "Error 00 : Argument impaire" << "\n";
			return -1;
		}
		else //: si nb d'arguments est impaire --> fonctionement normale
		{
			for (int i = 1; i < argc; i += 2)
			{
				std::cout << "Analyse du couple d'arguments :  " << argv[i] << " " << argv[i + 1] << "\n";
				if (strcmp(argv[i], Str_xFix) == 0)
				{
					if (Config->X.isFix == 0) //Si premier config
					{
						Config->X.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xFix << ": value is not float " << "\n";
							return -1;
						}
						Config->X.start = value;
						Config->X.end = value;
						Config->X.NbStep = 0;
						Config->X.step = 0.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_xFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wFix) == 0)
				{
					if (Config->W.isFix == 0) //Si premier config
					{
						Config->W.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wFix << ": value is not float " << "\n";
							return -1;
						}
						Config->W.start = value;
						Config->W.end = value;
						Config->W.NbStep = 0;
						Config->W.step = 0.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_wFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yFix) == 0)
				{
					if (Config->Y.isFix == 0) //Si premier config
					{
						Config->Y.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yFix << ": value is not float " << "\n";
							return -1;
						}
						Config->Y.start = value;
						Config->Y.end = value;
						Config->Y.NbStep = 0;
						Config->Y.step = 0.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_yFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zFix) == 0)
				{
					if (Config->Z.isFix == 0) //Si premier config
					{
						Config->Z.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zFix << ": value is not float " << "\n";
							return -1;
						}
						Config->Z.start = value;
						Config->Z.end = value;
						Config->Z.NbStep = 0;
						Config->Z.step = 0.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_zFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xMax) == 0)
				{
					if (Config->X.isFix == 0 || Config->X.isFix == 2)
					{
						Config->X.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xMax << ": value is not float " << "\n";
							return -1;
						}
						Config->X.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wMax) == 0)
				{
					if (Config->W.isFix == 0 || Config->W.isFix == 2)
					{
						Config->W.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wMax << ": value is not float " << "\n";
							return -1;
						}
						Config->W.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yMax) == 0)
				{
					if (Config->Y.isFix == 0 || Config->Y.isFix == 2)
					{
						Config->Y.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yMax << ": value is not float " << "\n";
							return -1;
						}
						Config->Y.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zMax) == 0)
				{
					if (Config->Z.isFix == 0 || Config->Z.isFix == 2)
					{
						Config->Z.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zMax << ": value is not float " << "\n";
							return -1;
						}
						Config->Z.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_zMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xMin) == 0)
				{
					if (Config->X.isFix == 0 || Config->X.isFix == 2)
					{
						Config->X.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xMin << ": value is not float " << "\n";
							return -1;
						}
						Config->X.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wMin) == 0)
				{
					if (Config->W.isFix == 0 || Config->W.isFix == 2)
					{
						Config->W.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wMin << ": value is not float " << "\n";
							return -1;
						}
						Config->W.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yMin) == 0)
				{
					if (Config->Y.isFix == 0 || Config->Y.isFix == 2)
					{
						Config->Y.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yMin << ": value is not float " << "\n";
							return -1;
						}
						Config->Y.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zMin) == 0)
				{
					if (Config->Z.isFix == 0 || Config->Z.isFix == 2)
					{
						Config->Z.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zMin << ": value is not float " << "\n";
							return -1;
						}
						Config->Z.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_zMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xNbStep) == 0)
				{
					if (Config->X.isFix == 0 || Config->X.isFix == 2)
					{
						Config->X.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_xNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config->X.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wNbStep) == 0)
				{
					if (Config->W.isFix == 0 || Config->W.isFix == 2)
					{
						Config->W.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_wNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config->W.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yNbStep) == 0)
				{
					if (Config->Y.isFix == 0 || Config->Y.isFix == 2)
					{
						Config->Y.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_yNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config->Y.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zNbStep) == 0)
				{
					if (Config->Z.isFix == 0 || Config->Z.isFix == 2)
					{
						Config->Z.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_yNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config->Z.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterFix) == 0)
				{
					if (Config->Iter.isFix == 0) //Si premier config
					{
						Config->Iter.isFix = 1;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_IterFix << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config->Iter.max = value;
						Config->Iter.min = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterMax) == 0)
				{
					if (Config->Iter.isFix == 0 || Config->Iter.isFix == 2) //Si premier config
					{
						Config->Iter.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_IterMax << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config->Iter.max = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterMin) == 0)
				{
					if (Config->Iter.isFix == 0 || Config->Iter.isFix == 2) //Si premier config
					{
						Config->Iter.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_IterMin << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config->Iter.min = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_dev) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_dev << ": value is not int" << "\n";
						return -1;
					}
					Config->dev = value;
				}
				else if (strcmp(argv[i], Str_filter) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_filter << ": value is not int" << "\n";
						return -1;
					}
					Config->filter = value;
				}
				else if (strcmp(argv[i], Str_power) == 0)
				{
					float value = (float)atof(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_filter << ": value is not float" << "\n";
						return -1;
					}
					Config->power = value;
				}
				else if (strcmp(argv[i], Str_IsShow) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_IsShow << ": value is not int" << "\n";
						return -1;
					}
					Config->isShow = value;
				}
				else if (strcmp(argv[i], Str_rMax) == 0)
				{
					float value = (float)atof(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_rMax << ": value is not float" << "\n";
						return -1;
					}
					Config->rMax = value;
				}
				else if (strcmp(argv[i], Str_Out) == 0)
				{
					if (strlen(argv[i + 1])<100)
						strcpy(Config->nameFile.root, argv[i + 1]);
					else
					{
						std::cout << "Error 04 strlen fileOutput must be inf to 100 signe \n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_mode) == 0)
				{
					if (strcmp(argv[i + 1], Str_modeAuto) == 0)
						Config->modeAM = true;
					else if (strcmp(argv[i + 1], Str_modeManu) == 0)
						Config->modeAM = false;
					else
					{
						std::cout << "Error 08 mode no found\n";
						return -1;
					}
				}
				else
				{
					std::cout << "Warning 05  Arg not know : " << argv[i] << " " << argv[i + 1] << "\n";
				}
			}
		}
	}
	bool IsErrors = false;
	// Verification W
	if (Config->W.isFix == 1 && Config->W.start == Config->W.end)
	{
		;//OK
	}
	else if (Config->W.isFix == 2 || Config->W.isFix == 0)
	{
		Config->W.isFix = 2;
		if (Config->W.start < Config->W.end && Config->W.NbStep>=1)
		{
			Config->W.step = (Config->W.end - Config->W.start) / ((float)Config->W.NbStep);
		}
		else
		{
			std::cout << "Error W  :  wmax must be sup wmin AND nbStep must be sup 1\n";
			IsErrors = true;
		}
	}
	else
	{
		std::cout << "Error W  :  Error unknow\n";
		IsErrors = true;
	}
	// Verification X
	if (Config->X.isFix == 1 && Config->X.start == Config->X.end)
	{
		;//OK
	}
	else if (Config->X.isFix == 2 || Config->X.isFix == 0)
	{
		Config->X.isFix = 2;
		if (Config->X.start < Config->X.end && Config->X.NbStep>=1)
		{
			Config->X.step = (Config->X.end - Config->X.start) / ((float)Config->X.NbStep);
		}
		else
		{
			std::cout << "Error X  :  xmax must be sup xmin AND nbStep must be sup 1\n";
			IsErrors = true;
		}
	}
	else
	{
		std::cout << "Error X  :  Error unknow\n";
		IsErrors = true;
	}
	// Verification Y
	if (Config->Y.isFix == 1 && Config->Y.start == Config->Y.end)
	{
		;//OK
	}
	else if (Config->Y.isFix == 2 || Config->Y.isFix == 0)
	{
		Config->Y.isFix = 2;
		if (Config->Y.start < Config->Y.end && Config->Y.NbStep>=1)
		{
			Config->Y.step = (Config->Y.end - Config->Y.start) / ((float)Config->Y.NbStep);
		}
		else
		{
			std::cout << "Error Y  :  ymax must be sup ymin AND nbStep must be sup 1\n";
			IsErrors = true;
		}
	}
	else
	{
		std::cout << "Error Y  :  Error unknow\n";
		IsErrors = true;
	}
	// Verification Z
	if (Config->Z.isFix == 1 && Config->Z.start == Config->Z.end)
	{
		;//OK
	}
	else if (Config->Z.isFix == 2 || Config->Z.isFix == 0)
	{
		Config->Z.isFix = 2;
		if (Config->Z.start < Config->Z.end && Config->Z.NbStep>=1)
		{
			Config->Z.step = (Config->Z.end - Config->Z.start) / ((float)Config->Z.NbStep);
		}
		else
		{
			std::cout << "Error Z  :  zmax must be sup zmin AND nbStep must be sup 1\n";
			IsErrors = true;
		}
	}
	else
	{
		std::cout << "Error Z  :  Error unknow\n";
		IsErrors = true;
	}
	// Verification Iter
	if (Config->Iter.isFix == 1)
	{
		if (Config->Iter.min == Config->Iter.max && Config->Iter.max > 0)
		{
			; //OK
		}
		else
		{
			std::cout << "Error Iter  :  value must be sup at 0\n";
			IsErrors = true;
		}
	}
	else if (Config->Iter.isFix == 2 || Config->Iter.isFix == 0)
	{
		Config->Iter.isFix = 2;
		if (Config->Iter.min < Config->Iter.max)
		{
			; // OK
		}
		else
		{
			std::cout << "Error Iter  :  max must be sup min\n";
			IsErrors = true;
		}
	}
	else
	{
		std::cout << "Error Iter  :  Error unknow\n";
		IsErrors = true;
	}
	// Vérification dev
	int count;
	cudaGetDeviceCount(&count);
	if (Config->dev >= 0 && Config->dev < count)
	{
		; //OK
	}
	else
	{
		std::cout << "Error dev :  dev must be between 0 and " << count - 1 << "\n";
		IsErrors = true;
	}
	//Verification filter
	if (Config->filter >= 0 && Config->filter <= 2)
	{
		if (Config->Iter.isFix == 1 && Config->filter == 2)
			Config->filter = 1;
	}
	else
	{
		std::cout << "Error filter :  filter must be between 0 and 2\n";
		IsErrors = true;
	}
	//Verification power
	if (Config->power >= 2.0f && Config->power <= 50.0f)
	{
		; //OK
	}
	else
	{
		std::cout << "Error power :  power must be between 2.0 and 50.0\n";
		IsErrors = true;
	}
	//Verification IsShow
	if (Config->isShow == 0 || Config->isShow == 1)
	{
		; //OK
	}
	else
	{
		std::cout << "Error isShow :  isShow must be between 0 and 1 \n";
		IsErrors = true;
	}
	//Verification rMax
	if (Config->rMax >0.0f)
	{
		; //OK
	}
	else
	{
		std::cout << "Error rMax :  rMax must be sup 0.0 \n";
		IsErrors = true;
	}
	if (IsErrors)
		return -1;

	// creation des fichiers :
	strcpy(Config->nameFile.csv, Config->nameFile.root);
	strcpy(Config->nameFile.histo, Config->nameFile.root);
	strcpy(Config->nameFile.stat, Config->nameFile.root);
	strcpy(Config->nameFile.stl, Config->nameFile.root);
	strcpy(Config->nameFile.txt, Config->nameFile.root);

	strcat(Config->nameFile.csv, ".csv");
	strcat(Config->nameFile.histo, ".histo");
	strcat(Config->nameFile.stat, ".stat");
	strcat(Config->nameFile.stl, ".stl");
	strcat(Config->nameFile.txt, ".txt");

	//Affichage de la config
	std::ofstream file;
	file.open(Config->nameFile.stat);
	file << "Parameters Current : " << "\n";
	file << "				W_start = " << Config->W.start << ", W_end = " << Config->W.end << ", W_Step = " << Config->W.step << ", W_NbStep = " << Config->W.NbStep << "\n";
	file << "				X_start = " << Config->X.start << ", X_end = " << Config->X.end << ", X_Step = " << Config->X.step << ", X_NbStep = " << Config->X.NbStep << "\n";
	file << "				Y_start = " << Config->Y.start << ", Y_end = " << Config->Y.end << ", Y_Step = " << Config->Y.step << ", Y_NbStep = " << Config->Y.NbStep << "\n";
	file << "				Z_start = " << Config->Z.start << ", Z_end = " << Config->Z.end << ", Z_Step = " << Config->Z.step << ", Z_NbStep = " << Config->Z.NbStep << "\n";
	file << "				Root FileOutput = " << Config->nameFile.root << "\n";
	file << "				iterMax = " << Config->Iter.max << "\n";
	file << "				iterMin = " << Config->Iter.min << "\n";
	file << "				rMax = " << Config->rMax << "\n";
	file << "				Filter = " << Config->filter << "\n";
	file << "				Power = " << Config->power << "\n";
	file << "				dev = " << Config->dev << "\n";
	file << "				IsShow = " << Config->isShow << "\n";
	file << "				NbPoints per step = " << Config->NbPointByStep << "\n";
	file << "				ouput File :  " << Config->nameFile.root << "\n";
	file << "cmd for use this configuration: " << "\n";
	file << "               " << argv[0] << "  ";
	if (Config->W.isFix == 1)
		file << Str_wFix << " " << Config->W.start << " ";
	else
	{
		file << Str_wMin << " " << Config->W.start << " ";
		file << Str_wMax << " " << Config->W.end << " ";
		file << Str_wNbStep << " " << Config->W.NbStep << " ";
	}
	if (Config->X.isFix == 1)
		file << Str_xFix << " " << Config->X.start << " ";
	else
	{
		file << Str_xMin << " " << Config->X.start << " ";
		file << Str_xMax << " " << Config->X.end << " ";
		file << Str_xNbStep << " " << Config->X.NbStep << " ";
	}
	if (Config->Y.isFix == 1)
		file << Str_yFix << " " << Config->Y.start << " ";
	else
	{
		file << Str_yMin << " " << Config->Y.start << " ";
		file << Str_yMax << " " << Config->Y.end << " ";
		file << Str_yNbStep << " " << Config->Y.NbStep << " ";
	}
	if (Config->Z.isFix == 1)
		file << Str_zFix << " " << Config->Z.start << " ";
	else
	{
		file << Str_zMin << " " << Config->Z.start << " ";
		file << Str_zMax << " " << Config->Z.end << " ";
		file << Str_zNbStep << " " << Config->Z.NbStep << " ";
	}
	if (Config->Iter.isFix == 1)
		file << Str_IterFix << " " << Config->Iter.max << " ";
	else
	{
		file << Str_IterMin << " " << Config->Iter.min << " ";
		file << Str_IterMax << " " << Config->Iter.max << " ";
	}
	file << Str_dev << " " << Config->dev << " ";
	file << Str_filter << " " << Config->filter << " ";
	file << Str_power << " " << Config->power << " ";
	file << Str_IsShow << " " << Config->isShow << " ";
	file << Str_rMax << " " << Config->rMax << " ";
	file << Str_Out << " " << Config->nameFile.root << " ";

	file << "\n";
	file.close();
	return 1;
}

__device__ void shunte(struct_QD_T *Q_Current)
{
	if (abs(Q_Current->w) <= 0.01)
		Q_Current->w = 0.0;
	if (abs(Q_Current->x) <= 0.01)
		Q_Current->x = 0.0;
	if (abs(Q_Current->y) <= 0.01)
		Q_Current->y = 0.0;
	if (abs(Q_Current->z) <= 0.01)
		Q_Current->z = 0.0;
}

// CUDA kernel to Compute itermax of quaternion
__global__ void kernel(const struct_P_Simulation_T *P_Simulation, int *Tab_Iter)
{
	//int Tempindex = 0;
	struct_QD_T Q_Current;
	double w, x, y, z;
	double N = 0.0;
	int iter = 0;
	//X
	x = __fma_rd(((double)blockIdx.x),P_Simulation->X.step , P_Simulation->X.start);
	//Y
	y = __fma_rd(((double)blockIdx.y),P_Simulation->Y.step , P_Simulation->Y.start);
	//Z
	z = __fma_rd(((double)blockIdx.z),P_Simulation->Z.step , P_Simulation->Z.start);
	//W
	w = __fma_rd(((double)threadIdx.x),P_Simulation->W.step , P_Simulation->W.start);

		CreateQ_By_float(&Q_Current, x, y, z, w);

		for (iter = 0; iter <= P_Simulation->Iter.max; iter++)
		{

			//if (Q_Current.x==NAN || Q_Current.y == NAN || Q_Current.z == NAN || Q_Current.w == NAN)
			//	goto Fin;
			shunte(&Q_Current);

			Get_QPow(&Q_Current, P_Simulation->power);

			shunte(&Q_Current);
			//if (Q_Current.x == NAN || Q_Current.y == NAN || Q_Current.z == NAN || Q_Current.w == NAN)
			//	goto Fin;
			N = Q_Current.x;
			Q_Current.x = __fadd_rn(x,N);

			N = Q_Current.y;
			Q_Current.y = __fadd_rn(y, N);

			N = Q_Current.z;
			Q_Current.z = __fadd_rn(z, N);

			N = Q_Current.w;
			Q_Current.w = __fadd_rn(w, N);

			shunte(&Q_Current);
			/*if (Q_Current.x == NAN || Q_Current.y == NAN || Q_Current.z == NAN || Q_Current.w == NAN)
			goto Fin;*/
			N = norm4d(Q_Current.x, Q_Current.y, Q_Current.z, Q_Current.w);
			//printf("%f , %f, %f, %f, %f\n", N, x, y, z, w);
			if (N == NAN)
				goto Fin;
			if (N> P_Simulation->rMax)
				goto Fin;
		}

Fin:
	if (iter > 0)
		iter--;
	
	int index = blockIdx.x*P_Simulation->X.coef + blockIdx.y*P_Simulation->Y.coef + blockIdx.z*P_Simulation->Z.coef + threadIdx.x*P_Simulation->W.coef;
	//if (index < P_Simulation->max)
		Tab_Iter[index] = iter;// index % 255;
	/*else
		printf("%d > %d", index, P_Simulation->max);*/
}

void printParameterer(struct_P_Simulation_T *P_Simulation)
{
	std::cout << "*****************  P_Simulation Config  *********************" << "\n";
	if (P_Simulation->W.step == 0.0f)
		std::cout << " W Fix = " << P_Simulation->W.start << "\n";
	else
		std::cout << " W start = " << P_Simulation->W.start << " to W end =" << P_Simulation->W.end << " by  W step =" << P_Simulation->W.step << "\n";

	if (P_Simulation->X.step == 0.0f)
		std::cout << " X Fix = " << P_Simulation->X.start << "\n";
	else
		std::cout << " X start = " << P_Simulation->X.start << " to X end =" << P_Simulation->X.end << " by  X step =" << P_Simulation->X.step << "\n";

	if (P_Simulation->Y.step == 0.0f)
		std::cout << " Y Fix = " << P_Simulation->Y.start << "\n";
	else
		std::cout << " Y start = " << P_Simulation->Y.start << " to Y end =" << P_Simulation->Y.end << " by  Y step =" << P_Simulation->Y.step << "\n";

	if (P_Simulation->Z.step == 0.0f)
		std::cout << " Z Fix = " << P_Simulation->Z.start << "\n";
	else
		std::cout << " Z start = " << P_Simulation->Z.start << " to Z end =" << P_Simulation->Z.end << " by  Z step =" << P_Simulation->Z.step << "\n";


	std::cout << "*****************  P_Simulation Config  *********************" << "\n";

}


struct_Stat_float_T compute(struct_P_All_T Config, bool state)
{
std::ofstream filetxt;
std::ofstream file;
	
int Tab_Histo[300];
int  Nbpoint_iter = 0;
int NoConfig = 0;
int NbConfig = max(1, Config.W.NbStep)* max(1, Config.X.NbStep)* max(1, Config.Y.NbStep)* max(1, Config.Z.NbStep);
//Stat
struct_Stat_float_T Stat;
for (int NoW = 0; NoW < max(1,Config.W.NbStep); NoW++)
{
	for (int NoX = 0; NoX < max(1, Config.X.NbStep); NoX++)
	{
		for (int NoY = 0; NoY < max(1, Config.Y.NbStep); NoY++)
		{
			for (int NoZ = 0; NoZ < max(1, Config.Z.NbStep); NoZ++)
			{
				NoConfig++;
				std::cout << "---------------------------------------------------\n";
				std::cout << "Config  " << NoConfig << " sur " << NbConfig << "\n";

				float W = NoW*Config.W.step + Config.W.start;
				float X = NoX*Config.X.step + Config.X.start;
				float Y = NoY*Config.Y.step + Config.Y.start;
				float Z = NoZ*Config.Z.step + Config.Z.start;

				int PasW = Config.NbPointByStep;
				if (Config.W.NbStep == 0)
					PasW = 1;

				int PasX = Config.NbPointByStep;
				if (Config.X.NbStep == 0)
					PasX = 1;

				int PasY = Config.NbPointByStep;
				if (Config.Y.NbStep == 0)
					PasY = 1;

				int PasZ = Config.NbPointByStep;
				if (Config.Z.NbStep == 0)
					PasZ = 1;

				//Taille de tableau
				int max =PasZ*PasY*PasX*PasW;

				if (Config.isShow)
					std::cout << "cudaMallocManaged Config  -->  Start" << "\n";

				struct_P_Simulation_T *P_Simulation;
				int *Tab_Iter;
				// Allocate Unified Memory -- accessible from CPU or GPU
				cudaMallocManaged(&P_Simulation, sizeof(struct_P_Simulation_T));
				cudaMallocManaged(&Tab_Iter, max * sizeof(int));
				if (Config.isShow)
					std::cout << "cudaMallocManaged Config  -->  End " << "\n";

				if (Config.isShow)
					std::cout << "P_Simulation Config  -->  Start" << "\n";
				// Pramatrage de W
				P_Simulation->W.start = W;
				P_Simulation->W.end = W + Config.W.step;
				P_Simulation->W.NbStep = PasW;
				if (PasW == 1)
					P_Simulation->W.step = 0.0f;
				else
					P_Simulation->W.step = (Config.W.step) / (PasW - 1);
				P_Simulation->W.coef = PasZ*PasY*PasX;

				// Pramatrage de X
				P_Simulation->X.start = X;
				P_Simulation->X.end = X + Config.X.step;
				P_Simulation->X.NbStep = PasX;
				if (PasX == 1)
					P_Simulation->X.step = 0.0f;
				else
					P_Simulation->X.step = (Config.X.step) / (PasX - 1);
				P_Simulation->X.coef = PasZ*PasY;

				// Pramatrage de Y
				P_Simulation->Y.start = Y;
				P_Simulation->Y.end = Y + Config.Y.step;
				P_Simulation->Y.NbStep = PasY;
				if (PasY == 1)
					P_Simulation->Y.step = 0.0f;
				else
					P_Simulation->Y.step = (Config.Y.step) / (PasY - 1);
				P_Simulation->Y.coef = PasZ;

				// Pramatrage de Z
				P_Simulation->Z.start = Z;
				P_Simulation->Z.end = Z + Config.Z.step;
				P_Simulation->Z.NbStep = PasZ;
				if (PasZ == 1)
					P_Simulation->Z.step = 0.0f;
				else
					P_Simulation->Z.step = (Config.Z.step) / (PasZ - 1);
				P_Simulation->Z.coef = 1;

				printParameterer(P_Simulation);

				//Stat Step

				Stat.Wstep = P_Simulation->W.step;
				Stat.Xstep = P_Simulation->X.step;
				Stat.Ystep = P_Simulation->Y.step;
				Stat.Zstep = P_Simulation->Z.step;


				//Parametrage Iter
				P_Simulation->Iter.max = Config.Iter.max;
				P_Simulation->Iter.min = Config.Iter.min;

				//Parametrage Power
				P_Simulation->power = Config.power;

				//Parametrage Rmax
				P_Simulation->rMax = Config.rMax;

				//Parametrage max
				P_Simulation->max = max;

				if (Config.isShow)
					std::cout << "P_Simulation Config  -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Tab_Iter and Tab_Histo Init  -->  Start" << "\n";
				for (int i = 0; i < max; i++)
					Tab_Iter[i] = 0;

				for (int i = 0; i <= Config.Iter.max; i++)
					Tab_Histo[i] = 0;
				if (Config.isShow)
					std::cout << "Tab_Iter and Tab_Histo Init -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Compude GPU -->  Start" << "\n";
				dim3 grid(PasX, PasY, PasZ);
				dim3 block(PasW, 1, 1);
				kernel << <grid, block >> >(P_Simulation, Tab_Iter);
				if (Config.isShow)
					std::cout << "Compude GPU -->  End" << "\n";


				if (Config.isShow)
					std::cout << "cudaDeviceSynchronize-->  Start" << "\n";
				cudaDeviceSynchronize();
				if (Config.isShow)
					std::cout << "cudaDeviceSynchronize -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Analyzer Simulation -->  Start" << "\n";
				Nbpoint_iter = 0;
				for (int i = 0; i < max; i++)
				{
					if (Tab_Iter[i] > 0)
						Nbpoint_iter++;
					Tab_Histo[Tab_Iter[i]]++;
				}
				if (Config.isShow)
				{
					std::cout << "Nb point Nbpoint_iter = " << Nbpoint_iter << "\n";
					std::cout << "Soit  :  " << (float)(Nbpoint_iter / ((float)max / 10000.0f)) / 100.0f << "%  soit " << Nbpoint_iter << "pt sur " << max << "pt \n";
					std::cout << "Analyzer Simulation -->  End" << "\n";
				}
				if (state == true)
				{
					if (Config.isShow)
						std::cout << "Write Histogram -->  Start" << "\n";
					file.open(Config.nameFile.histo, std::ofstream::out | std::ofstream::app);
					file << NoConfig << ";";
					for (int i = 0; i <= Config.Iter.max; i++)
						file << Tab_Histo[i] << ";";
					file << "\n";
					file.close();
					if (Config.isShow)
						std::cout << "Write Histogram -->  End" << "\n";
				}
				


				file.open(Config.nameFile.csv, std::ofstream::out | std::ofstream::app);
				filetxt.open(Config.nameFile.txt, std::ofstream::out | std::ofstream::app);


				for (int i = 0; i < max; i++)
				{
					int j = i;

					//W
					int iW = 0;
					if (PasW > 1)
					{
						iW = j / P_Simulation->W.coef;
					}
					//printf("index = %d  - Z Tempindex = %d \n", i, Tempindex);
					float w = (float)iW*P_Simulation->W.step + P_Simulation->W.start;
					// on retranche 
					j -= iW*P_Simulation->W.coef;

					//X
					int iX = 0;
					if (PasX > 1)
					{
						iX = j / P_Simulation->X.coef;
					}
					float x = (float)iX*P_Simulation->X.step + P_Simulation->X.start;
					// on retranche 
					j -= iX*P_Simulation->X.coef;

					//Y
					int iY = 0;
					if (PasY > 1)
					{
						iY = j / P_Simulation->Y.coef;
					}
					float y = (float)iY*P_Simulation->Y.step + P_Simulation->Y.start;
					// on retranche 
					j -= iY*P_Simulation->Y.coef;

					//Z
					int iZ = 0;
					if (PasZ > 1)
					{
						iZ = j / P_Simulation->Z.coef;
					}
					float z = (float)iZ*P_Simulation->Z.step + P_Simulation->Z.start;
					// on retranche 
					j -= iZ*P_Simulation->Z.coef;


					int iter = Tab_Iter[i];
					if ((iter >= Config.Iter.min && Config.Iter.isFix == 2) || (iter == Config.Iter.min && Config.Iter.isFix == 1))
					{
						int filter = Config.filter;
						if (state == true)
						{
							if (FilterQ_H(&filter, &iX, &iY, &iZ, &iW, iter, P_Simulation))
							{
								file << x << ";" << y << ";" << z << ";" << w << ";" << iter << ";\n";
								int NbDim = 0;
								if (PasW > 1)
								{
									filetxt << w << ";";
									NbDim++;
								}

								if (PasX > 1)
								{
									filetxt << x << ";";
									NbDim++;
								}

								if (PasY > 1)
								{
									filetxt << y << ";";
									NbDim++;
								}

								if (PasZ > 1)
								{
									filetxt << z << ";";
									NbDim++;
								}

								if (NbDim >= 3)
									filetxt << iter << "\n";
								else
									filetxt << ((float)iter) / ((float)Config.Iter.max)*3.0f << "\n";
							}
						}
						if (w > Stat.Wmax)
							Stat.Wmax = w;
						if (w < Stat.Wmin)
							Stat.Wmin = w;

						if (x > Stat.Xmax)
							Stat.Xmax = x;
						if (x < Stat.Xmin)
							Stat.Xmin = x;

						if (y > Stat.Ymax)
							Stat.Ymax = y;
						if (y < Stat.Ymin)
							Stat.Ymin = y;

						if (z > Stat.Zmax)
							Stat.Zmax = z;
						if (z < Stat.Zmin)
							Stat.Zmin = z;

						Stat.NbPoint++;
					}

				}
				file.close();
				filetxt.close();

				if (Config.isShow)
					std::cout << "Clear Mem + Reste  -->  Start" << "\n";
				cudaFree(P_Simulation);
				cudaFree(Tab_Iter);
				//cudaDeviceReset();
				if (Config.isShow)
					std::cout << "Clear Mem + Reste  -->  End" << "\n";
			}
		}
	}
}

return Stat;

}


int main(int argc, char *argv[])
{
	struct_P_All_T Config;
    int state= menu( argc, argv, &Config);

	if (state != 1)
		return state;

	/********  Clear File ************/
	std::ofstream filetxt;
	filetxt.open(Config.nameFile.txt);
	filetxt.close();

	std::ofstream file;
	file.open(Config.nameFile.csv);
	file << "X;Y;Z;W;iter;\n";
	file.close();


	file.open(Config.nameFile.histo);
	file << "index;";
	for (int i = 0; i <= Config.Iter.max; i++)
		file << i << ";";
	file << "\n";
	file.close();

	//Affiche Stat
	FILE *fileman;
	char line[1000];
	fileman = fopen(Config.nameFile.stat, "r");
	if (fileman != NULL)
	{
		while (std::fgets(line, 1000, fileman))
		{
			std::cout << line;
		}
		fclose(fileman);
	}
	else
		std::cout << "Error  : "<< Config.nameFile.stat <<" no found" << "\n";
	cudaSetDevice(Config.dev);


	struct_Stat_float_T Stat;

	if(Config.modeAM==false)
		Stat =compute(Config, true);
	else
	{
		bool Dim[] = {false,false ,false ,false };
		
		if (Config.W.isFix == 2)
			Dim[0] = true;
		if (Config.X.isFix == 2)
			Dim[1] = true;
		if (Config.Y.isFix == 2)
			Dim[2] = true;
		if (Config.Z.isFix == 2)
			Dim[3] = true;

		int NbDim = 0;

		for (int i = 0; i <= 3; i++)
		{
			if (Dim[i]==true)
			{
				NbDim++;
				switch (i)
				{
				case 0:
					Config.W.isFix = 2;
					Config.W.start = -20.0f;
					Config.W.end = 20.0f;
					Config.W.NbStep = 2;
					Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
					break;
				case 1:
					Config.X.isFix = 2;
					Config.X.start = -20.0f;
					Config.X.end = 20.0f;
					Config.X.NbStep = 2;
					Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
					break;
				case 2:
					Config.Y.isFix = 2;
					Config.Y.start = -20.0f;
					Config.Y.end = 20.0f;
					Config.Y.NbStep = 2;
					Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
					break;
				case 3:
					Config.Z.isFix = 2;
					Config.Z.start = -20.0f;
					Config.Z.end = 20.0f;
					Config.Z.NbStep = 2;
					Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
					break;
				default:
					break;
				}
			}	
		}
		//si et selement si il y a 3 dimentions
		if (NbDim == 3)
		{
			float pas = 10;
			float OldPas = 20;
			int nb_iter = 10;
			for (int iter = 0; iter <= nb_iter; iter++)
			{
				Stat = compute(Config,false);
				// rétro action
				for (int i = 0; i <= 3; i++)
				{
					if (Dim[i] == true)
					{
						switch (i)
						{
						case 0:
							if (Stat.Wmin == Config.W.start)
							{
								Config.W.start -= OldPas;
								if (nb_iter != i)
									Config.W.start += pas;
							}
							else if(nb_iter!=i)
								Config.W.start +=pas;

							if (Stat.Wmax == Config.W.end)
							{
								Config.W.end += OldPas;
								if (nb_iter != i)
									Config.W.end -= pas;
							}								
							else if(nb_iter != i)
								Config.W.end -= pas;
							Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
							break;
						case 1:
							if (Stat.Xmin == Config.X.start)
							{
								Config.X.start -= OldPas;
								if (nb_iter != i)
									Config.X.start += pas;
							}
							else if (nb_iter != i)
								Config.X.start += pas;
							if (Stat.Xmax == Config.X.end)
							{
								Config.X.end += OldPas;
								if (nb_iter != i)
									Config.X.end -= pas;
							}
							else if (nb_iter != i)
								Config.X.end -= pas;
							Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
							break;
						case 2:
							if (Stat.Ymin == Config.Y.start)
							{
								Config.Y.start -= OldPas;
								if (nb_iter != i)
									Config.Y.start += pas;
							}
							else if (nb_iter != i)
								Config.Y.start += pas;
							if (Stat.Ymax == Config.Y.end)
							{
								Config.Y.end += OldPas;
								if (nb_iter != i)
									Config.Y.end -= pas;
							}
							else if (nb_iter != i)
								Config.Y.end -= pas;
							Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
							break;
						case 3:
							if (Stat.Zmin == Config.Z.start)
							{
								Config.Z.start -= OldPas;
								if (nb_iter != i)
									Config.Z.start += pas;
							}
							else if (nb_iter != i)
								Config.Z.start += pas;
							if (Stat.Zmax == Config.Z.end)
							{
								Config.Z.end += OldPas;
								if (nb_iter != i)
									Config.Z.end -= pas;
							}
							else if (nb_iter != i)
								Config.Z.end -= pas;
							Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
							break;
						default:
							break;
						}
					}

				}

				file.open(Config.nameFile.stat, std::ofstream::out | std::ofstream::app);
				file << "Statistiques : \n";
				file << "				X min = " << Stat.Xmin << "\n";
				file << "				X max = " << Stat.Xmax << "\n";
				file << "				Y min = " << Stat.Ymin << "\n";
				file << "				Y max = " << Stat.Ymax << "\n";
				file << "				Z min = " << Stat.Zmin << "\n";
				file << "				Z max = " << Stat.Zmax << "\n";
				file << "				W min = " << Stat.Wmin << "\n";
				file << "				W max = " << Stat.Wmax << "\n";
				file << "				X step = " << Stat.Xstep << "\n";
				file << "				Y step = " << Stat.Ystep << "\n";
				file << "				Z step = " << Stat.Zstep << "\n";
				file << "				W step = " << Stat.Wstep << "\n";
				file << "				NbPoint plot = " << Stat.NbPoint << "\n";
				file.close();


				OldPas = pas;
				pas = pas / 2.0f;
			}

			for (int i = 0; i <= 3; i++)
			{
				if (Dim[i] == true)
				{
					switch (i)
					{
					case 0:
						Config.W.NbStep = 5;
						Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
						break;
					case 1:
						Config.X.NbStep = 5;
						Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
						break;
					case 2:
						Config.Y.NbStep = 5;
						Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
						break;
					case 3:
						Config.Z.NbStep = 5;
						Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
						break;
					default:
						break;
					}
				}

			}
			if(Stat.NbPoint!=0)
				Stat = compute(Config, true);
		}
		else
		{
			std::cout << "Error  mode: mode auto only 3 dim\n";
			return -1;
		}

	}
	
	file.open(Config.nameFile.stat, std::ofstream::out | std::ofstream::app);
	file << "Statistiques : \n";
	file << "				X min = " << Stat.Xmin << "\n";
	file << "				X max = " << Stat.Xmax << "\n";
	file << "				Y min = " << Stat.Ymin << "\n";
	file << "				Y max = " << Stat.Ymax << "\n";
	file << "				Z min = " << Stat.Zmin << "\n";
	file << "				Z max = " << Stat.Zmax << "\n";
	file << "				W min = " << Stat.Wmin << "\n";
	file << "				W max = " << Stat.Wmax << "\n";
	file << "				X step = " << Stat.Xstep << "\n";
	file << "				Y step = " << Stat.Ystep << "\n";
	file << "				Z step = " << Stat.Zstep << "\n";
	file << "				W step = " << Stat.Wstep << "\n";
	file << "				NbPoint plot = " << Stat.NbPoint << "\n";
	file.close();




	return 0;//Fin du programme
}