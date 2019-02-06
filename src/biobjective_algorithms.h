
using namespace std;


#ifdef _WIN32
	#ifdef _WIN64
		//WINDOWS 64 bits
	#elif
		//Exlusivo windows 32 bits
	#endif
#include <direct.h>
class TIEMPO {
	clock_t  t;
	int tipo;

public:
	TIEMPO() {	//Crea objeto tiempo y lo inicializa al reloj
		t = clock();
		tipo = 0;
	} //Constructor
	
	void TIEMPO::init() { //Reinicializa objeto a cero
		t = clock() - clock();
	}

	void TIEMPO::acum() { //Acumula tiempo desde la ultima llamada al objeto
		t += clock() - t;
	}
	void TIEMPO::interval() { //Tiempo transcurrido entre las dos ultimas llamadas al objeto. Ese es el nuevo valor del objeto.
		t = clock() - t;
	}

	double TIEMPO::value() { //Transforma a segundos (tipo double) la ultima llamada al objeto
		return ((double)(t / double(CLOCKS_PER_SEC)));
	}

	int gettype(){
			return tipo;
	}

	void chgtype(int newt){
		//No hace nada este metodo bajo windows
		tipo = newt;
	}

	~TIEMPO() { } //Destructor
};

#elif __linux__
#include <sys/stat.h> //Para mkdir en Linux

class TIEMPO {
	timespec t;
	int tipo;

private:
	timespec diff(timespec start, timespec end){
		timespec temp;
		if ((end.tv_nsec - start.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec - start.tv_sec - 1;
			temp.tv_nsec = 1000000000.0 + end.tv_nsec - start.tv_nsec;
		}
		else {
			temp.tv_sec = end.tv_sec - start.tv_sec;
			temp.tv_nsec = end.tv_nsec - start.tv_nsec;
		}
		return temp;
	}

	timespec add(timespec accum, timespec now) {
		timespec temp;
		if (accum.tv_nsec + now.tv_nsec >= 1E9) {
			temp.tv_sec = accum.tv_sec + now.tv_sec + 1;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec - 1E9;
		}
		else {
			temp.tv_sec = accum.tv_sec + now.tv_sec;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec;
		}
		return temp;
	}


public:

	TIEMPO() {	//Crea objeto tiempo y lo inicializa a uno de los cuatro tipos
		tipo = CLOCK_THREAD_CPUTIME_ID;
		clock_gettime(tipo, &t);
	} //Constructor


	void init() { //Reinicializa objeto a cero
		t.tv_sec = t.tv_sec - t.tv_sec;
		t.tv_nsec = t.tv_nsec - t.tv_nsec;
		//printf("\nt_ini_s = %0.1f  t_init_nsec = %0.1f" , (double) t.tv_sec , (double) t.tv_nsec);
	}


	void acum() { //Acumula tiempo desde la �ltima llamada al objeto
		timespec newt , aux;
		//printf("\nnewt = %2.5f, aux = %2.5f , t=%2.5f" , (double) newt.tv_sec, (double) aux.tv_sec , (double) t.tv_sec);
		//printf("\nt time %0.2f, %0.2f " , (double) t.tv_sec, (double) t.tv_nsec );
		clock_gettime(tipo, &newt);
		//printf("\nclock get time %0.2f, %0.2f " , (double) newt.tv_sec, (double) newt.tv_nsec );
		aux = diff(t , newt);
		//printf("\nclock get time dif %0.2f, %0.2f " , (double) aux.tv_sec, (double) aux.tv_nsec );
		t = add(t, aux);

	}
	void interval() { //Tiempo transcurrido entre las dos �ltimas llamadas al objeto. Ese es el nuevo valor del objeto.
		timespec newt;
		clock_gettime(tipo, &newt);
		t = diff(t , newt);
	}

	double value() { //Transforma a segundos (tipo double) la �ltima llamada al objeto
		//printf("\n%2.5f Sec" , (double) t.tv_sec);
		return (t.tv_sec + t.tv_nsec / 1E9);
	}

	int gettype(){
		//printf("\n Tipo = %d" , tipo);
		return tipo;
	}

	void chgtype(int newt){
		//0 = CLCOK_REALTIME
		//1 = CLOCK_MONOTONIC
		//2 = CLOCK_PROCESS_CPUTIME_ID
		//3 = CLOCK_THREAD_CPUTIME_ID
		tipo = newt;
		//printf("\n nuevo tipo %d" , tipo);
	}



	~TIEMPO() { } //Destructor
};

#elif __APPLE__
	//APPLE
#elif __unix__
	//Unix
#elif defined(_POSIX_VERSION)
	// POSIX
#else
	#   error "Unknown compiler"
#endif

int Obtener_no_dominados_lexicograficos(BOILP *P1, int *indice_col, double **z1, double **z2) {
	CPXENVptr env = *P1->env;
	CPXLPptr lp = *P1->lp;
	list<solution> *Efic = P1->Efic;
	double **F = P1->F;
	double hipervolumen = P1->hypervolume;
	int n = P1->n;
	int numiteraciones = P1->n_iterations;
	double tiempo_max = P1->max_time;
	int m = P1->m;
	int matbeg = 0;

	double Lex1[2], Lex2[2];

	int sta;
	CPXLPptr lp2 = CPXcloneprob(env, lp, &sta);

	//FIRST LEXICOGRAPHICAL OPTIMAL SOLUTION
	int status = CPXchgobj(env, lp, n, indice_col, F[0]);
	status = CPXmipopt(env, lp); //min (f1)
	if (!CPXgetstat(env, lp)) { printf("Problema infactible");		return 0; }
	status = CPXgetobjval(env, lp, &(Lex1[0])); //save value

	status = CPXaddrows(env, lp, 0, 1, n, &(Lex1[0]), "L", &matbeg, indice_col, F[0], NULL, NULL); //Add constraint f1(z1) <= e
	status = CPXchgobj(env, lp, n, indice_col, F[1]);
	status = CPXmipopt(env, lp); //min (f2) s.t. f1(x)<= f1(z1)
	if (!CPXgetstat(env, lp)) { printf("Problema infactible");		return 0; }
	status = CPXgetobjval(env, lp, &(Lex1[1]));
	status = CPXdelrows(env, lp, m, m);

	//SECOND LEXICOGRAPHICAL OPTIMAL SOLUTION
	status = CPXchgobj(env, lp2, n, indice_col, F[1]);
	status = CPXmipopt(env, lp2); //min (f2)
	status = CPXgetobjval(env, lp2, &(Lex2[1]));
	status = CPXaddrows(env, lp2, 0, 1, n, &(Lex2[1]), "L", &matbeg, indice_col, F[1], NULL, NULL);
	status = CPXchgobj(env, lp2, n, indice_col, F[0]);
	status = CPXmipopt(env, lp2); //min (f1) s.t. f2(x)<= f2(z2)
	status = CPXgetobjval(env, lp2, &(Lex2[0]));
	status = CPXdelrows(env, lp2, m, m);

	if (comparar_vectores(Lex1, Lex2, 2) == 1) { //In case of the optimum is the ideal point
		(*z1)[0] = (*z2)[0] = round(Lex1[0]);
		(*z1)[1] = (*z2)[1] = round(Lex1[1]);
		return 1;
	}
	else { //There are two lexicographical optimal values
		(*z1)[0] = round(Lex1[0]);		(*z1)[1] = round(Lex1[1]);
		(*z2)[0] = round(Lex2[0]);		(*z2)[1] = round(Lex2[1]);
		return 2;
	}
}

int  Obtener_lexicograficos(BOILP *P1, int *indice_col, double **x1, double **x2, double **z1, double **z2, double *t1 , double *t2, TIEMPO t_ref) {
	CPXENVptr env = *P1->env;
	CPXLPptr lp = *P1->lp;
	list<solution> *Efic = P1->Efic;
	double **F = P1->F;
	double hipervolumen = P1->hypervolume;
	int n = P1->n;
	int numiteraciones = P1->n_iterations;
	double tiempo_max = P1->max_time;
	int m = P1->m;
	int i;
	int matbeg = 0;

	double Lex1[2], Lex2[2];

	int sta;
	CPXLPptr lp2 = CPXcloneprob(env, lp, &sta);

	//FIRST LEXICOGRAPHICAL OPTIMAL POINT
	int status = CPXchgobj(env, lp, n, indice_col, F[0]);
	status = CPXmipopt(env, lp); //min (f1)
	if (!CPXgetstat(env, lp)) { printf("Problema infactible");		return 0; }
	status = CPXgetobjval(env, lp, &(Lex1[0])); //save value

	status = CPXaddrows(env, lp, 0, 1, n, &(Lex1[0]), "L", &matbeg, indice_col, F[0], NULL, NULL); //Add constraint f1(z1) <= e
	status = CPXchgobj(env, lp, n, indice_col, F[1]);
	status = CPXmipopt(env, lp); //min (f2) s.t. f1(x)<= f1(z1)
	if (!CPXgetstat(env, lp)) { printf("Problema infactible");		return 0; }
	status = CPXgetobjval(env, lp, &(Lex1[1]));

	//Save solution
	for (i = 0; i < n; i++) status = CPXgetx(env, lp, *x1, 0, n - 1);
	status = CPXdelrows(env, lp, m, m);
	t_ref.acum();
	*t1 = t_ref.value();

	//SECOND LEXICOGRAPHICAL OPTIMAL POINT
	status = CPXchgobj(env, lp2, n, indice_col, F[1]);
	status = CPXmipopt(env, lp2); //min (f2)
	status = CPXgetobjval(env, lp2, &(Lex2[1]));
	status = CPXaddrows(env, lp2, 0, 1, n, &(Lex2[1]), "L", &matbeg, indice_col, F[1], NULL, NULL);
	status = CPXchgobj(env, lp2, n, indice_col, F[0]);
	status = CPXmipopt(env, lp2); //min (f1) s.t. f2(x)<= f2(z2)
	status = CPXgetobjval(env, lp2, &(Lex2[0]));

	//Save solution
	for (i = 0; i < n; i++) status = CPXgetx(env, lp2, *x2, 0, n - 1);
	status = CPXdelrows(env, lp2, m, m);
	t_ref.acum();
	*t2 = t_ref.value();

	if (comparar_vectores(Lex1, Lex2, 2) == 1) { //The optimum is the ideal point
		(*z1)[0] = (*z2)[0] = round(Lex1[0]);
		(*z1)[1] = (*z2)[1] = round(Lex1[1]);
		return 1;
	}
	else { //There are two lexicographical optimal solutions
		(*z1)[0] = round(Lex1[0]);		(*z1)[1] = round(Lex1[1]);
		(*z2)[0] = round(Lex2[0]);		(*z2)[1] = round(Lex2[1]);
		return 2;
	}
}

void TuneProblem(BOILP *P1) {
	int tunestat = 1;
	double d1, d2, d3, d4, d5;
	int i1;

	int dp[6] = { CPX_PARAM_EPINT, CPX_PARAM_EPGAP , CPX_PARAM_EPAGAP , CPX_PARAM_EPRHS, CPX_PARAM_EPOPT };
	int ip[1] = { CPX_PARAM_ADVIND };

	CPXgetdblparam(*P1->env, CPX_PARAM_EPINT, &d1);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPGAP, &d2);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPAGAP, &d3);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPRHS, &d4);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPOPT, &d5);
	CPXgetintparam(*P1->env, CPX_PARAM_ADVIND, &i1);

	double dpv[6] = {d1,d2,d3,d4,d5};
	int ipv[1] = {i1};
	CPXtuneparam(*P1->env, *P1->lp, 1, ip, ipv, 5, dp, dpv, 0, 0, 0, &tunestat);
	
}

int CALCULA_LEXICOGRAFICOS(BOILP *P1, int *indice_col, solution *sol, TIEMPO t_ref) {
	list<solution> *Efic = P1->Efic;
	int n = P1->n;

	double *x_1 = (double *)malloc(n * sizeof(double)); //First solution (NULL if empty)
	double *y_1 = (double *)malloc(2 * sizeof(double));
	double *x_2 = (double *)malloc(n * sizeof(double)); //Second solution (NULL if empty)
	double *y_2 = (double *)malloc(2 * sizeof(double));

	double t1, t2;
	int num_opt_lexicograficos = Obtener_lexicograficos(P1, indice_col, &x_1, &x_2, &y_1, &y_2, &t1 , &t2 , t_ref);
	P1->n_iterations += num_opt_lexicograficos * 2;

	if (num_opt_lexicograficos == 0) {
		free(x_1); free(x_2); free(y_1); free(y_2);
		return 0;
	}
	else if (num_opt_lexicograficos == 1) {
		free(x_2); free(y_2);
		sol->dimension_x = n; sol->dimension_y = 2; sol->vector_x = x_1; sol->vector_y = y_1; sol->time = t2;
		Efic->push_front(*sol);
		return 0;
	}
	else {
		sol->dimension_x = n; sol->dimension_y = 2;		
		sol->vector_x = x_1; sol->vector_y = y_1;		sol->time = t1;		sol->acum_hyper = 0.0;			Efic->push_front(*sol);
		sol->vector_x = x_2; sol->vector_y = y_2;		sol->time = t2;		sol->acum_hyper = 0.0;			Efic->push_back(*sol);
		if ((y_1[0] + 1 == y_2[0]) || (y_2[1] + 1 == y_1[1]))	return 0; 
	}
	return 1;
}

void Crear_dos_nuevas_cajas_a_explorar(caja_bidimensional *box, double *nuevo_y, caja_bidimensional *newb1, caja_bidimensional *newb2, int opcion) {
	double y1, x1, area1, area2;
	double z11 = box->z1[0];	double z12 = box->z1[1];
	double z21 = box->z2[0];	double z22 = box->z2[1];
	double z31 = nuevo_y[0];	double z32 = nuevo_y[1];

	if (opcion == 0) { //We don't take into account the area of a box
		area1 = area2 = 0.0;
	}
	else if (opcion == 1) { //We take into account the area of the box
		area1 = (z12 - z32) * (z31 - z11);
		area2 = (z32 - z22) * (z21 - z31);
	}
	else if (opcion == 2) { //We take into account the area of the non explored zone of the box (only for hybrid)
		double lambda1 = z12 - z22;	double lambda2 = z21 - z11;
		double z0 = lambda1 * z11 + lambda2 * z12; //Tambien z0 = lambda1 * z21 + lambda2 * z22
		 //Line joining extremes of the box:     y = Mx + NN
		double M = (z22 - z12) / (z21 - z11);	//Also M = -lamba1 / lambda2
		double zz = lambda1 * z31 + lambda2 * z32;
		double NN = z32 - M * z31;

		if (zz >= z0) { //The two new regions will be triangles
			x1 = 1.0 / M * (z12 - NN);
			y1 = M * z21 + NN;
			area1 = (z31 - x1) * (z12 - z32) / 2.0;
			area2 = (z21 - z31) * (z32 - y1) / 2.0;
		}
		else { //The two new regions will be trapezoids
			y1 = M * z11 + NN;
			x1 = 1.0 / M * (z22 - NN);
			area1 = ((z12 - z32) + (z12 - y1)) / 2.0 * (z31 - z11);
			area2 = ((z21 - z31) + (z21 - x1)) / 2.0 * (z32 - z22);
		}
	}
	else if (opcion == 3) { //We take into account the lenght of x-axis
		area1 = z31 - z11;
		area2 = z21 - z31;
	}

	//Create two new boxes
	newb1->value = area1;
	newb2->value = area2;

	if ((z31 - z11 == 1) || (z12 - z32 == 1))	newb1->z1 = newb1->z2 = NULL;
	else										newb1->z1 = box->z1;		newb1->z2 = nuevo_y;
	
	if ((z32 - z22 == 1) || (z21 - z31 == 1)) 	newb2->z1 = newb2->z2 = NULL;
	else										newb2->z1 = nuevo_y;		newb2->z2 = box->z2;

}

void crear_problema_tipo_hybrid(CPXENVptr env, CPXLPptr lp, CPXLPptr *newproblem, int n, int *indice_col, double **F) {
	int status;
	*newproblem = CPXcloneprob(env, lp, &status);

	int m = CPXgetnumrows(env, lp);
	int matbeg = 0;
	double rhs = 0.0;
	CPXaddrows(env, *newproblem, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL);
	CPXaddrows(env, *newproblem, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL);
}

void crear_problema_tipo_tchebycheff(CPXENVptr env, CPXLPptr lp, CPXLPptr *newproblem, int n, int *indice_col, double **F) {
	int status;
	*newproblem = CPXcloneprob(env, lp, &status);

	int m = CPXgetnumrows(env, lp);
	int matbeg = 0;
	double rhs = 0.0;

	int ind = n;
	double val = 1.0;

	CPXaddrows(env, *newproblem, 1, 1, 0, &rhs, "L", &matbeg, &matbeg, &rhs, NULL, NULL); //constraint lambda >= w1(z1 - s1) Incomplete yet. Add lambda variable
	CPXaddrows(env, *newproblem, 0, 1, 0, &rhs, "L", &matbeg, &matbeg, &rhs, NULL, NULL); //constraint lambda >= w2(z2 - s2) Incomplete yet
	
	CPXaddrows(env, *newproblem, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL); //constraint  f1 <= z21 Incomplete yet
	CPXaddrows(env, *newproblem, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL); //constraint f2 <= z12 Incomplete yet
	
	CPXchgcoef(env, *newproblem, m, n, -1.0); //lambda coefficient in the first constraint
	CPXchgcoef(env, *newproblem, m + 1, n, -1.0); //lambda coefficient in the second constraint 
	CPXchgobj(env, *newproblem, 1, &ind, &val); //lambda coefficient in the objective function
	
}

void cargar_problemas_lp(BOILP *P1, CPXLPptr *lp1, CPXLPptr *lp2, std::string algoritmo1, std::string algoritmo2, int *indice_col) {

	if (algoritmo1 == "hybrid")
		crear_problema_tipo_hybrid(*P1->env, *P1->lp, lp1, P1->n, indice_col, P1->F);
	else if (algoritmo1 == "tchebycheff")
		crear_problema_tipo_tchebycheff(*P1->env, *P1->lp, lp1, P1->n, indice_col, P1->F);

	if (algoritmo2 == "hybrid")
		crear_problema_tipo_hybrid(*P1->env, *P1->lp, lp2, P1->n, indice_col, P1->F);
	else if (algoritmo2 == "tchebycheff")
		crear_problema_tipo_tchebycheff(*P1->env, *P1->lp, lp2, P1->n, indice_col, P1->F);
}

bool ordenar_lex_2(solution v1, solution v2) { //Arrange lexicographically 2-dimension vectors
	double y1 = v1.vector_y[0];
	double y2 = v2.vector_y[0];

	if (y1 < y2)
		return true;
	else if (y1 == y2) {
		if (v1.vector_y[1] < v2.vector_y[1])
			return true;
	}
	return false;
}

void busca_siguiente_algoritmo_a_utilizar2(caja_bidimensional *Box, caja_bidimensional *newb1, caja_bidimensional *newb2, double *nuevo_y) {
	double lambda1 = Box->z1[1] - Box->z2[1];
	double lambda2 = Box->z2[0] - Box->z1[0];
	double objcaja = lambda1 * Box->z1[0] + lambda2 * Box->z1[1]; 
	double objsol = lambda1 * nuevo_y[0] + lambda2 * nuevo_y[1];

	if (objsol <= objcaja) {
		newb1->algorithm_to_use = 1;
		newb2->algorithm_to_use = 1;
		return;
	}

	double dx1 = nuevo_y[0] - Box->z1[0];
	double dx2 = Box->z2[0]-nuevo_y[0];
	double dy1 = Box->z1[1] - nuevo_y[1];
	double dy2 = nuevo_y[1] - Box->z2[1];

	double xx = Box->z2[0] - Box->z1[0];
	double yy = Box->z1[1] - Box->z2[1];


	if (newb1->z1 != NULL) {
		if ((dx1 < 0.25 * xx) || (dx2 < 0.25 *xx)) {
			newb1->algorithm_to_use = 2;
		}
		else {
			newb1->algorithm_to_use = 1;
		}
	}
	if (newb2->z1 != NULL) {
		if ((dy1 < 0.25 * yy) || (dy2 < 0.25 * yy)) {
			newb1->algorithm_to_use = 2;
		}
		else {
			newb1->algorithm_to_use = 1;
		}
		
	}

}

void Introducir_caja_en_lista_segun_valor(list<caja_bidimensional> *L, caja_bidimensional box) {
	//Arrange boxes in the list according value
	list<caja_bidimensional>::iterator it;

	if ((box.z2[0] - box.z1[0] > 1) && (box.z1[1] - box.z2[1] > 1)) {
		if (L->size() == 0) {
			L->push_front(box);
			return;
		}
		it = L->begin();
		while (it != L->end()) {
			if (it->value > box.value)
				++it;
			else {
				L->insert(it, box);
				return;
			}
		}
		L->push_back(box);
	}

}

void Inicializar_caja_a_explorar(caja_bidimensional *box, double *elem1, double *elem2, int opcion, int algorithm_to_use) {
	box->z1 = elem1;
	box->z2 = elem2;

	if (opcion == 0) { //No value for the box
		box->value = 0.0;
	}
	else if ((opcion == 1) || (opcion == 2)) { //The value is the area
		box->value = (box->z2[0] - box->z1[0]) * (box->z1[1] - box->z2[1]);
	}
	box->algorithm_to_use = 1;
}

void Inicializar_caja_a_explorar(caja_bidimensional *box, double *elem1, double *elem2, int opcion) {
	box->z1 = elem1;
	box->z2 = elem2;

	if (opcion == 0) { //No value for the box
		box->value = 0.0;
	}
	else if ((opcion == 1) || (opcion == 2)) { //The value is the area
		box->value = (box->z2[0] - box->z1[0]) * (box->z1[1] - box->z2[1]);
	}
}

bool resolver_problema_segun_parte_convexa_o_concava(int parte_a_resolver, BOILP *P, CPXLPptr lp1, CPXLPptr lp2, caja_bidimensional Box, int *indice_col,	double **nuevo_x, double **nuevo_y, double *Ideal, std::string algoritmo1, std::string algoritmo2, int *numiteraciones) {

	*nuevo_y = (double *)malloc(2 * sizeof(double));
	if (parte_a_resolver == 1) {
		if (algoritmo1 == "hybrid") {
			*nuevo_x = (double *)malloc(P->n * sizeof(double));
			return (Solve_Hibrido(Box.z1, Box.z2, P->n, P->m, P->F, *P->env, lp1, indice_col, numiteraciones, nuevo_x, nuevo_y));
		}
		else if (algoritmo1 == "tchebycheff") {
			*nuevo_x = (double *)malloc((P->n + 1) * sizeof(double));
			return (Solve_Tchebycheff(Box.z1, Box.z2, Ideal, P->n, P->m, P->F, *P->env, lp1, indice_col, numiteraciones, nuevo_x, nuevo_y));
		}
	}
	else {
		if (algoritmo2 == "hybrid") {
			*nuevo_x = (double *)malloc(P->n * sizeof(double));
			return (Solve_Hibrido(Box.z1, Box.z2, P->n, P->m, P->F, *P->env, lp2, indice_col, numiteraciones, nuevo_x, nuevo_y));
		}
		else if (algoritmo2 == "tchebycheff") {
			*nuevo_x = (double *)malloc((P->n + 1) * sizeof(double));
			return (Solve_Tchebycheff(Box.z1, Box.z2, Ideal, P->n, P->m, P->F, *P->env, lp2, indice_col, numiteraciones, nuevo_x, nuevo_y));
		}
	}
	return false;
}

void Ejecutar_algoritmo_mixed(BOILP *P1, std::string algoritmo1, std::string algoritmo2) {

	CPXENVptr env = *P1->env;				CPXLPptr lp = *P1->lp;				list<solution> *Efic = P1->Efic;
	double **F = P1->F;						int n = P1->n;						double tiempo_max = P1->max_time;		
	char sense = P1->sense;					int m = P1->m;
	
	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;
	
	bool intime;
	solution sol;
	TIEMPO t_ref;
	double hipervolumen = 0.0;
	int numiteraciones = 0;

	t_ref.init();	

	CALCULA_LEXICOGRAFICOS(P1, indice_col , &sol, t_ref);
	double S[2] = { Efic->front().vector_y[0] - 0.1 , Efic->back().vector_y[1] - 0.1 }; //Cota inferior del punto ideal

	CPXLPptr lp1, lp2;
	cargar_problemas_lp(P1, &lp1, &lp2, algoritmo1, algoritmo2, indice_col);
		
	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;

	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, 1, 1);	
	Introducir_caja_en_lista_segun_valor(&CAJAS, InitialBox); 

	double *nuevo_x = NULL, *nuevo_y = NULL;
	bool solutionfound = false;
	

	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {

		Box = CAJAS.front(); CAJAS.pop_front();

		solutionfound = resolver_problema_segun_parte_convexa_o_concava(Box.algorithm_to_use, P1, lp1, lp2, Box, indice_col, &nuevo_x, &nuevo_y, S, algoritmo1, algoritmo2, &numiteraciones);

		if (solutionfound) { //New solution found
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				hipervolumen += (Box.z1[1] - nuevo_y[1]) * (Box.z2[0] - nuevo_y[0]);
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				caja_bidimensional *newb1 = new(caja_bidimensional);
				caja_bidimensional *newb2 = new(caja_bidimensional);

				Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, 1);

				busca_siguiente_algoritmo_a_utilizar2(&Box, newb1, newb2, nuevo_y);
				
				if (newb1->z1 != NULL)	Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
				if (newb2->z1 != NULL) 	Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);

				delete(newb1);	delete(newb2);
			}
		}
		else { //No more solutions in the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		}

	}
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
	CPXfreeprob(env, &lp1);		CPXfreeprob(env, &lp2);
}

void Ejecutar_algoritmo_tchebycheff(BOILP *P1, int tipo_valor_caja) {
	// min (max{w1(z1-s1) + w2(z2-s2)} + ro*(z1-s1 + z2-s2) = min ( lambda + ro*(z1-s1+z2-s2) )
	// s.a.
	// lambda >= wi + ro * (zi - si)		i = 1,2
	// x en X
	//

	CPXENVptr env = *P1->env;					CPXLPptr lp = *P1->lp;					list<solution> *Efic = P1->Efic;
	double **F = P1->F;							int n = P1->n;							double tiempo_max = P1->max_time;
	char sense = P1->sense;						int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	int numiteraciones = 0, matbeg = 0, ind;
	bool intime;
	solution sol;
	double localnadir[2], hipervolumen = 0.0, rhs = 0.0, val;
	TIEMPO t_ref;

	t_ref.init();	

	if (!CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref)) return;
	
	double S[2] = { Efic->front().vector_y[0] - 0.1 , Efic->back().vector_y[1] - 0.1 }; //Lower bound for the ideal point (utopian point)
	
	ind = n;
	val = 1.0;

	CPXaddrows(env, lp, 1, 1, 0, &rhs, "L", &matbeg, &matbeg, &rhs, NULL, NULL); //constraint lambda >= w1(z1 - s1) Incomplete yet. Add lambda variable
	CPXaddrows(env, lp, 0, 1, 0, &rhs, "L", &matbeg, &matbeg, &rhs, NULL, NULL); //constraint lambda >= w2(z2 - s2) Incomplete yet
	
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL); //restriccion  f1 <= z21 Incomplete yet
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL); //restriccion  f2 <= z12 Incomplete yet
	
	CPXchgcoef(env, lp, m, n, -1.0); //lambda coefficient in the first constraint
	CPXchgcoef(env, lp, m + 1, n, -1.0); //lambda coefficient in the second constraint 
	CPXchgobj(env, lp, 1, &ind, &val); //lambda coefficient in the objective function

	hipervolumen = 0.0;

	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;

	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, tipo_valor_caja);	
	Introducir_caja_en_lista_segun_valor(&CAJAS, InitialBox); 

	if (P1->Tune) TuneProblem(P1); 

	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x = (double *)malloc((n + 1) * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();
		if (Solve_Tchebycheff(Box.z1, Box.z2, S, n, m, F, env, lp, indice_col, &numiteraciones, &nuevo_x, &nuevo_y)) { //New solution found
			
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				
				localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
				hipervolumen += (localnadir[1] - nuevo_y[1]) * (localnadir[0] - nuevo_y[0]);
				sol.time = t_ref.value();		sol.acum_hyper = hipervolumen;
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				Efic->push_front(sol);

				caja_bidimensional *newb1 = new(caja_bidimensional);
				caja_bidimensional *newb2 = new(caja_bidimensional);

				Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, tipo_valor_caja);
				if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
				if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);

				delete(newb1);	delete(newb2);
			}
		}
		else { //No more solutions in the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		}
	}
	CPXdelcols(env, lp, n, n); //Eliminate lambda variable
	CPXdelrows(env, lp, m, m + 1);	//Eliminate the previous two constraints
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_mixed_with_SPF(BOILP *P1, std::string algoritmo1, std::string algoritmo2) {

	CPXENVptr env = *P1->env;					CPXLPptr lp = *P1->lp;				list<solution> *Efic = P1->Efic;
	double **F = P1->F;							int n = P1->n;						double tiempo_max = P1->max_time;
	char sense = P1->sense;						int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	bool intime;
	solution sol;
	TIEMPO t_ref;
	double hipervolumen = 0.0;
	int numiteraciones = 0;
	
	t_ref.init();	

	CALCULA_LEXICOGRAFICOS(P1, indice_col , &sol, t_ref);
	double S[2] = { Efic->front().vector_y[0] - 0.1 , Efic->back().vector_y[1] - 0.1 }; //Cota inferior del punto ideal

	CPXLPptr lp1, lp2;
	cargar_problemas_lp(P1, &lp1, &lp2, algoritmo1, algoritmo2, indice_col);

	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS1, CAJAS2;

	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, 1, 1);	
	Introducir_caja_en_lista_segun_valor(&CAJAS1, InitialBox); 

	double *nuevo_x = NULL, *nuevo_y = NULL;
	bool solutionfound = false;
	
	
	intime = true;
	while ((intime) && (CAJAS1.size() > 0)) {
		
		Box = CAJAS1.front(); CAJAS1.pop_front();

		solutionfound = resolver_problema_segun_parte_convexa_o_concava(Box.algorithm_to_use, P1, lp1, lp2, Box, indice_col, &nuevo_x, &nuevo_y, S, algoritmo1, algoritmo2, &numiteraciones);

		if (solutionfound) { //New solution found
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				hipervolumen += (Box.z1[1] - nuevo_y[1]) * (Box.z2[0] - nuevo_y[0]);
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				caja_bidimensional *newb1 = new(caja_bidimensional);
				caja_bidimensional *newb2 = new(caja_bidimensional);

				Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, 1);

				double lambda1 = Box.z1[1] - Box.z2[1];
				double lambda2 = Box.z2[0] - Box.z1[0];
				double objcaja = lambda1 * Box.z1[0] + lambda2 * Box.z1[1]; 
				double objsol = lambda1 * nuevo_y[0] + lambda2 * nuevo_y[1];

				busca_siguiente_algoritmo_a_utilizar2(&Box, newb1, newb2, nuevo_y);
				
				if (objsol <= objcaja) {
					if (newb1->z1 != NULL)		Introducir_caja_en_lista_segun_valor(&CAJAS1, *newb1);
					if (newb2->z1 != NULL) 		Introducir_caja_en_lista_segun_valor(&CAJAS1, *newb2);
				}
				else {
					if (newb1->z1 != NULL) 		Introducir_caja_en_lista_segun_valor(&CAJAS2, *newb1);
					if (newb2->z1 != NULL) 		Introducir_caja_en_lista_segun_valor(&CAJAS2, *newb2);
				}
				delete(newb1);	delete(newb2);

				
			}
		}
		else { //No more solutions into the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		}
	}

	while ((intime) && (CAJAS2.size() > 0)) {

		Box = CAJAS2.front(); CAJAS2.pop_front();

		solutionfound = resolver_problema_segun_parte_convexa_o_concava(Box.algorithm_to_use, P1, lp1, lp2, Box, indice_col, &nuevo_x, &nuevo_y, S, algoritmo1, algoritmo2, &numiteraciones);

		if (solutionfound) { //New solution found
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				hipervolumen += (Box.z1[1] - nuevo_y[1]) * (Box.z2[0] - nuevo_y[0]);
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();				sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				caja_bidimensional *newb1 = new(caja_bidimensional);
				caja_bidimensional *newb2 = new(caja_bidimensional);

				Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, 1);

				busca_siguiente_algoritmo_a_utilizar2(&Box, newb1, newb2, nuevo_y);
				if (newb1->z1 != NULL)	Introducir_caja_en_lista_segun_valor(&CAJAS2, *newb1);
				if (newb2->z1 != NULL) 	Introducir_caja_en_lista_segun_valor(&CAJAS2, *newb2);
				delete(newb1);	delete(newb2);
	
				
			}
		}
		else { //No more solutions into the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		}
	}
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
	CPXfreeprob(env, &lp1);
	CPXfreeprob(env, &lp2);
}

void Crear_tres_nuevas_cajas_a_explorar_cuadrantes(caja_bidimensional *box, double *z3, double *z4,	caja_bidimensional *newb1, caja_bidimensional *newb2, caja_bidimensional *newb3, double eps1, double eps2) {

	double area1, area2, area3;
	double z11 = box->z1[0];	double z12 = box->z1[1];
	double z21 = box->z2[0];	double z22 = box->z2[1];
	double z31 = z3[0];	double z32 = z3[1];
	double z41 = z4[0];	double z42 = z4[1];

	double lambda1 = z12 - z22;	double lambda2 = z21 - z11;
	
	double z0 = lambda1 * z11 + lambda2 * z12; //
	double obj3 = lambda1 * z31 + lambda2 * z32;
	double obj4 = lambda1 * z41 + lambda2 * z42;
	
	double M = (z22 - z12) / (z21 - z11);	
	double N1 = z32 - M * z31; 
	double N2 = z42 - M * z41;

	double B, b, h; 
	double bp, hp; 

	//FIRST REGION
	B = z12 - z32;	
	if (obj3 >= z0) {
		b = 0;						h = z31 - (z12 - N1) * 1 / M;
	}
	else {
		b = z12 - (M*z11 + N1);		h = z31 - z11;
	}
	area1 = (B + b) / 2 * h;

	//SECOND REGION
	if (obj3 <= obj4) {
		B = z41 - z31;				h = z32 - z42;		hp = eps2 - z42;
		if (obj3 < obj4) {		
			bp = b = z41 - (z42 - N1) * 1 / M;		}
		else {		
			bp = b = hp = 0;	}
	}
	else {
		B = z32 - z42;			bp = b = z32 - (z31*M + N2);		h = z41 - z31;		hp = eps1 - z31;
	}
	area2 = (B + b) / 2 * h - bp * hp; 

	//THIRD REGION
	B = z21 - z41;
	if (obj4 >= z0) {
		b = 0;								h = z42 - (M * z21 + N2);
	}
	else {
		b = z21 - (z22 - N2) * 1/ M;		h = z42 - z22;
	}
	area3 = (B + b) / 2 * h;

	//CREATE NEW BOXES
	newb1->value = area1;
	newb2->value = area2;
	newb3->value = area3;

	if ((z31 - z11 == 1) || (z12 - z32 == 1)) {
		newb1->z1 = newb1->z2 = NULL;
	}
	else {
		newb1->z1 = box->z1;		newb1->z2 = z3;
	}

	if ((z41 - z31 == 1) || (z32 - z42 == 1)) {
		newb2->z1 = newb2->z2 = NULL;
	}
	else {
		newb2->z1 = z3;		newb2->z2 = z4;
	}

	if ((z21 - z41 == 1) || (z42 - z22 == 1)) {
		newb3->z1 = newb3->z2 = NULL;
	}
	else {
		newb3->z1 = z4;		newb3->z2 = box->z2;
	}

}

void Crear_dos_nuevas_cajas_a_explorar_cuadrantes(caja_bidimensional *box, double *z3, 	caja_bidimensional *newb1, caja_bidimensional *newb2, double eps1 , double eps2 , int cuadrante) {
	
	double aux, x1, y1, area1, area2;
	double z11 = box->z1[0];	double z12 = box->z1[1];
	double z21 = box->z2[0];	double z22 = box->z2[1];
	double z31 = z3[0];	double z32 = z3[1];

	double lambda1 = z12 - z22;	double lambda2 = z21 - z11;
	double z0 = lambda1 * z11 + lambda2 * z12; 
	double M = (z22 - z12) / (z21 - z11);	
	
	double zz = lambda1 * z31 + lambda2 * z32;
	double NN = z32 - M * z31;
		
	if (zz >= z0) { 
		x1 = 1.0 / M * (z12 - NN);				area1 = (z31 - x1) * (z12 - z32) / 2.0; 
		y1 = M * z21 + NN;						area2 = (z21 - z31) * (z32 - y1) / 2.0; 
		
		if (cuadrante == 2) {
			aux = 1.0 / M * (eps2 - NN);
			area2 -= (z21 - aux) * (eps2 - y1) / 2.0;
		}
		else if (cuadrante == 4) {
			aux = M * eps1 + NN;
			area1 -= (z12 - aux) * (eps1 - x1) / 2.0;
		}
	}
	else { 
		y1 = M * z11 + NN;
		x1 = 1.0 / M * (z22 - NN);
		area1 = ((z12 - z32) + (z12 - y1)) / 2.0 * (z31 - z11);
		area2 = ((z21 - z31) + (z21 - x1)) / 2.0 * (z32 - z22);

		if (cuadrante == 2) {
			aux = 1.0 / M * (eps2 - NN);
			area2 = ((z21 - z31) + (z21 - aux)) / 2.0 * (z32 - eps2);
		}
		else if (cuadrante == 4) {
			aux = M * eps1 + NN;
			area1 = ((z12 - z32) + (z12 - aux)) / 2.0 * (z31 - eps1);
		}
	}
	
	//CREATE NEW BOXES
	newb1->value = area1;
	newb2->value = area2;

	if ((z31 - z11 == 1) || (z12 - z32 == 1)) {
		newb1->z1 = newb1->z2 = NULL;
	}
	else {
		newb1->z1 = box->z1;		newb1->z2 = z3;
	}

	if ((z32 - z22 == 1) || (z21 - z31 == 1)) {
		newb2->z1 = newb2->z2 = NULL;
	}
	else {
		newb2->z1 = z3;		newb2->z2 = box->z2;
	}
}

void Ejecutar_algoritmo_econstraint_1ILP(BOILP *P1, int f_obj){
	//Execute algorithm econst1. In the end, we filter the weak efficient solutions
	//We begin by min f(f_obj)
	
	// min (f1)
	// s.a. f2(x) <= eps2

	CPXENVptr env = *P1->env;		CPXLPptr lp = *P1->lp;			list<solution> *Efic = P1->Efic;					double **F = P1->F;
	int n = P1->n;					double tiempo_max = P1->max_time;
	
	int status, f_rest = 1 - f_obj, stat, matbeg = 0, m, numiteraciones = 0 ;
	double eps, objval, localnadir[2], hipervolumen = 0.0;
	solution sol;
	bool intime;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;
	double *z1 = (double *)malloc(2 * sizeof(double));
	double *z2 = (double *)malloc(2 * sizeof(double));
	m = CPXgetnumrows(env, lp);

	//WE CALCULATE THE NADIR POINT
	if (Obtener_no_dominados_lexicograficos(P1, indice_col, &z1, &z2) != 2) return;
	localnadir[0] = z2[0];	localnadir[1] = z1[1];

	double *x_ini = (double *)malloc(n * sizeof(double)); 
	double *y_ini = (double *)malloc(2 * sizeof(double));

	status = CPXchgobj(env, lp, n, indice_col, F[f_obj]); //min f1 
	status = CPXaddrows(env, lp, 0, 1, n, &eps, "L", &matbeg, indice_col, F[f_rest], NULL, NULL); //Add constraint   f2 <= eps2.	

	TIEMPO t_ref;
	t_ref.init();

	
	if (P1->Tune) TuneProblem(P1);
	
	//Obtaining the first solution
	eps = localnadir[f_rest];
	Solve_Econstraint_1(m, eps, &env, &lp, &stat, &objval, &numiteraciones);

	if (stat == CPXMIP_INFEASIBLE) { //It is supposed that the algorithm never enters here
		free(x_ini); free(y_ini);
		return;
	}
	else {
		status = CPXgetx(env, lp, x_ini, 0, n - 1);
		y_ini[f_obj] = round(objval);
		y_ini[f_rest] = round(mult(x_ini, F[f_rest], n));
		sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = x_ini; sol.vector_y = y_ini;
		t_ref.acum();
		sol.time = t_ref.value();
		Efic->push_back(sol);
		
	}
	eps = y_ini[f_rest];

	intime = true;
	while (intime) {
		double *nuevo_x = (double *)malloc(n * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		eps--;
	
		Solve_Econstraint_1(m, eps, &env, &lp, &stat, &objval, &numiteraciones);
		
		t_ref.acum();
		if (t_ref.value() >= tiempo_max) {
			intime = false;
			free(nuevo_x); free(nuevo_y);
			numiteraciones--;
		}
		else {
			if (stat == CPXMIP_INFEASIBLE) { //Infeasible problem
				free(nuevo_x); free(nuevo_y);
				intime = false;
			}
			else {
				status = CPXgetx(env, lp, nuevo_x, 0, n - 1);
				nuevo_y[f_obj] = round(objval);
				nuevo_y[f_rest] = round(mult(nuevo_x, F[f_rest], n));

				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();
				Efic->push_back(sol);
				eps = nuevo_y[f_rest];
				
			}
		}
	}
	status = CPXdelrows(env, lp, m, m); //Restore to the initial problem

	//Filtering the dominated points and calculating the hypervolume
	std::list<solution>::iterator it1, it2;
	int num_filtradas = 0;

	it1 = it2 = Efic->begin();
	std::advance(it2, 1); //it2 points to the second element of the list

	while (it2 != Efic->end()) {
		if (round(it1->vector_y[f_obj]) == round(it2->vector_y[f_obj])) {  
			it1 = Efic->erase(it1);
			++it2;
			num_filtradas++;
		}
		else {
			hipervolumen += (localnadir[f_rest] - it1->vector_y[f_rest]) * (localnadir[f_obj] - it1->vector_y[f_obj]);
			it1->acum_hyper = hipervolumen;
			localnadir[f_rest] = it1->vector_y[f_rest];
			++it1; ++it2;
		}
	}
	P1->hypervolume = hipervolumen;
	if (it1 != Efic->end())		it1->acum_hyper = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col); free(z1); free(z2);
}

void Ejecutar_algoritmo_econstraint_2ILP(BOILP *P1, int f_obj) {
	//Execute algorithm econst2. In every iteration, we do two calls to CPLEX
	//We begin by min f(f_fob)
	
	//z = min f2
	// s.a. f1 <= eps1
	// x en X  
	//
	//min f1
	// s.a. f2 <= z
	// x en X
	

	CPXENVptr env = *P1->env;				CPXLPptr lp = *P1->lp;					list<solution> *Efic = P1->Efic;
	double **F = P1->F;						int n = P1->n;								
	double tiempo_max = P1->max_time;		char sense = P1->sense;					int m = P1->m;
	
	int status, f_rest = 1 - f_obj, matbeg = 0, j, numiteraciones = 0;
	double eps = 0.0, eps_fin, objval, hipervolumen = 0;
	solution sol;
	bool intime;
	TIEMPO t_ref;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;
	
	t_ref.init();	
	
	if (P1->Tune) TuneProblem(P1);

	//Lexicographical optimal solutions and nadir point
	double localnadir[2];
	if (!CALCULA_LEXICOGRAFICOS(P1, indice_col,  &sol, t_ref)) return;
			
 	localnadir[0] = Efic->back().vector_y[0];	localnadir[1] = Efic->front().vector_y[1];
	
	eps = localnadir[f_rest] - 1;
	eps_fin = (f_obj == 1) ? Efic->front().vector_y[0] : Efic->back().vector_y[1];

	status = CPXaddrows(env, lp, 0, 1, 1, &eps, "L", &matbeg, &matbeg, &eps, NULL, NULL); //Anadimos nueva restriccion VACIA

	status = CPXchgobj(env, lp, n, indice_col, F[f_obj]); //min f1
	for (j = 0; j < n; j++) CPXchgcoef(env, lp, m, j, F[f_rest][j]);  //s.a f2 <= eps

	status = CPXchgobj(env, lp, n, indice_col, F[f_rest]); //min f2
	for (j = 0; j < n; j++)	CPXchgcoef(env, lp, m+1, j, F[f_obj][j]);  //Modify (m+1)-constraint in order to create constraint f1 <= f2(x*)
	

	intime = true;
	while (intime) {
		double *nuevo_x = (double *)malloc(n * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));
				
		Solve_Econstraint_22(m, eps, env, lp, indice_col, n, f_obj, F, &objval, &numiteraciones, &nuevo_y);
	
		t_ref.acum();

		if (t_ref.value() >= tiempo_max) {
			free(nuevo_x); free(nuevo_y);
			numiteraciones -= 2;
			intime = false;
		}
		else {
			eps = round(objval - 1);
			if (eps < eps_fin) {
				free(nuevo_x); free(nuevo_y);
				intime = false;
			}
			else {
				CPXgetx(env, lp, nuevo_x, 0, n - 1);
	
				hipervolumen += (localnadir[f_rest] - nuevo_y[f_rest]) * (localnadir[f_obj] - nuevo_y[f_obj]);
				localnadir[f_rest] = nuevo_y[f_rest];

				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
				Efic->push_back(sol);
				
			}
		}
	}
	CPXdelrows(env, lp, m, m); //Restore to the original problem
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_SPF_exact(BOILP *P1) {
	//Supported Pareto front

	// min (lambda1 * f1(x) + lambda2 * f2)   
	// s.t.  x en X
	//	f1 <= k1
	//  f2 <= k2
	// If the solution is not supported, it is discarded

	CPXENVptr env = *P1->env;				CPXLPptr lp = *P1->lp;				list<solution> *Efic = P1->Efic;
	double **F = P1->F;						int n = P1->n;						double tiempo_max = P1->max_time;		
	char sense = P1->sense;					int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	bool intime;
	solution sol;
	double localnadir[2], hipervolumen = 0.0;
	int matbeg = 0, numiteraciones = 0;
	double rhs = 0.0;
	
	TIEMPO t_ref;

	t_ref.init();	

	if (!CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref)) return;
	
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL);
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL);

	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;

	InitialBox.z1 = Efic->front().vector_y;	InitialBox.z2 = Efic->back().vector_y;
	CAJAS.push_back(InitialBox);

	if (P1->Tune) TuneProblem(P1);

	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x = (double *)malloc(n * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();

		if (Solve_SPF(Box.z1, Box.z2, n, m, F, env, lp, indice_col, &numiteraciones, &nuevo_x, &nuevo_y)) { //New solution found
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				caja_bidimensional newb;
				
				localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
				hipervolumen += (localnadir[1] - nuevo_y[1]) * (localnadir[0] - nuevo_y[0]);
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();		sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				newb.z1 = Box.z1;	newb.z2 = nuevo_y;				CAJAS.push_back(newb);
				newb.z1 = nuevo_y; newb.z2 = Box.z2;				CAJAS.push_back(newb);
				
			}
		}
		else { //No more supported solutions in the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		};
	}
	CPXdelcols(env, lp, m, m + 1);
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Hallar_cuadrante_solucion_encontrada(double *z1, double *z2, double *eps1 , double *eps2, double *z3 , int *cuadrante) {
	*eps1 = int ((z1[0] + z2[0]) / 2);
	*eps2 = int((z1[1] + z2[1]) / 2);

	if ((z3[0] >= *eps1) && (z3[1] >= *eps2)) *cuadrante = 1;
	else if ((z3[0] < *eps1) && (z3[1] >= *eps2)) *cuadrante = 2;
	else if ((z3[0] < *eps1) && (z3[1] < *eps2)) *cuadrante = 3;
	else *cuadrante = 4;
}

void Ejecutar_algoritmo_hibrido_por_cuadrantes(BOILP *P1) {

	//THIS OPTION IS NOT USED IN THE EXPERIMENTS
	//Executing hypbrid algorithm by quadrants. 

	// min (lambda1 * f1(x) + lambda2 * f2)   
	// s.t.  x en X
	//	f1 <= k1
	//  f2 <= k2


	CPXENVptr env = *P1->env;					CPXLPptr lp = *P1->lp;						list<solution> *Efic = P1->Efic;
	double **F = P1->F;							int n = P1->n;								double tiempo_max = P1->max_time;
	char sense = P1->sense;						int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;
	
	bool intime;
	solution sol;
	TIEMPO t_ref;
	int cuadrante, numiteraciones = 0, matbeg = 0;
	double eps1, eps2 ,localnadir[2], hipervolumen = 0.0, rhs = 0.0;


	t_ref.init();	

	CALCULA_LEXICOGRAFICOS(P1, indice_col , &sol, t_ref);

	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL); 
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL); 

	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;

	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, 2);	
	Introducir_caja_en_lista_segun_valor(&CAJAS, InitialBox); 

	
	if (P1->Tune) TuneProblem(P1);
	
	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x3 = (double *)malloc(n * sizeof(double));
		double *nuevo_y3 = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();
				
		if (Solve_Hibrido(Box.z1, Box.z2, n, m, F, env, lp, indice_col, &numiteraciones, &nuevo_x3, &nuevo_y3)) { 

			Hallar_cuadrante_solucion_encontrada(Box.z1, Box.z2, &eps1, &eps2, nuevo_y3, &cuadrante);

			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x3); free(nuevo_y3);
			}
			else {
				sol.dimension_x = n + 1; sol.dimension_y = 2; sol.vector_x = nuevo_x3; sol.vector_y = nuevo_y3;
				
				if ((cuadrante == 1) || (cuadrante == 3)) {
					caja_bidimensional *newb1 = new(caja_bidimensional);
					caja_bidimensional *newb2 = new(caja_bidimensional);
					Crear_dos_nuevas_cajas_a_explorar_cuadrantes(&Box, nuevo_y3, newb1, newb2, eps1, eps2, cuadrante);
					if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
					if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);
					delete(newb1);	delete(newb2);
					localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
					hipervolumen += (localnadir[1] - nuevo_y3[1]) * (localnadir[0] - nuevo_y3[0]);
					
					sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
					Efic->push_back(sol);
				}
				else if (cuadrante == 2) {
					double *nuevo_x4 = (double *)malloc(n * sizeof(double));
					double *nuevo_y4 = (double *)malloc(2 * sizeof(double));
					double EPS[2] = { Box.z2[0] , eps2 };
					if (Solve_Hibrido_cuadrantes(EPS[0], EPS[1], n, m, F, env, lp, &numiteraciones, &nuevo_x4, &nuevo_y4)) { //Encontrada nueva solucion											
						t_ref.acum();
						if (t_ref.value() >= tiempo_max) {
							intime = false;
							free(nuevo_x4); free(nuevo_y4);
						}
						else {
							caja_bidimensional *newb1 = new(caja_bidimensional);
							caja_bidimensional *newb2 = new(caja_bidimensional);
							caja_bidimensional *newb3 = new(caja_bidimensional);
							Crear_tres_nuevas_cajas_a_explorar_cuadrantes(&Box, nuevo_y3, nuevo_y4, newb1, newb2, newb3, eps1, eps2);
							if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
							if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);
							if (newb3->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb3);
							delete(newb1);	delete(newb2);
							localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
							hipervolumen += (localnadir[1] - nuevo_y3[1]) * (localnadir[0] - nuevo_y3[0]);
							sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
							Efic->push_back(sol);
							localnadir[0] = Box.z2[0];			localnadir[1] = nuevo_y3[1];
							hipervolumen += (localnadir[1] - nuevo_y4[1]) * (localnadir[0] - nuevo_y4[0]);
							sol.dimension_x = n + 1; sol.dimension_y = 2; sol.vector_x = nuevo_x4; sol.vector_y = nuevo_y4;
							sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
							Efic->push_back(sol);
						}
					}
					else { 
						caja_bidimensional *newb1 = new(caja_bidimensional);
						caja_bidimensional *newb2 = new(caja_bidimensional);
						Crear_dos_nuevas_cajas_a_explorar_cuadrantes(&Box, nuevo_y3, newb1, newb2, eps1, eps2, cuadrante);
						if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
						if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);
						delete(newb1);	delete(newb2);
						localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
						hipervolumen += (localnadir[1] - nuevo_y3[1]) * (localnadir[0] - nuevo_y3[0]);
					}
				}
				else if (cuadrante == 4) {
					double *nuevo_x4 = (double *)malloc(n * sizeof(double));
					double *nuevo_y4 = (double *)malloc(2 * sizeof(double));
					double EPS[2] = { eps1 , Box.z1[1] };
					if (Solve_Hibrido_cuadrantes(EPS[0], EPS[1], n, m, F, env, lp, &numiteraciones, &nuevo_x4, &nuevo_y4)) {
																					
						t_ref.acum();
						if (t_ref.value() >= tiempo_max) {
							intime = false;
							free(nuevo_x4); free(nuevo_y4);
						}
						else {
							caja_bidimensional *newb1 = new(caja_bidimensional);
							caja_bidimensional *newb2 = new(caja_bidimensional);
							caja_bidimensional *newb3 = new(caja_bidimensional);
							Crear_tres_nuevas_cajas_a_explorar_cuadrantes(&Box, nuevo_y4, nuevo_y3, newb1, newb2, newb3, eps1, eps2);
							if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
							if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);
							if (newb3->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb3);
							delete(newb1);	delete(newb2);
							localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
							hipervolumen += (localnadir[1] - nuevo_y4[1]) * (localnadir[0] - nuevo_y4[0]);
							sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
							Efic->push_back(sol);
							localnadir[0] = Box.z2[0];			localnadir[1] = nuevo_y4[1];
							hipervolumen += (localnadir[1] - nuevo_y3[1]) * (localnadir[0] - nuevo_y3[0]);
							sol.dimension_x = n + 1; sol.dimension_y = 2; sol.vector_x = nuevo_x4; sol.vector_y = nuevo_y4;
							sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
							Efic->push_back(sol);
							
						}
					}
					else { 
						caja_bidimensional *newb1 = new(caja_bidimensional);
						caja_bidimensional *newb2 = new(caja_bidimensional);
						Crear_dos_nuevas_cajas_a_explorar_cuadrantes(&Box, nuevo_y3, newb1, newb2, eps1, eps2, cuadrante);
						if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
						if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);
						delete(newb1);	delete(newb2);
						localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
						hipervolumen += (localnadir[1] - nuevo_y3[1]) * (localnadir[0] - nuevo_y3[0]);
						
					}
				}
			}
		}
		else { 
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x3); free(nuevo_y3);
		}
	}
	CPXdelrows(env, lp, m, m + 1); 
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_hibrido(BOILP *P1, int opcion_caja) {
	//Execute hybrid algorithm
	// min (lambda1 * f1(x) + lambda2 * f2)   
	// s.t.  x en X
	//	f1 <= k1
	//  f2 <= k2
	
	CPXENVptr env = *P1->env;			CPXLPptr lp = *P1->lp;				list<solution> *Efic = P1->Efic;
	double **F = P1->F;					int n = P1->n;						double tiempo_max = P1->max_time;	
	char sense = P1->sense;				int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	bool intime;
	solution sol;
	double localnadir[2], hipervolumen = 0.0, rhs = 0.0;
	int matbeg = 0, numiteraciones = 0;
	TIEMPO t_ref;

	t_ref.init();	

	if (!CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref)) return;
		
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[0], NULL, NULL);
	CPXaddrows(env, lp, 0, 1, n, &rhs, "L", &matbeg, indice_col, F[1], NULL, NULL);
	
	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;
	
	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, opcion_caja);	
	Introducir_caja_en_lista_segun_valor(&CAJAS, InitialBox); 

	
	if (P1->Tune) TuneProblem(P1);
	
	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x = (double *)malloc(n * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();

		if (Solve_Hibrido(Box.z1, Box.z2, n, m, F, env, lp, indice_col, &numiteraciones, &nuevo_x, &nuevo_y)) { //New solution found

			t_ref.acum();
			if (t_ref.value() >= tiempo_max) {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
			else {
				solution sol;
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
				hipervolumen += (localnadir[1] - nuevo_y[1]) * (localnadir[0] - nuevo_y[0]);

				caja_bidimensional *newb1 = new(caja_bidimensional);
				caja_bidimensional *newb2 = new(caja_bidimensional);

				Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, opcion_caja);
				if (newb1->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS , *newb1);
				if (newb2->z1 != NULL) Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);

				delete(newb1);	delete(newb2);
			
			}
		}
		else { //No more solutions inside the box
			t_ref.acum();
			if (t_ref.value() >= tiempo_max) intime = false;
			free(nuevo_x); free(nuevo_y);
		};
	}
	CPXdelrows(env, lp, m, m + 1); //Restore
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_SPF_approximate(BOILP *P1) {
	//Supported Pareto front
	// min (lambda1 * f1(x) + lambda2 * f2)   
	// s.t.  x en X
	
	CPXENVptr env = *P1->env;						CPXLPptr lp = *P1->lp;					list<solution> *Efic = P1->Efic;
	double **F = P1->F;								int n = P1->n;							double tiempo_max = P1->max_time;	
	char sense = P1->sense;					int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	
	int matbeg = 0, numiteraciones = 0;
	bool intime;
	solution sol;
	double localnadir[2], hipervolumen = 0.0;
	TIEMPO t_ref;

	t_ref.init();	

	//Lexicographical optimal solutions
	CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref);
	
	hipervolumen = 0.0;
	caja_bidimensional InitialBox, Box;
	InitialBox.z1 = Efic->front().vector_y;
	InitialBox.z2 = Efic->back().vector_y;

	list<caja_bidimensional> CAJAS;
	CAJAS.push_back(InitialBox);

	if (P1->Tune) TuneProblem(P1);

	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x = (double *)malloc(n * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();

		if (Solve_SPF_approximate(Box.z1, Box.z2, n, F, env, lp, indice_col, &numiteraciones, &nuevo_x, &nuevo_y)) { //Encontrada nueva solucion
			t_ref.acum();
			if (t_ref.value() < tiempo_max) {
				solution sol;
				caja_bidimensional newb;
				
				localnadir[0] = Box.z2[0];			localnadir[1] = Box.z1[1];
				hipervolumen += (localnadir[1] - nuevo_y[1]) * (localnadir[0] - nuevo_y[0]);
				sol.dimension_x = n; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
				sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
				Efic->push_front(sol);

				newb.z1 = Box.z1;	newb.z2 = nuevo_y;				CAJAS.push_front(newb);
				newb.z1 = nuevo_y; newb.z2 = Box.z2;				CAJAS.push_front(newb);
				
			}
			else {
				intime = false;
				free(nuevo_x); free(nuevo_y);
			}
		}
		else { //The solution was previously found
			free(nuevo_x); free(nuevo_y);
		}
	}
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_augmecon_anytime(BOILP *P1,int f_obj, double lambda_) {
	//Algorithm Augmented e-constraint. Boxes with greater area are explored first
	// We begin by min f(f_obj).
	// Lambda might be fixed (between  [10^-6 , 10^-3]) or variable
	// min (f1(x) - lambda * s2)   
	// s.t.  f2(x) + s2 <= eps2
	//			s2 >= 0
	//			x en X
	// (we can change f1 <--> f2)
	
	CPXENVptr env = *P1->env;						CPXLPptr lp = *P1->lp;							list<solution> *Efic = P1->Efic;
	double **F = P1->F;								int n = P1->n;									double tiempo_max = P1->max_time;
	char sense = P1->sense;							int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;
	
	int f_rest = 1 - f_obj, matbeg = 0, numiteraciones = 0;
	double eps = 0.0, eps_fin, objval, lambda = 0.0, hipervolumen = 0.0;
	solution sol;
	bool intime;
	TIEMPO t_ref;
	
	t_ref.init();	

	//Lexicographical optimal solutions
	CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref);

	caja_bidimensional InitialBox, Box;
	list<caja_bidimensional> CAJAS;

	Inicializar_caja_a_explorar(&InitialBox, Efic->front().vector_y, Efic->back().vector_y, 1);
	Introducir_caja_en_lista_segun_valor(&CAJAS, InitialBox);

	eps_fin = (f_obj == 1) ? Efic->front().vector_y[0] : Efic->back().vector_y[1];

	CPXchgobj(env, lp, n, indice_col, F[f_obj]);	
	CPXaddrows(env, lp, 1, 1, n, &eps, "L", &matbeg, indice_col, F[f_rest], NULL, NULL); //Add variable "s" and create constraint   f1 <= eps   <-- incomplete
	CPXchgcoef(env, lp, m, n, 1.0); //f1 + s <= eps  (Modify constraint)

	bool lambdaisfixed = (lambda_ == 0) ? false : true;

	if (P1->Tune) TuneProblem(P1);

	intime = true;
	while ((intime) && (CAJAS.size() > 0)) {
		double *nuevo_x = (double *)malloc((n + 1) * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		Box = CAJAS.front(); CAJAS.pop_front();
		eps = (Box.z1[f_rest] + Box.z2[f_rest]) / 2;

		if (lambdaisfixed) lambda = lambda_;
		else {
			if (f_obj == 1)  
				lambda = 1.0 / (Box.z2[f_rest] - eps_fin + 0.001);
			else			
				lambda = 1.0 / (Box.z1[f_rest] - eps_fin + 0.001);
		}
		

		Solve_Augmecon(F[f_obj], F[f_rest], lambda, m, eps, env, lp, indice_col, n, &objval, &numiteraciones);

		CPXgetx(env, lp, nuevo_x, 0, n);
		nuevo_y[f_obj] = round(mult(F[f_obj], nuevo_x, n));  
		nuevo_y[f_rest] = round(mult(F[f_rest], nuevo_x, n));

		t_ref.acum();

		if ((t_ref.value() < tiempo_max) && (nuevo_y[0] <= Box.z1[0]) || (nuevo_y[1] <= Box.z2[1])) { //Solution previously found
			caja_bidimensional *newb3 = new(caja_bidimensional);
			double *aa = (double *)malloc(2 * sizeof(double));

			if (f_obj == 0) {		aa[0] = Box.z2[0];			aa[1] = eps;				Inicializar_caja_a_explorar(newb3, Box.z1, aa, 1);			}
			else			{		aa[0] = eps;				aa[1] = Box.z1[1];			Inicializar_caja_a_explorar(newb3, aa, Box.z2, 1);			}

			free(nuevo_x); free(nuevo_y);
			Introducir_caja_en_lista_segun_valor(&CAJAS, *newb3);
		}
		else if (t_ref.value() >= tiempo_max) {
			free(nuevo_x); free(nuevo_y);	numiteraciones--;
			intime = false;
		}
		else {
			hipervolumen += (Box.z2[0] - nuevo_y[0]) * (Box.z1[1] - nuevo_y[1]);
			sol.dimension_x = n + 1; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y;
			sol.time = t_ref.value();		sol.acum_hyper = hipervolumen;
			Efic->push_back(sol);

			caja_bidimensional *newb1 = new(caja_bidimensional);
			caja_bidimensional *newb2 = new(caja_bidimensional);

			Crear_dos_nuevas_cajas_a_explorar(&Box, nuevo_y, newb1, newb2, 1);
			if (newb1->z1 != NULL)			Introducir_caja_en_lista_segun_valor(&CAJAS, *newb1);
			if (newb2->z1 != NULL)			Introducir_caja_en_lista_segun_valor(&CAJAS, *newb2);

			delete(newb1);	delete(newb2);
		}
	}
	CPXdelcols(env, lp, n, n); //Eliminate variable s
	CPXdelrows(env, lp, m, m); //Eliminate constraint  f2(x) + s2 <= eps2
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}

void Ejecutar_algoritmo_augmecon_exact(BOILP *P1, int f_obj, double lambda_) {
	//Algorithm Augmented E-constraint.
	// We begin by min f(f_obj)
	// The lambda value might be fixed (between [10^-6 , 10^-3]) or variable
	// min (f1(x) - lambda * s2)   
	// s.t.  f2(x) + s2 <= eps2
	//			s2 >= 0
	//			x en X
	// (it can be f1 <--> f2)

	CPXENVptr env = *P1->env;			CPXLPptr lp = *P1->lp;			list<solution> *Efic = P1->Efic;
	double **F = P1->F;					int n = P1->n;					double tiempo_max = P1->max_time;		
	char sense = P1->sense;				int m = P1->m;

	int *indice_col = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) indice_col[i] = i;

	int f_rest = 1 - f_obj, matbeg = 0, numiteraciones = 0;
	double eps = 0.0, eps_fin, objval, lambda = 0.0, hipervolumen = 0.0, localnadir[2];
	solution sol;
	bool intime;
	TIEMPO t_ref;

	t_ref.init();	

	if (!CALCULA_LEXICOGRAFICOS(P1, indice_col, &sol, t_ref)) return;

	localnadir[0] = Efic->back().vector_y[0];		localnadir[1] = Efic->front().vector_y[1];

	eps = localnadir[f_rest] - 1;
	eps_fin = (f_obj == 1) ? Efic->front().vector_y[0] : Efic->back().vector_y[1];

	CPXchgobj(env, lp, n, indice_col, F[f_obj]);	//min f2  <--incomplete
	CPXaddrows(env, lp, 1, 1, n, &eps, "L", &matbeg, indice_col, F[f_rest], NULL, NULL); //Add variable "s" and create constraint  f1 <= eps   <-- incomplete
	CPXchgcoef(env, lp, m, n, 1.0); //f1 + s <= eps  (Modify constraint)

	bool lambdaisfixed = (lambda_ == 0) ? false : true;
		
	if (P1->Tune) TuneProblem(P1);

	intime = true;
	while (intime) {
		double *nuevo_x = (double *)malloc((n + 1) * sizeof(double));
		double *nuevo_y = (double *)malloc(2 * sizeof(double));

		lambda = (lambdaisfixed) ? lambda_ : 1.0 / (eps + 1 - eps_fin + 0.001); //fix or variable lambda

		Solve_Augmecon(F[f_obj], F[f_rest], lambda, m, eps, env, lp, indice_col, n, &objval, &numiteraciones);

		CPXgetx(env, lp, nuevo_x, 0, n);
		
		nuevo_y[f_obj] = round(objval + lambda * nuevo_x[n]);	
		nuevo_y[f_rest] = round(mult(F[f_rest], nuevo_x, n));

		eps = nuevo_y[f_rest] - 1;
		
		t_ref.acum();
		if ((t_ref.value() >= tiempo_max) || (eps < eps_fin)) {
			if (t_ref.value() >= tiempo_max) numiteraciones--;
			free(nuevo_x); free(nuevo_y);
			intime = false;
		}
		else {
			hipervolumen += (localnadir[f_rest] - nuevo_y[f_rest]) * (localnadir[f_obj] - nuevo_y[f_obj]);
			localnadir[f_rest] = nuevo_y[f_rest];

			sol.dimension_x = n + 1; sol.dimension_y = 2; sol.vector_x = nuevo_x; sol.vector_y = nuevo_y; 
			sol.time = t_ref.value();	sol.acum_hyper = hipervolumen;
			Efic->push_back(sol);

			if (eps == eps_fin) 
				intime = false;
		}
	

	}
	CPXdelcols(env, lp, n, n); //Eliminate variable s
	CPXdelrows(env, lp, m, m); //Eliminate constraint f2(x) + s2 <= eps2
	P1->hypervolume = hipervolumen;
	P1->n_iterations += numiteraciones;
	free(indice_col);
}
