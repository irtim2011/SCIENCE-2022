#include <iostream>
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <cassert>
#define pi 3.1415926

using namespace std; //1112

typedef vector <double> vec_t;
vector <double> progonka(vec_t a,vec_t b, vec_t c, vec_t d, int n) {
	assert(a.size() == n);
	assert(b.size() == n);
	assert(c.size() == n);
	assert(d.size() == n);
	vector<double> C(n);
	vector<double> D(n + 1); // tut nol ne usaem
	a.insert(a.begin(), 1, 0);
	b.insert(b.begin(), 1, 0);
	c.insert(c.begin(), 1, 0);
	d.insert(d.begin(), 1, 0);
	C[1] = c[1] / b[1];
	for (int i = 2; i <= n - 1; i++) { 
		C[i] = c[i] / (b[i] - a[i] * C[i - 1]);
	}
	D[1] = d[1] / b[1];
	for (int i = 2; i <= n; i++) { 
		D[i] = (d[i] - a[i] * D[i - 1]) / (b[i] - a[i] * C[i - 1]);
	}
	vector < double> result(n + 2);
	result[n] = D[n];

	for (int i = n - 1; i >= 1; i--) {
		result[i] = D[i] - C[i] * result[i + 1];
	
	}
	vec_t result_1(n);
	for (int i = 0; i < n; i++){
		result_1[i] = result[i + 1];
	}

	return result_1;
}

bool check_progonka(vec_t a, vec_t b, vec_t c, vec_t d, vec_t x, int n) {
	double constexpr eps = 1e-7;
	for (int i = 0; i < n; i++) {
		double res = x[i] * b[i];
		if (i > 0) res += a[i] * x[i - 1];
		if (i < (n - 1)) res += c[i] * x[i + 1];
		if (abs(res - d[i]) > eps) 
		{ return false; }
	}
	return true;
}

void test_progonka() {
	vec_t a = { 0,1,1 };
	vec_t b = { 1,2,1 };
	vec_t c = {0, 1, 0};
	vec_t d = { 1,8,5 };
	vec_t x = progonka(a, b, c, d, 3);
	assert(x[0] == 1);
	assert(x[1] == 2);
	assert(x[2] == 3);
}

double distance_L2(vector<double> a, vector<double> b) {
	assert(a.size() == b.size());
	double result = 0;
	for (int i = 0; i < a.size(); i++) {
		result += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(result);
}




vector <double> init_q(double U_0, double H, vector <double> eta) {
	const int N = eta.size() - 1;
	vector <double> F(eta.size());
	vector <double> q(N);
	for (int j = 0; j < N + 1; j++) {
		F[j] = 3 * U_0 / 2 * (eta[j] * eta[j] / H - eta[j] * eta[j] * eta[j] / 3 / H / H);
	}
	for (int j = 0; j < N; j++) {
		q[j] = F[j + 1] - F[j];
	}
	return q;
}

vector <double> init_q_alternative(double xi, vector <double> u, vector <double> h){
vector <double> q(h.size() - 1);
double sum_q = 0;
for (int i = 0; i < h.size() - 1; i++) {
	q[i] = (xi * (u[i + 1] + u[i]) * (h[i + 1] - h[i])) / 2.;
	sum_q += q[i];
}
return q;
}
//L -- u_x = L(u) * u_yy
double L(double u) {
	return 1/u;
}

double left_(double eta, double H, double U_0) {
	return (3 * U_0 / 2) * (2 * eta / H - (eta / H) * (eta / H));
}

// uk -- prevedushiy slice 
// uk_1_prev --sledyshiy slice s prevedushei iterazii
// u_bottom, u_top --kraevie znacheniya
//h_setka -- сетка на следующем слое полученная на преведущей итерации(т.е h_i+1 ^ old)
//result --uk+1 -new 

//C_x = -dh/dx + alpha*x  (это константа для каждого слоя по иксу)
//
vector<double> recompute(int N_y, double dx, double C_x, vector <double> uk, vector <double> uk_1_prev, double u_bottom, vector <double> h_setka) {
	vector<double> a(N_y + 2), b(N_y + 2), c(N_y + 2), d(N_y + 2);
	for (int i = 1; i < N_y + 1; i++) {
		double dyi = h_setka[i] - h_setka[i - 1];
		double dyi1 = h_setka[i + 1] - h_setka[i];
		double dd = dyi + dyi1;
		double A = 1 / dyi / dd;
		double B = 1 / dyi1 / dd;
		double A_old = A * L(uk[i]);	
		double B_old = B * L(uk[i]);
		double A_new = A * L(uk_1_prev[i]);
		double B_new = B * L(uk_1_prev[i]);
		a[i] = -A_new;
		b[i] = (1 / dx + A_new + B_new);
		c[i] = -B_new;
		d[i] = uk[i] / dx + A_old * (uk[i - 1] - uk[i]) + B_old * (uk[i + 1] - uk[i]) - C_x; //если не реботает прога первое что надо делать меня знак C_x
		//if (i == N_y) { d[i] += B_new * u_top; }
	}

	a[N_y + 1] = 1; //если не работает то проверить эти граничные условие !!!
	b[N_y + 1] = -1;//старая версия
	c[N_y + 1] = 0;
	d[N_y + 1] = 0;


	a[0] = 0;
	b[0] = 1;
	c[0] = 0;
	d[0] = 0;

	//d[1] = d[1] + 0.5 * 1 / (dy) * 1 / (dy)*u[k + 1][0];
	//d[N_y] = d[N_y] + 0.5 * 1 / (dy) * 1 / (dy)*u[k + 1][N_y + 1];



	//vector<double> C(N_y);
	//vector<double> D(N_y + 1);
	//C[1] = c[1] / b[1];
	//for (int i = 2; i <= N_y - 1; i++) { //why N_y
	//	C[i] = c[i] / (b[i] - a[i] * C[i - 1]);
	//}
	//D[1] = d[1] / b[1];
	//for (int i = 2; i <= N_y; i++) { //N_y + 1
	//	D[i] = (d[i] - a[i] * D[i - 1]) / (b[i] - a[i] * C[i - 1]);
	//}
	//vector < double> result(N_y + 2);
	//result[N_y] = D[N_y]; //N_y + 1

	//for (int i = N_y - 1; i >= 1; i--) {//N_y...1
	//	result[i] = D[i] - C[i] * result[i + 1];
	//}
	////result[N_y + 1] = result[N_y]; //dirty hack!!!!!!!!!!!!
	vector <double> result = progonka(a, b, c, d, N_y + 2);
	bool variable_check_progonka = check_progonka(a, b, c, d, result, N_y + 2);
	if (!variable_check_progonka) { cout << "ERROR PROGONKA" << endl; exit(-1);}

	for (int i = 1; i < N_y + 1; i++) {//Zatichka Jutkaia
		if (result[i] < 0)
		{
			result[i] = 0;
		}
	}

	return result;
}
//x_min *(u[0][4]*h[0][4] + u[0][3]*h[0][3])*(h[0][4] - h[0][3])
// x_i = ksi (на каждом i-ом слое)
vector <double> h_recompute(double x_i, vector <double> u_new,vector <double> q) {
	int N = q.size();
	double q_sum = 0;
	vector <double> h_new(N + 1);// to do
	h_new[0] = 0;
	for (int n = 0; n < N; n++) {
		//double a = x_i * u_new[n + 1];
		//double b = -x_i * h_new[n] * (u_new[n + 1] - u_new[n]);
		//double c = -x_i * u_new[n] * h_new[n] * h_new[n] - 2 * q[n];
		//double D = b * b - 4 * a * c;
		//if (a < 0) { h_new[n + 1] = -c / b;}
		//else { h_new[n + 1] = (-b + sqrt(D)) / 2 / a; }
		q_sum += q[n];
		h_new[n + 1] = 2. * q[n] / (x_i * (u_new[n + 1] + u_new[n])) + h_new[n];
	}
	return h_new;
 }

bool check_invariant(double x_i, vector <double> u_new, vector <double> q, vector <double> h_new) {
	int n = q.size();
	for (int i = 0; i  < n; i++) {
		double tmp = x_i * (u_new[i + 1] + u_new[i]) * (h_new[i + 1] - h_new[i]) / 2.;
		if (abs(tmp - q[i]) > 1.0e-6) return false;
	}
	return true;
}

vector <double> h_recompute_alternative(double x_i, vector <double> u_old, vector <double> q) {
	int N = q.size();
	vector <double> h_new(N + 1);// to do
	h_new[0] = 0;	
	for (int n = 0; n < N; n++) {
		double a = x_i * u_old[n + 1];
		double b = -x_i * h_new[n] * (u_old[n + 1] - u_old[n]);
		double c = -x_i * u_old[n] * h_new[n] * h_new[n] - 2 * q[n];
		double D = b * b - 4 * a * c;
		assert(D >= 0);
		assert(a != 0);
		//if (a < 0) { h_new[n + 1] = -c / b; }
		h_new[n + 1] = (-b + sqrt(D)) / 2 / a; 
	}
	return h_new;
}

void Solve_equation() {
	const double x_min = 0.1;//last dot in x_setka could not equal(less) than x_max
	const double x_max = 1;
	const double dx = 0.005;
	const int N_y = 200; //На один меньше чем число интервалов и на два меньше чем число точек
	//число внутренних точек разбиения не учитывая концы //
	const double y_min = 0;

	int N_x = round((x_max - x_min) / dx - 1); //На один меньше чем число интервалов и на два меньше чем число точек
	//число внутренних точек разбиения не учитывая концы

	vector<double> x_setka(N_x + 2);

	for (int i = 0; i < N_x + 2; i++) {
		x_setka[i] = x_min + dx * i;
	}

	const double U_0 = 10; //Average speed at x_0
	const double H = 1 / U_0 / x_min;

	const double y_max = H;

	double  dy = H / (N_y + 1);

	vector <double> eta(N_y + 2);
	for (int i = 0; i < N_y + 2; i++) {
		//y_setka[i] = (1 / 2.) - 1 / 2. * cos(pi * i / (N_y + 1));
		eta[i] = y_min + dy * i;
	}

	ifstream F("f_values");
	vector <double> f_values(N_y + 2);
	for (int i = 0; i < N_y + 2; i++) {
		F >> f_values[i];
	}

	vector<vector<double> > u(N_x + 2, vector<double>(N_y + 2));

	for (int i = 0; i < N_x + 2; i++) {
		u[i][0] = 0;
	}

	for (int i = 1; i < N_y + 2; i++) {
		u[0][i] = U_0 * f_values[i];
	}

	const vector <double> q = init_q_alternative(x_min, u[0], eta);

	vector <double> u_old(N_y + 2);
	vector <double> u_new(N_y + 2);

	vector< vector <double> > h(N_x + 2, vector<double>(N_y + 2));
	//double D = -1 / U_0 + 1.814 * x_min * x_min * x_min / 3;
	
	
	//double D = -1 / U_0 - 1.814 * x_min * x_min * x_min / 3;
	//double h_0 = -1.814 * x_min * x_min / 3 - D / x_min;

	double c = -1.814;
	double D = -1 / U_0 - c * x_min * x_min * x_min / 3;
	double h_0 =  -c * x_min * x_min / 3 - D / x_min;

	cout << "h_0 correct=" << 1 / U_0 / x_min;

	cout << "D= " << D << " h_0= " << h_0 << endl;
	
	for (int i = 0; i < N_y + 2; i++) {
		h[0][i] = eta[i] * h_0;
	}
	
	//h[0] = eta;
	vector <double> h_old(N_y + 2);
	vector <double> h_new(N_y + 2);

	const double alpha = 0.2; // Uznat

	for (int k = 0; k < N_x + 1; k++) {// v etom meste u[k], h[k](...) uzhe poschitani
	//for (int k = 0; k < N_x; k++) { //!!!ЭТО ВРЕМЕННАЯ ПРОВЕРКА БЕЗ ПОСЛЕДНЕГО ШАГА ТК РАЗВАЛИВАЛОСЬ!!!
		double discrepancy = 1;
		u_old = u[k];
		u_old[N_y + 1] = u_old[N_y]; //DIRTY HACK
		//h_old = h[k];
		for (int i = 0; i < N_y + 2; i++) {
			h_old[i] = (1.814 * x_setka[k + 1] * x_setka[k + 1] / 3 - D / x_setka[k + 1]) * eta[i];
			if (h_old[i] < 0) { cout << "ERROR line 305" << endl;
			exit(-2); }
		}
		int iter_whl = 0;
		while (discrepancy > 0.00000001) {
			//-dh / dx + alpha * x
			double dh = h_old[N_y + 1] - h[k][N_y + 1];
			//double C_x = -dh / dx + alpha * x_setka[k];//первые два слагаемых из правой части уравнения в красной рамке
			double C_x = 0; //решение Вотсона
			u_new = recompute(N_y, dx, C_x, u[k], u_old, 0, h_old);//to do u_new[N_y + 1] не вычисляется 
			u_new[0] = u[k + 1][0];

			//cout << "pon" << endl;
			for (int i = 1; i < N_y + 1; i++) {//Zatichka Jutkaia
				if (u_new[i] < 0) 
				{
					cout << "zatichka" << endl;
					u_new[i] = 9999;}
			}
			//u_new[N_y + 1] = u[k + 1][N_y + 1]; //to do  zasynyt v recompute i razobratsya
			h_new = h_recompute(x_setka[k + 1], u_old, q);
			//h_new = h_recompute(x_setka[k], u_new, q); //there are doubts about what recompute first u_new or h_new
			//h_new = h_recompute_alternative(x_setka[k], u_old, q); // funcan <u_old, h_old> --> <u_new, h_new> 
																   // Почему оператор сжимающий? 
			double dist1 = distance_L2(u_old, u_new);
			double dist2 = distance_L2(h_old, h_new);
			//discrepancy = distance_L2(u_old, u_new);
			discrepancy = dist1 + dist2;

			u_old = u_new;
			h_old = h_new;
			iter_whl += 1;
		}
		//cout << ' ' << u_new[10] << ' ' << u_new[25] << endl;
		if (k == 50) { for (double a : u_new) cout << a << endl; }
		int aa = 0;
		/*bool check = check_invariant(x_setka[k + 1], u_new, q, h_old);
		if (!check)
		{
			cout << "BADBADBAD";
		}*/
		//end of whle cycle
		for (int i_ = 0; i_ < N_y; i_++) {
			if (isnan(u_new[i_])) {
				int nenn = 1;
			}
		}
		u[k + 1] = u_new;
		h[k + 1] = h_old;



	} //for k cycle end
	ofstream X("X.csv");
	ofstream Y("Y.csv");
	ofstream Z("Z.csv");
	ofstream Q("Q.csv");
	ofstream H_dump("H.csv"); //linii poverhnostei ravnih rashodov. each line corresponds to one row in csv 
	for (int i = 0; i <= N_y + 1; i++) {
		for (int k = 0; k <= N_x + 1; k++) {
			X << x_setka[k] << ";";
			Y << h[k][i] << ";";
			H_dump << h[k][i] << ";";
			Z << u[k][i] << ";";
		}

	

			
		X << endl;
		Y << endl;
		Z << endl;
		H_dump << endl;
	}

	for (int k = 0; k < N_y + 1; k++) {
		Q << q[k] << ";";
	}
	Q << endl;

}


int main() {
	cout << "Water started";
	//test_progonka();
	Solve_equation();
	return 0;
}	