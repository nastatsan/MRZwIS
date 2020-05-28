#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

int p, m, q, n, r;
double p_double, m_double, q_double, n_double;

vector <vector<double>> A, B, E, G, Cij;
vector <vector<vector <double>>> Fijk, Dijk;

int count_of_calls_sum = 0, count_of_calls_diffference = 0, count_of_calls_multiplicity = 0, count_of_calls_comparing = 0;

int count_of_calls_implication_A_to_B = 0, count_of_calls_implication_B_to_A = 0, count_of_calls_composition_Fijk = 0,
count_of_calls_composition_Dijk = 0, count_of_calls_min_Fijk_and_Dijk = 0;						//???

int time_of_sum = 0, time_of_difference = 0, time_of_multiplicity = 0, time_of_comparing = 0;

double Tn = 0, Lsum = 0, T1 = 0, Ky = 0, e = 0, Lavg = 0;

vector <vector<double>> generateMatrix(int rows, int columns);
void printMatrix(int rows, int columns, vector <vector<double>> matrix, string s);
vector <vector<vector <double>>> calculate_Fijk();
vector <vector<vector <double>>> calculate_Dijk();
vector <vector<double>> calculate_Cij();
double implication_A_to_B(double a, double b);
double implication_B_to_A(double a, double b);
double composition_Fijk(int i, int j);
double composition_Dijk(int i, int j);
double min_Fijk_and_Dijk(int i, int j);
double random(double min, double max);


int main()
{
	system("cls");
	cout << "Enter p\n";
	cin >> p;
	p_double = (double)p;
	cout << "Enter m\n";
	cin >> m;
	m_double = (double)m;
	cout << "Enter q\n";
	cin >> q;
	q_double = (double)q;
	cout << "Enter n\n";
	cin >> n;
	n_double = (double)n;
	cout << "Enter time_of_sum\n";
	cin >> time_of_sum;
	cout << "Enter time_of_difference\n";
	cin >> time_of_difference;
	cout << "Enter time_of_multiplicity\n";
	cin >> time_of_multiplicity;
	cout << "Enter time_of_comparing\n";
	cin >> time_of_comparing;

	A = generateMatrix(p, m);
	B = generateMatrix(m, q);
	E = generateMatrix(1, m);
	G = generateMatrix(p, q);
	Dijk = calculate_Dijk();
	Fijk = calculate_Fijk();
	Cij = calculate_Cij();


	printMatrix(p, m, A, "Matrix A");
	printMatrix(m, q, B, "Matrix B");
	printMatrix(1, m, E, "Matrix E");
	printMatrix(p, q, G, "Matrix G");
	printMatrix(p, q, Cij, "Matrix C");


	T1 = count_of_calls_comparing * time_of_comparing + count_of_calls_diffference * time_of_difference
		+ count_of_calls_multiplicity * time_of_multiplicity + count_of_calls_sum * time_of_sum;

	Lsum = Tn;
	
	Ky = T1 / Tn;
	e = Ky / n;

	r = p * q * m;
	//calculate_Fijk
	Lavg += (7 * time_of_multiplicity + 2 * time_of_sum + 3 * time_of_difference) * (p * q * m);
	//calculate_Dijk
	Lavg += time_of_multiplicity * (p * q * m);
	//calculate_Cij
	Lavg += (7 * time_of_multiplicity + 2 * time_of_sum + 3 * time_of_difference) * (p * q);
	//implication_A_to_B
	Lavg += (time_of_multiplicity + time_of_sum + time_of_difference) * count_of_calls_implication_A_to_B;
	//implication_B_to_A
	Lavg += (time_of_multiplicity + time_of_sum + time_of_difference) * count_of_calls_implication_B_to_A;
	//composition_Fijk
	Lavg += time_of_multiplicity * m * count_of_calls_composition_Fijk;
	//composition_Dijk
	Lavg += (time_of_difference * (m + 1) + time_of_multiplicity * m) * count_of_calls_composition_Dijk;
	//min_Fijk_and_Dijk
	Lavg += (time_of_comparing) * count_of_calls_min_Fijk_and_Dijk;

	//итоговое значение
	Lavg /= r;

	double D = Tn / Lavg;
	cout << "\nT1 = " << T1 << "\nTn = " << Tn << "\nKy = " << Ky << "\ne = " << e << "\nLsum = " << Lsum
		<< "\nLavg = " << Lavg << "\nD = " << D << endl;
	return 0;
}


vector <vector<double>> generateMatrix(int rows, int columns) {
	vector <vector<double>> matrix;
	vector <double> vector_row;
	for (int i = 0; i < columns; i++)  vector_row.push_back(0);

	for (int row = 0; row < rows; row++) {
		matrix.push_back(vector_row);
		for (int column = 0; column < columns; column++)
			matrix[row][column] = random(-1, 1);
	}
	return matrix;
}

double random(double min, double max) {
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

void printMatrix(int rows, int columns, vector<vector<double>> matrix, string s) {
	cout << s << endl;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			cout << matrix[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}

vector <vector<vector <double>>> calculate_Fijk() {

	double time_of_operarion = 0;

	vector <vector<vector <double>>> Fijk;
	vector <vector<double>> Fijk_twodimensional;
	vector <double> Fijk_onedimensional;

	for (int k = 0; k < m; k++) Fijk_onedimensional.push_back(0);
	for (int j = 0; j < q; j++) Fijk_twodimensional.push_back(Fijk_onedimensional);
	for (int i = 0; i < p; i++) Fijk.push_back(Fijk_twodimensional);

	for (int i = 0; i < p; i++) {
		for (int j = 0; j < q; j++) {
			for (int k = 0; k < m; k++) {
				Fijk[i][j][k] = implication_A_to_B(A[i][k], B[k][j]) * (2 * E[0][k] - 1) * E[0][k] + implication_B_to_A(A[i][k], B[k][j])
					* (1 + (4 * implication_A_to_B(A[i][k], B[k][j]) - 2) * E[0][k]) * (1 - E[0][k]);
				count_of_calls_multiplicity += 7;
				count_of_calls_diffference += 3;
				count_of_calls_sum += 2;
			}
		}
	}

	time_of_operarion = (7 * time_of_multiplicity + 2 * time_of_sum + 3 * time_of_difference) * ceil(p_double * q_double * m_double / n_double);
	Tn += time_of_operarion;
	return Fijk;
}

vector <vector<vector <double>>> calculate_Dijk() {
	double time_of_operation = 0;

	vector <double> Dijk_onedimensional;
	vector < vector < double >> Dijk_twodimensional;
	vector <vector<vector <double>>> Dijk;


	for (int k = 0; k < m; k++) Dijk_onedimensional.push_back(0);
	for (int j = 0; j < q; j++) Dijk_twodimensional.push_back(Dijk_onedimensional);
	for (int i = 0; i < p; i++) Dijk.push_back(Dijk_twodimensional);

	for (int i = 0; i < p; i++) {
		for (int j = 0; j < q; j++) {
			for (int k = 0; k < m; k++) {
				Dijk[i][j][k] = A[i][k] * B[k][j];
				count_of_calls_multiplicity++;
			}
		}
	}

	time_of_operation = time_of_multiplicity * ceil(p_double * q_double * m_double / n_double);
	Tn += time_of_operation;
	return Dijk;
}

vector <vector<double>> calculate_Cij() {
	double time_of_operation = 0;

	vector <vector<double>> Cij;
	vector <double> Cij_onedimensional;

	for (int j = 0; j < q; j++) Cij_onedimensional.push_back(0);
	for (int i = 0; i < p; i++) Cij.push_back(Cij_onedimensional);


	for (int i = 0; i < p; i++) {
		for (int j = 0; j < q; j++) {
			Cij[i][j] = composition_Fijk(i, j) * (3 * G[i][j] - 2) * G[i][j] + (composition_Dijk(i, j) +
				(4 * min_Fijk_and_Dijk(i, j) - 3 * composition_Dijk(i, j)) * G[i][j]) * (1 - G[i][j]);
			count_of_calls_multiplicity += 7;
			count_of_calls_diffference += 3;
			count_of_calls_sum += 2;
		}
	}
	time_of_operation = (7 * time_of_multiplicity + 2 * time_of_sum + 3 * time_of_difference) * ceil(p_double * q_double / n_double);

	Tn += time_of_operation;
	return Cij;
}

double implication_A_to_B(double a, double b) {
	double time_of_operation = 0;
	time_of_operation = time_of_multiplicity + time_of_sum + time_of_difference;
	Tn += time_of_operation;

	count_of_calls_diffference++;
	count_of_calls_sum++;
	count_of_calls_multiplicity++;

	count_of_calls_implication_A_to_B++;
	return a * (b - 1) + 1;
}

double implication_B_to_A(double a, double b) {
	double time_of_operation = 0;
	time_of_operation = time_of_multiplicity + time_of_sum + time_of_difference;
	Tn += time_of_operation;

	count_of_calls_diffference++;
	count_of_calls_sum++;
	count_of_calls_multiplicity++;

	count_of_calls_implication_B_to_A++;

	return b * (a - 1) + 1;
}

double composition_Fijk(int i, int j) {
	double time_of_operation = 0;
	double result = 1;

	for (int k = 0; k < m; k++) {
		result *= Fijk[i][j][k];
		count_of_calls_multiplicity++;
	}

	time_of_operation = time_of_multiplicity * m;
	Tn += time_of_operation;
	
	count_of_calls_composition_Fijk++;
	return result;
}

double composition_Dijk(int i, int j) {
	double time_of_operation = 0;
	double result = 1;

	for (int k = 0; k < m; k++) {
		result *= 1 - Dijk[i][j][k];
		count_of_calls_diffference++;
		count_of_calls_multiplicity++;
	}

	count_of_calls_diffference++;
	time_of_operation = time_of_difference * (m + 1) + time_of_multiplicity * m;
	Tn += time_of_operation;

	count_of_calls_composition_Dijk++;

	return 1 - result;
}

double min_Fijk_and_Dijk(int i, int j) {
	double time_of_operation = 0;

	count_of_calls_comparing++;

	time_of_operation = time_of_comparing;
	Tn += time_of_operation;

	count_of_calls_min_Fijk_and_Dijk++;

	if (composition_Fijk(i, j) > composition_Dijk(i, j)) return composition_Dijk(i, j);
	else return composition_Fijk(i, j);
}
