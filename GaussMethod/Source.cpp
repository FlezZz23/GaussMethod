#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <list>
#include <vector>
using namespace std;

void divide(vector < vector < double>>& A, int n, double x) {
	for (int i = 0; i < A[n].size(); i++) {
		if (x == 0) {
			cout << "Error, division on null" << endl;
		}
		else {
			A[n][i] /= x;
		}
	}
}

void substract(vector < vector < double>>& A, int n, int m, double x) {

	for (int j = 0; j < A[m].size(); j++) {
		A[m][j] -= A[n][j] * x;
	}
}

void swapr(vector < vector < double>>& A, int n, int m) {
	for (int i = 0; i < A[n].size(); i++) {
		swap(A[n][i], A[m][i]);
	}
}

void transform(vector < vector < double>>& A) {
	for (int i = 0; i < A.size(); i++) {
		divide(A, i, A[i][i]);
#pragma omp parallel for schedule(guided)
		for (int j = 0; j < A.size(); j++) {
			if (j == i) {
				continue;
			}
			substract(A, i, j, A[j][i]);
		}
	}
}
void transformNew(vector < vector < double>>& A) {
	double Max;
	int NumMax = 0;
	for (int i = 0; i < A.size(); i++) {
		Max = 0;
		for (int k = i; k < A.size(); k++) {
			if (A[k][i] > Max) {
				Max = A[k][i];
				NumMax = k;
			}
		}
		swapr(A, NumMax, i);
		divide(A, i, A[i][i]);
#pragma omp parallel for schedule(guided)
		for (int j = 0; j < A.size(); j++) {
			if (j == i) {
				continue;
			}
			substract(A, i, j, A[j][i]);
		}
	}
}


void Gauss(vector < vector < double>>& A) {
	if (A.size() + 1 != A[0].size()) {
		cout << "Error" << endl;
	}
	else
	{
		transform(A);
		for (int i = 0; i < A.size(); i++) {
			cout << "x" << i + 1 << "= " << A[i][A.size()] << endl;
		}
	}
}
void GaussNew(vector < vector < double>>& A) {
	if (A.size() + 1 != A[0].size()) {
		cout << "Error" << endl;
	}
	else
	{
		transformNew(A);
		for (int i = 0; i < A.size(); i++) {
			cout << "x" << i + 1 << "= " << A[i][A.size()] << endl;
		}
	}
}
//Последовательная без выбора главного элемента: 135ms
//Параллельная без выбора главного элемента: 50 ms
//Последовательная c выбором главного элемента по столбцам: 120 ms
//Параллельная c выбором главного элемента по столбцам: 50 ms
int main() {
	int n = 100;

	vector < vector < double > > A(n, vector <double>(n + 1));
	//for (int i = 0; i < A.size(); i++) {
	//	for (int j = 0; j < A[0].size(); j++) {
	//		cin >> A[i][j];
	//	}
	//}
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[0].size() - 1; j++) {
			A[A.size() - i - 1][i] = (double)(i + 1);
			A[A.size() - i - 1][A.size()] = (double)i;
			//A[i][i] = (double)(i + 1);
			//A[i][A.size()] = (double)i;
			if ((A.size() - i - 1) != j) {
				A[i][j] = 0;
			}
		}
	}

	double t = omp_get_wtime();
	GaussNew(A);
	cout << "Time: " << (omp_get_wtime() - t) * 1000 << " ms" << endl;
	return 0;
}