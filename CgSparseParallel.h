
#include "MatrixData.h"

using namespace std;


class  CgSparseParallel {

public:
	const double EPSS = 1E-5;
	double **A;
	int **Ind;
	double *a, *x0;
	int n, nnzz;

	MatrixData matrixData;

	CgSparseParallel(MatrixData matrixData, double* a, double* x0);

	~CgSparseParallel();

	void arrayCpy(double *x, double *y, const int n);
	void intArrayCpy(int *x, int *y, const int n);

	long fillSparseData();
	long fillSparseDataLongWay();
	long fillSparseDataParrallel();
	void readMatrixdata(char * data, double **A, int **Ind, int begin, int end, int n, int nnzz);
	void readMatrixParallel(char * data, long numbytes, double **A, int **Ind, int begin, int end, int chunckSize, int n, int nnzz);

	//==============
	//დამხმარე ფუნქცია მასივების სწრაფი გამრავლებისთვის
	inline	double vecProd(double *x, double *y, const int n);

	//დამხმარე ფუნქცია სიმეტრიული მატრიცის და მასივის სწრაფი გამრავლებისთვის
	void MatrixByVector(double **m, int **index, double *x, double* res);

	//დამხმარე ფუნქცია სიმეტრიული სრული მატრიცის და მასივის სწრაფი გამრავლებისთვის
	void fullMatrixByVector(double **m, int **index, double *x, double* res, int begin, int end);
	void fullMatrixByVectorParallel(double **m, int **index, double *x, double* res, int begin, int end,int chunckSize);

	long getMinimal();
	long getMinimalParallel();
};