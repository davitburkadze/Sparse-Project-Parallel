#include <iostream>
#include <string>
#include<random>
#include <fstream>
#include <chrono>
#include<future>
#include "MatrixData.h"

#include "CgSparseParallel.h"

using namespace std;


CgSparseParallel::CgSparseParallel(MatrixData matrixData, double* a, double* x0) {
	this->matrixData = matrixData;
	this->a = a;
	this->x0 = x0;

}
CgSparseParallel::~CgSparseParallel() {}


void CgSparseParallel::arrayCpy(double *x, double *y, const int n)
{
	int i, n5;
	if (n <= 0) return;
	n5 = n % 5;
	for (i = 0; i < n5; i++)  y[i] = x[i];
	for (; i < n; i += 5)
	{
		y[i] = x[i]; y[i + 1] = x[i + 1]; y[i + 2] = x[i + 2];  y[i + 3] = x[i + 3]; y[i + 4] = x[i + 4];
	}
}
void CgSparseParallel::intArrayCpy(int *x, int *y, const int n)
{
	int i, n5;
	if (n <= 0) return;
	n5 = n % 5;
	for (i = 0; i < n5; i++)  y[i] = x[i];
	for (; i < n; i += 5)
	{
		y[i] = x[i]; y[i + 1] = x[i + 1]; y[i + 2] = x[i + 2];  y[i + 3] = x[i + 3]; y[i + 4] = x[i + 4];
	}
}

//დამხმარე ფუნქცია მასივების სწრაფი გამრავლებისთვის
inline	double CgSparseParallel::vecProd(double *x, double *y, const int n)
{
	int i, n5;
	double sum(0.0);
	if (n <= 0) return sum;
	n5 = n % 5;
	for (i = 0; i < n5; i++) sum += x[i] * y[i];
	for (; i < n; i += 5)
	{
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
			+ x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}



long CgSparseParallel::fillSparseData()
{
	auto st = chrono::high_resolution_clock::now();

	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char * buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		cout << "the file does not exist!" << endl;

	/* Get the number of bytes */
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);

	/* reset the file position indicator to
	the beginning of the file */
	fseek(infile, 0L, SEEK_SET);

	/* grab sufficient memory for the
	buffer to hold the text */
	buffer = (char*)calloc(numbytes, sizeof(char));

	/* memory error */
	if (buffer == NULL)
		cout << "memory error!" << endl;

	/* copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);

	// Ignore comment section
	size_t pos = 0;
	char *data = buffer;
	while (data[pos] == '%')
	{
		++pos;
		while (data[pos] != '\n')
			++pos;
		data += (pos + 1);
		pos = 0;
	}

	//რაოდენობები წავიკითხეთ
	while (data[pos] != ' ') ++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	// There is second n in the matrix file
	pos = 0;
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	pos = 0;
	while (data[pos] != '\n') ++pos;
	data[pos] = '\0';
	nnzz = (int)atoi(data);
	++pos;
	data += pos;
	pos = 0;
	//დავიწყოთ მეჩხერი არანულოვანი ქვემატრიცის შევსება
	//მეხსიერების გამოყოფა
	A = (double **)malloc(sizeof(double *)*n);
	Ind = (int**)malloc((n + 1) * sizeof(int*));
	Ind[0] = (int*)malloc((n + 1) * sizeof(int));
	Ind[0][0] = n;

	//მაქს. სიგრძის ორი ვექტორი, სტრიქონში ინდექსის და მნიშვნელობებისთვის
	int* ind = (int *)malloc(sizeof(int)*n);
	double*	val = (double *)malloc(sizeof(double)*n);

	int row = 1;
	int c = 0;		//მთვლელი
					//ფაილიდან სტრიქონის წამოსაღები 
	double v; 	int i, j;
	//ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < nnzz - 1; ++ii)
	{
		//წავიკითხოთ ასეთი რიგით: j,i,v
		//j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		while (data[pos] != '\n') ++pos;
		data[pos] = '\0';
		v = (double)atof(data);
		++pos;
		data += pos;
		pos = 0;
		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));

			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			////////cout << c << endl;////////
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
		/**/
	}

	{//ბოლო სტრიქონი
	 //წავიკითხოთ ასეთი რიგით: j,i,v
	 //j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		v = (double)atof(data);
		//მოვრჩით კითხვას..

		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));
			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
	}

	A[n - 1] = (double*)malloc(c * sizeof(double));
	Ind[n] = (int*)malloc(c  * sizeof(int));
	arrayCpy(val, A[n - 1], c);
	intArrayCpy(ind, Ind[n], c);
	Ind[0][n] = c;

	delete[] ind;
	delete[] val;
	free(buffer);

	/*//მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
	a = (double *)malloc(sizeof(double)*n);
	x0 = (double *)malloc(sizeof(double)*n);
	default_random_engine dre;
	uniform_real_distribution<double> di(0, 20);
	for (int i = 0; i < n; i++)
	{
	a[i] = di(dre);			//ფაილში რომ წერია?
	x0[i] = 0.;
	}*/
	auto diff = chrono::high_resolution_clock::now() - st;
	return chrono::duration_cast<chrono::milliseconds>(diff).count();
}


long CgSparseParallel::fillSparseDataLongWay()
{
	auto st = chrono::high_resolution_clock::now();
	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char * buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		std::cout << "the file does not exist!" << endl;

	/* Get the number of bytes */
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);

	/* reset the file position indicator to
	the beginning of the file */
	fseek(infile, 0L, SEEK_SET);

	/* grab sufficient memory for the
	buffer to hold the text */
	buffer = (char*)calloc(numbytes, sizeof(char));

	/* memory error */
	if (buffer == NULL)
		std::cout << "memory error!" << endl;

	/* copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);

	// Ignore comment section
	size_t pos = 0;
	char *data = buffer;
	while (data[pos] == '%')
	{
		++pos;
		while (data[pos] != '\n')
			++pos;
		data += (pos + 1);
		pos = 0;
	}

	//რაოდენობები წავიკითხეთ
	while (data[pos] != ' ') ++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	// There is second n in the matrix file
	pos = 0;
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	pos = 0;
	while (data[pos] != '\n') ++pos;
	data[pos] = '\0';
	nnzz = (int)atoi(data);
	++pos;
	data += pos;

	//დავიწყოთ მეჩხერი არანულოვანი ქვემატრიცის შევსება
	//მეხსიერების გამოყოფა
	A = (double **)malloc(sizeof(double *)*n);
	Ind = (int**)malloc((n + 1) * sizeof(int*));
	Ind[0] = (int*)malloc((n + 1) * sizeof(int));
	Ind[0][0] = n;

	//მაქს. სიგრძის ორი ვექტორი, სტრიქონში ინდექსის და მნიშვნელობებისთვის
	int* ind = (int *)malloc(sizeof(int)*n);
	double*	val = (double *)malloc(sizeof(double)*n);
	//სტრიქონების დიაგონალამდე ელემენტების რაოდენობის მთვლელების მასივი შევქმნათ
	//და გავუნულოთ მნიშვნელობები
	int* rcounter = (int *)malloc(sizeof(int)*n);
	for (pos = 0; pos < n; ++pos)
		rcounter[pos] = 0;
	pos = 0;

	int row = 1;
	int c = 0;		//მთვლელი
					//ფაილიდან სტრიქონის წამოსაღები 
	double v; 	int i, j;
	int shifted;	//ქვედა სამკუთხას ჩასმისთვის ადგილების დასატოვებლად
					//ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < nnzz - 1; ++ii)
	{
		//წავიკითხოთ ასეთი რიგით: j,i,v
		//j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++rcounter[j - 1];
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		while (data[pos] != '\n') ++pos;
		data[pos] = '\0';
		v = (double)atof(data);
		++pos;
		data += pos;
		pos = 0;
		if (row != i)
		{
			--row;
			shifted = rcounter[row];
			--shifted;		//მთვლელს ყველა სტრიქონში ზედმეტად დაემატა დიაგონალური
							//std::cout << shifted << endl;			//////////////
			A[row] = (double*)malloc((c + shifted) * sizeof(double));
			Ind[row + 1] = (int*)malloc((c + shifted) * sizeof(int));
			arrayCpy(val, A[row] + shifted, c);
			for (int k = 0; k < shifted; ++k)   A[row][k] = 0;	/////////////
			intArrayCpy(ind, Ind[row + 1] + shifted, c);
			for (int k = 0; k < shifted; ++k)   Ind[row + 1][k] = 0;	/////////////
			Ind[0][row + 1] = c + shifted;
			row = i;
			////////cout << c << endl;////////
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
		/**/
	}

	{//ბოლო ელემენტი - რიგით: j,i,v
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++rcounter[j - 1];
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		v = (double)atof(data);
		//მოვრჩით კითხვას..

		if (row != i)
		{
			--row;
			shifted = rcounter[row];
			--shifted;		//მთვლელს ყველა სტრიქონში ზედმეტად დაემატა დიაგონალური
							//std::cout << shifted << endl;			//////////////
			A[row] = (double*)malloc((c + shifted) * sizeof(double));
			Ind[row + 1] = (int*)malloc((c + shifted) * sizeof(int));
			arrayCpy(val, A[row] + shifted, c);
			for (int k = 0; k < shifted; ++k)   A[row][k] = 0;	/////////////
			intArrayCpy(ind, Ind[row + 1] + shifted, c);
			for (int k = 0; k < shifted; ++k)   Ind[row + 1][k] = 0;	/////////////
			Ind[0][row + 1] = c + shifted;
			row = i;
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
	}
	shifted = rcounter[n - 1];
	--shifted;	//მთვლელს ბოლო სტრიქონშიც ზედმეტად დაემატა დიაგონალური
				//std::cout << shifted << endl;			//////////////
	A[n - 1] = (double*)malloc((c + shifted) * sizeof(double));
	Ind[n] = (int*)malloc((c + shifted) * sizeof(int));
	arrayCpy(val, A[n - 1] + shifted, c);
	for (int k = 0; k < shifted; ++k)   A[n - 1][k] = 0;	/////////////
	intArrayCpy(ind, Ind[n] + shifted, c);
	for (int k = 0; k < shifted; ++k)   Ind[n][k] = 0;	/////////////
	Ind[0][n] = c + shifted;

	delete[] ind;
	delete[] val;
	free(buffer);

	//მთვლელი, რომელიც Ind-ის ქვედა სამუკუთხა ნაწილის შევსებისთვისაა საჭირო

	int* which = (int *)malloc(sizeof(int)*n);
	for (pos = 0; pos < n; ++pos)
		which[pos] = 0;

	//ქვედა სამკუთხა მატრიცა ჩავსვათ  
	for (int i = 0; i < n; ++i)
	{
		for (int j = rcounter[i]; j < Ind[0][i + 1]; ++j)
		{
			A[Ind[i + 1][j]][which[Ind[i + 1][j]]] = A[i][j];
			Ind[1 + Ind[i + 1][j]][which[Ind[i + 1][j]]] = i;
			++which[Ind[i + 1][j]];
		}
	}

	/*//მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
	a = (double *)malloc(sizeof(double)*n);
	x0 = (double *)malloc(sizeof(double)*n);
	default_random_engine dre;
	uniform_real_distribution<double> di(0, 20);
	for (int i = 0; i < n; i++)
	{
	a[i] = di(dre);			//ფაილში რომ წერია?
	x0[i] = 0.;
	}*/
	auto diff = chrono::high_resolution_clock::now() - st;
	return chrono::duration_cast<chrono::milliseconds>(diff).count();
}


//დამხმარე ფუნქცია სიმეტრიული მატრიცის და მასივის სწრაფი გამრავლებისთვის
void CgSparseParallel::MatrixByVector(double **m, int **index, double *x, double* res)
{
	const int k(Ind[0][0]);
	
	int i, j;
	double *p;
	int *q;
	int size;

	for (i = 0; i < k; ++i)
		res[i] = 0.;

	for (i = 0; i < k; ++i)
	{
		p = m[i];
		q = index[i + 1];
		size = index[0][i + 1];
		res[i] += p[0] * x[q[0]];

		for (j = 1; j < size; ++j)
		{
			res[i] += p[j] * x[q[j]];
			res[q[j]] += p[j] * x[i];
		}
	}
}

void CgSparseParallel::fullMatrixByVector(double **m, int **index, double *x, double* res, int begin, int end)
{
	const int k(Ind[0][0]);

	int i, j;
	double *p;
	int *q;
	int size;

	for (i = begin; i < end; ++i)
		res[i] = 0.;

	for (i = begin; i < end; ++i)
	{
		p = m[i];
		q = index[i + 1];
		size = index[0][i + 1];

		for (j = 0; j < size; ++j)
		{
			res[i] += p[j] * x[q[j]];
		}
	}
}

void CgSparseParallel::fullMatrixByVectorParallel(double **m, int **index, double *x, double* res, int begin, int end, int chunckSize)
{
	unsigned int const length = end - begin;

	if (length <= chunckSize)
	{
		return fullMatrixByVector(m,index,x,res,begin,end);
	}
	else 
	{
		int mid = (begin + end) / 2;
	
		std::future<void> future = std::async(&CgSparseParallel::fullMatrixByVectorParallel,this, m, index, x, res, begin, mid, chunckSize);

		fullMatrixByVectorParallel(m, index, x, res, mid, end, chunckSize);
		future.get();
	}
}

long CgSparseParallel::getMinimal()
{

	auto st = chrono::high_resolution_clock::now();

	double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*n);
	d = (double *)malloc(sizeof(double)*n);
	q = (double *)malloc(sizeof(double)*n);

	for (int i = 0; i < n; i++)
		r[i] = d[i] = q[i] = 0.;

	for (int i = 0; i < n; i++)
		d[i] = r[i] = a[i];

	double deltaNew(vecProd(r, r, n));

	int count = 0;
	while (deltaNew >= EPSS*EPSS)
	{

		MatrixByVector(A, Ind, d, q);

		double alpha = deltaNew / vecProd(d, q, n);
		for (int i = 0; i < n; i++)
			x0[i] += alpha*d[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = vecProd(r, r, n);
		double beta(deltaNew / deltaOld);
		for (int i = 0; i < n; i++)
			d[i] = r[i] + (beta*d[i]);
		count++;
	}

	// Free memory
	free(r); r = NULL;
	free(d); d = NULL;
	free(q); q = NULL;

	auto diff = chrono::high_resolution_clock::now() - st;
	return chrono::duration_cast<chrono::milliseconds>(diff).count();
}

long CgSparseParallel::getMinimalParallel()
{

	auto st = chrono::high_resolution_clock::now();

	int chunckSize;
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();

	chunckSize = ((n / (concurentThreadsSupported)) < 1000) ? 1000 : (n /concurentThreadsSupported+ 1);

	double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*n);
	d = (double *)malloc(sizeof(double)*n);
	q = (double *)malloc(sizeof(double)*n);

	for (int i = 0; i < n; i++)
		r[i] = d[i] = q[i] = 0.;

	for (int i = 0; i < n; i++)
		d[i] = r[i] = a[i];

	double deltaNew(vecProd(r, r, n));

	int count = 0;
	while (deltaNew >= EPSS*EPSS)
	{
		fullMatrixByVectorParallel(A, Ind, d, q, 0, n, chunckSize);

		double alpha = deltaNew / vecProd(d, q, n);
		for (int i = 0; i < n; i++)
			x0[i] += alpha*d[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = vecProd(r, r, n);
		double beta(deltaNew / deltaOld);
		for (int i = 0; i < n; i++)
			d[i] = r[i] + (beta*d[i]);
		count++;
	}

	auto diff = chrono::high_resolution_clock::now() - st;
	return chrono::duration_cast<chrono::milliseconds>(diff).count();
}


long CgSparseParallel::fillSparseDataParrallel() {

	auto st = chrono::high_resolution_clock::now();

	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char * buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		cout << "the file does not exist!" << endl;

	/* Get the number of bytes */
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);

	/* reset the file position indicator to
	the beginning of the file */
	fseek(infile, 0L, SEEK_SET);

	/* grab sufficient memory for the
	buffer to hold the text */
	buffer = (char*)calloc(numbytes, sizeof(char));

	/* memory error */
	if (buffer == NULL)
		cout << "memory error!" << endl;

	/* copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);


	size_t pos = 0;
	char *data = buffer;

	// Ignore comment section
	while (data[pos] == '%')
	{
		++pos;
		while (data[pos] != '\n')
			++pos;
		data += (pos + 1);
		pos = 0;
	}

	//რაოდენობები წავიკითხეთ
	while (data[pos] != ' ') ++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	// There is second n in the matrix file
	pos = 0;
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	pos = 0;
	while (data[pos] != '\n') ++pos;
	data[pos] = '\0';
	nnzz = (int)atoi(data);
	++pos;
	data += pos;
	pos = 0;
	//დავიწყოთ მეჩხერი არანულოვანი ქვემატრიცის შევსება
	//მეხსიერების გამოყოფა
	A = (double **)malloc(sizeof(double *)*n);
	Ind = (int**)malloc((n + 1) * sizeof(int*));
	Ind[0] = (int*)malloc((n + 1) * sizeof(int));
	Ind[0][0] = n;
	int chunckSize;
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	chunckSize = nnzz / concurentThreadsSupported + 1;
	readMatrixParallel(data, numbytes, A, Ind, 0, nnzz, chunckSize, n, nnzz);

	/*//მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
	a = (double *)malloc(sizeof(double)*n);
	x0 = (double *)malloc(sizeof(double)*n);
	default_random_engine dre;
	uniform_real_distribution<double> di(0, 20);
	for (int i = 0; i < n; i++)
	{
	a[i] = di(dre);         //ფაილში რომ წერია?
	x0[i] = 0.;
	}*/
	auto diff = chrono::high_resolution_clock::now() - st;
	return chrono::duration_cast<chrono::milliseconds>(diff).count();
}




void CgSparseParallel::readMatrixParallel(char * data, long numbytes, double **A, int **Ind, int begin, int end, int chunckSize, int n, int nnzz) {
	unsigned int const length = end - begin;

	if (length <= chunckSize)
	{
		return readMatrixdata(data, A, Ind, begin, end, n, nnzz);
	}
	else
	{
		int mid = (begin + end) / 2;

		std::future<void> future = std::async(&CgSparseParallel::readMatrixParallel, this, data, numbytes, A, Ind, begin, mid, chunckSize, n, nnzz);

		char * dataCopy = (char*)calloc(numbytes, sizeof(char));

		memcpy(dataCopy, data, sizeof(char)*numbytes);

		readMatrixParallel(dataCopy, numbytes, A, Ind, mid, end, chunckSize, n, nnzz);
		future.get();
	}
}


void CgSparseParallel::readMatrixdata(char * data, double **A, int **Ind, int begin, int end, int n, int nnzz) {

	size_t pos = 0;
	int checkRow = 1;
	for (int ii = 0; ii < begin; ++ii)
	{

		if (ii == begin - 1) {
			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			int j = (int)atoi(data);
			++pos;
			data += pos;
			pos = 0;
			//i:
			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			int i = (int)atoi(data);
			++pos;
			data += pos;
			pos = 0;

			++pos;
			while (data[pos] != '\n')
				++pos;
			data += (pos + 1);
			pos = 0;

			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			j = (int)atoi(data);
			++pos;
			data += pos;
			pos = 0;
			//i:
			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			checkRow = (int)atoi(data);
			++pos;
			data += pos;

			//cout << "&&&" << i << "  " << checkRow << "&&&" << endl;
			if (i == checkRow) {
				++begin;
			}
			pos = 0;

			while (data[pos] != '\0') --pos;
			data[pos] = ' ';
			while (data[pos] != '\0') --pos;
			data[pos] = ' ';

			while (data[pos] != '\n')
				--pos;
			data += (pos + 1);
			pos = 0;


		}
		else {
			++pos;
			while (data[pos] != '\n')
				++pos;
			data += (pos + 1);
			pos = 0;
		}
	}


	if (end == nnzz) {
		end -= 1;
	}

	//მაქს. სიგრძის ორი ვექტორი, სტრიქონში ინდექსის და მნიშვნელობებისთვის
	int* ind = (int *)malloc(sizeof(int)*n);
	double* val = (double *)malloc(sizeof(double)*n);

	int row = checkRow;
	int c = 0;      //მთვლელი
					//ფაილიდან სტრიქონის წამოსაღები 
	double v;   int i, j;
	//ბოლო სტრიქონის გარდა
	for (int ii = begin; ii < end; ++ii)
	{
		//წავიკითხოთ ასეთი რიგით: j,i,v
		//j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		while (data[pos] != '\n') ++pos;
		data[pos] = '\0';
		v = (double)atof(data);

		//cout << "###" << i<<"  "<<j << "###" << endl;

		++pos;
		data += pos;
		pos = 0;
		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));

			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			////////cout << c << endl;////////
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;

		if (ii == end - 1 && end != nnzz - 1) {
			++pos;
			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			int j = (int)atoi(data);
			++pos;
			data += pos;
			pos = 0;
			//i:
			while (data[pos] != ' ') ++pos;
			data[pos] = '\0';
			int checkRow = (int)atoi(data);
			++pos;
			data += pos;
			pos = 0;

			//cout << "###" << i << "  " << checkRow << "###" << endl;

			if (i == checkRow) {
				++end;
			}
			else {
				--row;
				A[row] = (double*)malloc(c * sizeof(double));
				Ind[row + 1] = (int*)malloc(c * sizeof(int));

				arrayCpy(val, A[row], c);
				intArrayCpy(ind, Ind[row + 1], c);
				Ind[0][row + 1] = c;
				row = i;
				////////cout << c << endl;////////
				c = 0;

			}

			while (data[pos] != '\0') --pos;
			data[pos] = ' ';
			while (data[pos] != '\0') --pos;
			data[pos] = ' ';

			while (data[pos] != '\n')
				--pos;
			data += (pos + 1);
			pos = 0;
		}
		/**/
	}

	if (end == nnzz - 1) {//ბოლო სტრიქონი
						  //წავიკითხოთ ასეთი რიგით: j,i,v
						  //j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		v = (double)atof(data);
		//მოვრჩით კითხვას..

		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));
			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
		A[n - 1] = (double*)malloc(c * sizeof(double));
		Ind[n] = (int*)malloc(c  * sizeof(int));
		arrayCpy(val, A[n - 1], c);
		intArrayCpy(ind, Ind[n], c);
		Ind[0][n] = c;
	}
	delete[] ind;
	delete[] val;
}