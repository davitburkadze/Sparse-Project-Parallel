#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include"chrono"
#include"fstream"
#include "CgSparseParallel.h" 
#include "MatrixData.h"
#include "Helper.h" 


//#include "newTest.h" 
#include <experimental/filesystem> // C++-standard header file name  
#include <filesystem> // Microsoft-specific implementation header file name  
#include <regex>
#include <algorithm>

using namespace std::experimental::filesystem::v1;
using namespace std;


int main()
{

	

	// Random number generator for Ys
	random_device rd; // Random at each execution
	default_random_engine dre(rd());


	// Choose the generation type of Ys
	int genType = -1;
	cout << "1 - Existing Ys\n2 - Random Ys\n";
	while (true) {
		cin >> genType;
		if (genType > 2 || genType < 1) {
			cout << "Wrong number! Try again." << endl;
			continue;
		}
		else {
			break;
		}
	}

	// JSON RESULT
	std::ofstream jsonOfs("data.json");
	jsonOfs << "[";

	// Get the files with extension .mtx (specify full or relative path)
	string patternPath = "matrices\\";
	path current_dir(patternPath +"test");
	std::regex pattern("(.*\\.mtx)");

	
	//Collect basic information about matrix files (path, name, n, nnz)
	
	std::vector<MatrixData> matrixDatas = std::vector<MatrixData>();

	for (recursive_directory_iterator iter(current_dir), end; iter != end; ++iter)
	{
		string name = iter->path().filename().string();
		if (regex_match(name, pattern))
		{

			string path = iter->path().string();

			// Construct matrix data object
			MatrixData tmpData;
			tmpData.name = name;
			tmpData.path = path;

			bool skipMatrix = false;
			setMatrixDataExtras(tmpData, skipMatrix);

			// n1 != n2
			if (skipMatrix) {
				continue;
			}

			uniform_real_distribution<double> di(tmpData.min, tmpData.max);

			// Create '*_Ys.mtx' file for matrix (if it doesn't exist).
			string partPath = path.substr(0, path.find(".mtx", 0));
			string yPath = partPath + "_Ys.txt";
			tmpData.yPath = yPath; // Save path in MatrixData object
			bool forceOverwrite = false; // Be careful with this parameter!
			if (!exists(yPath) || forceOverwrite) {
				ofstream ofs(yPath);
				for (int i = 0; i < tmpData.n; i++) {
					ofs << di(dre) << '\n';
				}
				ofs.close();
			}

			matrixDatas.push_back(tmpData);
		}
	}

	// Sort by N
	sort(matrixDatas.begin(), matrixDatas.end(), MatrixData::compareByN);

	cout << matrixDatas.size() << " matrices..." << endl;

	try {
		// Loop through matrices and do the calculations
		for (unsigned int i = 0; i < matrixDatas.size(); i++)
		{
			string name = matrixDatas[i].name;
			string path = matrixDatas[i].path;
			int n = matrixDatas[i].n;
			int nnz = matrixDatas[i].nnz;

			double *a, *x0;

			// მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
			a = (double *)malloc(sizeof(double)*n);
			x0 = (double *)calloc(n, sizeof(double)); // allocates memory and fills with zeroes.

			if (genType == 1) {
				// FIXME Change to fast reading!!!!!
				std::ifstream ifs(matrixDatas[i].yPath);
				for (int i = 0; i < n; i++)
				{
					// Globally available static variables
					ifs >> a[i];
				}
				ifs.close();
			}
			else {
				uniform_real_distribution<double> di(matrixDatas[i].min, matrixDatas[i].max);
				for (unsigned int h = 0; h < n; h++) {
					a[h] = di(dre);
				}
			}

			cout << (i + 1) << endl;
			CgSparseParallel* abs = nullptr;
			pair<long, long> workingTimes;

			jsonOfs << "{";
			jsonOfs << "\"name\":\"" << name << "\",";
			jsonOfs << "\"n\":" << n << ",";
			jsonOfs << "\"nnz\":" << nnz << ",";

			jsonOfs << "\"results\":{";

			// sparse
			cout << path << ", N: " << n << ", NNZ: " << nnz << endl;
		
			
			abs = new CgSparseParallel(matrixDatas[i], a, x0);

			workingTimes.first = abs->fillSparseData();
			workingTimes.second = abs->getMinimal();

			cout << "[CG_SPARSE] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_sparse\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";

			delete abs;
			free(x0);
			x0 = (double *)calloc(n, sizeof(double)); // allocates memory and fills with zeroes.

			abs = new CgSparseParallel(matrixDatas[i], a, x0);

			workingTimes.first = abs->fillSparseDataLongWay();
			workingTimes.second = abs->getMinimalParallel();

			cout << "[CG_SPARSE_PARRALLEL] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_sparse_parrallel\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";

			delete abs;

			abs = new CgSparseParallel(matrixDatas[i], a, x0);
			workingTimes.first = abs->fillSparseDataParrallel();


			cout << "[FILL_PARRALLEL] : " << workingTimes.first << endl;
			jsonOfs << "\"fill_parrallel\":" <<workingTimes.first;

			cout << endl;

			jsonOfs << "}}";
			if (i != matrixDatas.size() - 1) {
				jsonOfs << ",";
			}

			jsonOfs.flush();

			// Free global variables for next iteration
			free(a);
			free(x0);
		}
		jsonOfs << "]";
		jsonOfs.close();

	}
	catch (std::overflow_error e) {
		cout << "ERROR!" << endl;
	}

	return 0;

}
