#include <iostream>
#include <fstream>
#include <stdio.h>
#include <complex>
#include <Python.h>

using namespace std;
int main(){

	/*
	Python rubbish
	*/
	char filename[]="matrixtest.py";
	FILE* fp;
	Py_Initialize();
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append(\".\")");
	fp= _Py_fopen(filename, "r");
	PyRun_SimpleFile(fp, filename);

	Py_Finalize();	
	// end of python rubbish
	// creates a file called temp.txt with data in
	// now reading it in
	
	cout << "\nPython stuff has finished\n" << endl;
	
	ifstream mfile;
	mfile.open("temp.txt");

	if ( !mfile ) exit( 1);

	int row1, col1;

	mfile >> row1 >> col1;
	int var[10][10];
	

	cout << row1 << endl;
	cout << col1 << endl;

	for (int i = row1 - 1; i >= 0; i--)
    		cout << a[i][0] << "\n";
	cout << a << "\n";
    	
	return 0;
}

