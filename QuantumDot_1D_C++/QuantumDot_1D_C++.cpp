#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <chrono>
#include "Overlaps.h"
#include "Parameters.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_sort_vector.h"
#include "FilesStrings.h"
#include <vector>

using namespace Parameters;
using namespace std::chrono;

void EigenProblem(double distance, std::vector<double> &values)
{

	steady_clock::time_point startTime = steady_clock::now();

	double* Vmat = new double[nmax * nmax * 2 * 2]();
	double* Smat = new double[nmax * nmax * 2 * 2]();

	const int matsize = 4 * nmax * nmax;

	gsl_matrix* Hmat = gsl_matrix_calloc(matsize, matsize);
	gsl_matrix* Emat = gsl_matrix_calloc(matsize, matsize);

	for (int ib = 0; ib < nmax; ib++)
	{
		for (int ik = 0; ik < nmax; ik++)
		{
			for (int alb = 0; alb < 2; alb++)
			{
				for (int alk = 0; alk < 2; alk++)
				{
					int alk_in = alk + 1;
					int alb_in = alb + 1;
					int count = Overlaps::GetCount(ib, ik, alb, alk, nmax, nmax, 2, 2);
					Vmat[count] = Overlaps::VIntegrate(alk_in, alb_in, ik, ib, distance);
					Smat[count] = Overlaps::SIntegrate(alk_in, alb_in, ik, ib, distance);
				}
			}
		}
	}

	for (int i1b = 0; i1b < nmax; i1b++)
	{
		for (int i1k = 0; i1k < nmax; i1k++)
		{
			for (int i2b = 0; i2b < nmax; i2b++)
			{
				for (int i2k = 0; i2k < nmax; i2k++)
				{
					for (int alb = 0; alb < 2; alb++)
					{
						for (int alk = 0; alk < 2; alk++)
						{
							for (int beb = 0; beb < 2; beb++)
							{
								for (int bek = 0; bek < 2; bek++)
								{
									int alb_in = alb + 1;
									int alk_in = alk + 1;
									int beb_in = beb + 1;
									int bek_in = bek + 1;

									int i = i1b * nmax * 2 * 2 +
											i2b * 2 * 2 +
											alb * 2 +
											beb;

									int j = i1k * nmax * 2 * 2 +
											i2k * 2 * 2 +
											alk * 2 +
											bek;

									double Cmat_el = Overlaps::CoulombMatel(alb_in, alk_in, beb_in, bek_in, i1b, i1k, i2b, i2k, distance);

									int count_a = Overlaps::GetCount(i1b, i1k, alb, alk, nmax, nmax, 2, 2);
									int count_b = Overlaps::GetCount(i2b, i2k, beb, bek, nmax, nmax, 2, 2);
									double ez1 = Overlaps::Eharmonic(i1k);
									double ez2 = Overlaps::Eharmonic(i2k);

									double Hmat_el = (ez1 + ez2) * Smat[count_a] * Smat[count_b] + Vmat[count_a] * Smat[count_b] + Smat[count_a] * Vmat[count_b];

									gsl_matrix_set(Hmat, i, j, Hmat_el + Cmat_el);

									double Emat_el = Smat[count_a] * Smat[count_b];

									gsl_matrix_set(Emat, i, j, Emat_el);

								}
							}
						}
					}
				}
			}
		}
	}


	delete[]Vmat;
	delete[]Smat;

	gsl_vector* eigenvalues = gsl_vector_alloc(matsize);
	gsl_eigen_gensymm_workspace* w = gsl_eigen_gensymm_alloc(matsize);

	gsl_eigen_gensymm(Hmat, Emat, eigenvalues, w);
	gsl_sort_vector(eigenvalues);

//	std::ofstream outputfile("EigenValuesSorted.txt");
//	if (!outputfile.is_open())
//	{
//		std::cout << "Output file could not be opened! Terminating!" << std::endl;
//		return 1;
//	}

	for (int i = 0; i < matsize; i++)
	{
		double eigenvalue = gsl_vector_get(eigenvalues, i);
	//	outputfile << std::setprecision(10) << std::fixed;
	//	outputfile << eigenvalue << '\n';
	//	printf("eigenvalue =%g\n", eigenvalue);
		values.push_back(eigenvalue);
	}

	gsl_eigen_gensymm_free(w);
	gsl_vector_free(eigenvalues);
	gsl_matrix_free(Hmat);
	gsl_matrix_free(Emat);

	steady_clock::time_point endTime = steady_clock::now();
	duration<float> elapsed = endTime - startTime;
	printf("Time elapsed since simulation start: %g\n", elapsed.count());
}

int main()
{
	std::ofstream outputfile("EigenValuesSorted.txt");
	if (!outputfile.is_open())
	{
		std::cout << "Output file could not be opened! Terminating!" << std::endl;
		return 1;
	}

	for (double distance = 1.0; distance <= 10.0; distance += 0.2) 
	{
		std::vector<double> values;
		EigenProblem(distance, values);

		for (int i = 0; i < values.size(); i++)
		{
			outputfile << std::setprecision(10) << std::fixed;
			outputfile << values[i] << '\t';
			printf("eigenvalue =%g\n", values[i]);
		}
		outputfile << '\n';
		values.clear();
	}
}
