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

void EigenProblem(double distance, std::vector<double> &values, bool hasInteraction)
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

									double Cmat_el = 0.0;

									if (hasInteraction)
									{
										Cmat_el = Overlaps::CoulombMatel(alb_in, alk_in, beb_in, bek_in, i1b, i1k, i2b, i2k, distance);
									}

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

	for (int i = 0; i < matsize; i++)
	{
		double eigenvalue = gsl_vector_get(eigenvalues, i);
		values.push_back(eigenvalue);
	}

	gsl_eigen_gensymm_free(w);
	gsl_vector_free(eigenvalues);
	gsl_matrix_free(Hmat);
	gsl_matrix_free(Emat);

	steady_clock::time_point endTime = steady_clock::now();
	duration<float> elapsed = endTime - startTime;
	printf("Time elapsed point: %g\n", elapsed.count());
}

int main()
{
	steady_clock::time_point started = steady_clock::now();

	std::ofstream outputfileInt("EigenValuesIntSorted.dat");
	if (!outputfileInt.is_open())
	{
		std::cout << "Output file could not be opened! Terminating!" << std::endl;
		return 1;
	}
	std::ofstream outputfileNoInt("EigenValuesNoIntSorted.dat");
	if (!outputfileNoInt.is_open())
	{
		std::cout << "Output file could not be opened! Terminating!" << std::endl;
		return 1;
	}

	for (double distance = 1.0; distance <= 10.0; distance += 0.1) 
	{
		std::vector<double> values;
		EigenProblem(distance, values, true);

		outputfileInt << std::setprecision(2) << std::fixed;
		outputfileInt << distance << '\t';

		for (int i = 0; i < values.size(); i++)
		{
			outputfileInt << std::setprecision(10) << std::fixed;
			outputfileInt << values[i] << '\t';
//			printf("eigenvalue =%g\n", values[i]);
		}
		outputfileInt << '\n';
		values.clear();
	}

	for (double distance = 1.0; distance <= 10.0; distance += 0.1)
	{
		std::vector<double> values;
		EigenProblem(distance, values, false);

		outputfileNoInt << std::setprecision(2) << std::fixed;
		outputfileNoInt << distance << '\t';

		for (int i = 0; i < values.size(); i++)
		{
			outputfileNoInt << std::setprecision(10) << std::fixed;
			outputfileNoInt << values[i] << '\t';
//			printf("eigenvalue =%g\n", values[i]);
		}
		outputfileNoInt << '\n';
		values.clear();
	}

	steady_clock::time_point finished = steady_clock::now();
	duration<float> elapsed = finished - started;
	printf("TOTAL TIME ELAPSED: %g\n", elapsed.count());
}
