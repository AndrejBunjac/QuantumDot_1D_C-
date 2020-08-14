#include "Overlaps.h"
#include <cmath>
#include <math.h>
#include "Parameters.h"
#include "gsl/gsl_sf_hermite.h"
#include "BasicFunctions.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_gamma.h"

using namespace Parameters;

struct integration_params { int alk; int alb; int ik; int ib; };

struct C_integration_params_inner { int alb; int alk; int beb; int bek; int n1b; int n1k; int n2b; int n2k; double z_2; };
struct C_integration_params_outer { int alb; int alk; int beb; int bek; int n1b; int n1k; int n2b; int n2k; };


double V_integration_function(double z, void* p)
{
	struct integration_params* params = (struct integration_params*)p;
	int alk = (params->alk);
	int alb = (params->alb);
	int ik = (params->ik);
	int ib = (params->ib);

	double zk = (alk == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double zb = (alb == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double phiket = Overlaps::phi(ik, zk);
	double phibra = Overlaps::phi(ib, zb);
	double dv = Overlaps::DV(alk, z);
	return phiket * dv * phibra;
}


double S_integration_function(double z, void* p)
{
	struct integration_params* params = (struct integration_params*)p;
	int alk = (params->alk);
	int alb = (params->alb);
	int ik = (params->ik);
	int ib = (params->ib);

	double zk = (alk == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double zb = (alb == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double phiket = Overlaps::phi(ik, zk);
	double phibra = Overlaps::phi(ib, zb);
	return phiket * phibra;
}

double C_integration_function_inner(double z_1, void* p)
{
	struct C_integration_params_inner* params = (struct C_integration_params_inner*)p;
	int alb = (params->alb);
	int alk = (params->alk);
	int beb = (params->beb);
	int bek = (params->bek);
	int n1b = (params->n1b);
	int n1k = (params->n1k);
	int n2b = (params->n2b);
	int n2k = (params->n2k);
	double z_2 = (params->z_2);

	double zk1 = (alk == 1) ? Overlaps::z1(z_1) : Overlaps::z2(z_1);
	double zb1 = (alb == 1) ? Overlaps::z1(z_1) : Overlaps::z2(z_1);
	double zk2 = (bek == 1) ? Overlaps::z1(z_2) : Overlaps::z2(z_2);
	double zb2 = (beb == 1) ? Overlaps::z1(z_2) : Overlaps::z2(z_2);

	double phiket1 = Overlaps::phi(n1k, zk1);
	double phibra1 = Overlaps::phi(n1b, zb1);
	double phiket2 = Overlaps::phi(n2k, zk2);
	double phibra2 = Overlaps::phi(n2b, zb2);

	double potential = 1 / sqrt(pow(z_1 - z_2, 2) + 2);

	return phiket1 * phiket2 * potential * phibra1 * phibra2;
}

double C_integration_function_outer(double z_2, void* p)
{
	struct C_integration_params_outer* params = (struct C_integration_params_outer*)p;
	int alb = (params->alb);
	int alk = (params->alk);
	int beb = (params->beb);
	int bek = (params->bek);
	int n1b = (params->n1b);
	int n1k = (params->n1k);
	int n2b = (params->n2b);
	int n2k = (params->n2k);

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	C_integration_params_inner params_inner = { alb, alk, beb, bek, n1b, n1k, n2b, n2k, z_2 };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &C_integration_function_inner;
	F.params = &params_inner;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;

}


double Overlaps::V(const double& z)
{
    return pow(omega, 2) * pow( pow( z, 2 ) - pow( distance / 2, 2 ), 2 ) / ( 2.0 * pow(distance, 2) );
}

double Overlaps::z1(const double& z)
{
    return distance / 2.0 + z;
}

double Overlaps::z2(const double& z)
{
    return distance / 2.0 - z;
}

double Overlaps::Vho1(const double& z)
{
	using namespace std;
	return pow(omega, 2) * pow(z1(z), 2) / 2;
}


double Overlaps::Vho2(const double& z)
{
	using namespace std;
	return pow(omega, 2) * pow(z2(z), 2) / 2;
}


double Overlaps::DV1(const double& z)
{
	return V(z) - Vho1(z);
}


double Overlaps::DV2(const double& z)
{
	return V(z) - Vho2(z);
}

double Overlaps::DV(const int& alpha, const double& z)
{
	return (alpha == 1) ? DV1(z) : DV2(z);
}

double Overlaps::phi(const int& nz, const double& z)
{
	return pow(omega / M_PI1, 1.0 / 4) * exp(-omega * pow(z, 2) / 2) *
		gsl_sf_hermite(nz, sqrt(omega) * z) /
		sqrt(pow(2, nz) * gsl_sf_fact(nz));
}

int Overlaps::GetCount(const int& i, const int& j, const int& k, const int& z, const int& i_max, const int& j_max, const int& k_max, const int& z_max)
{
	return j_max * k_max * z_max * i + k_max * z_max * j + z_max * k + z;
}

int Overlaps::GetCount(const int& i, const int& j, const int& i_max, const int& j_max)
{
	return j_max * i + j;
}

double Overlaps::Eharmonic(const int& n_principal)
{
	return omega * (n_principal + 1.0 / 2);
}

double Overlaps::VIntegrate(const int& alk, const int& alb, const int& ik, const int& ib)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	integration_params params = { alk, alb, ik, ib };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &V_integration_function;
	F.params = &params;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}


double Overlaps::SIntegrate(const int& alk, const int& alb, const int& ik, const int& ib)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	integration_params params = { alk, alb, ik, ib };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &S_integration_function;
	F.params = &params;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

double Overlaps::CoulombMatel(const int& alb, const int& alk, const int& beb, const int& bek, const int& n1b, const int& n1k, const int& n2b, const int& n2k)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	C_integration_params_outer params_outer = { alb, alk, beb, bek, n1b, n1k, n2b, n2k };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &C_integration_function_outer;
	F.params = &params_outer;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}
