#pragma once
class Overlaps
{
public:
	static double V(const double& z, const double & distance);
	static double z1(const double& z, const double & distance);
	static double z2(const double& z, const double & distance);
	static double Vho1(const double& z, const double & distance);
	static double Vho2(const double& z, const double & distance);
	static double DV1(const double& z, const double & distance);
	static double DV2(const double& z, const double & distance);
	static double DV(const int& alpha, const double& z, const double & distance);

	static double phi(const int& nzk, const double& z_ket);

	static int GetCount(const int& i, const int& j, const int& k, const int& z, const int& i_max, const int& j_max, const int& k_max, const int& z_max);
	static int GetCount(const int& i, const int& j, const int& i_max, const int& j_max);

	static double Eharmonic(const int& n_principal);

	static double VIntegrate(const int& alk, const int& alb, const int& ik, const int& ib, const double & distance);
	static double SIntegrate(const int& alk, const int& alb, const int& ik, const int& ib, const double & distance);

	static double CoulombMatel(const int& alb, const int& alk, const int& beb, const int& bek, const int& n1b, const int& n1k, const int& n2b, const int& n2k, const double & distance );
};

