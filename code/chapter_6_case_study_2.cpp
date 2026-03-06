// chapter_6_case_study_2.cpp
//
// Using this code, the reliability of the continuous bridge with two spans can be calculated, given
// two load testing strategies and two load test vehicle configurations (tipper truck and semi-low
// trailer). If the axle loads are set to zero, the code can also be used to calculate the in-
// service proven strength from one year of traffic loading. The calculated reliability indices are 
// provided in Table 6.5 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/random_variable.hpp"
#include "rdv/vector_ext.hpp"
#include <Eigen/Dense>
#include <omp.h>
#include "influence_line.hpp"
#undef M_E

// Convert double to string
std::string to_string(double x) {
	std::ostringstream ss;
	ss << x;
	return ss.str();
}

// Define Eigen matrix and vector types
namespace Eigen {
	using Matrix7d = Matrix<double, 7, 7>;
	using Vector7d = Vector<double, 7>;
}

// Calculate load effect on two-span bridge given axle distances and loads
void calc_load_effect(const std::vector<double>& axle_dists, const std::vector<double>& axle_loads, 
		const std::vector<double>& x_M, const std::vector<double>& x_V, double L_1, double L_2, 
		std::vector<double>& M_E, std::vector<double>& V_E) {
	ASSERT(size(axle_dists) == size(axle_loads));
	M_E.assign(size(x_M), 0.0);
	V_E.assign(size(x_V), 0.0);
	for (size_t i = 0; i < size(axle_dists); i++) {
		for (size_t j = 0; j < size(x_M); j++)
			M_E[j] += two_spans_M(x_M[j], axle_dists[i], L_1, L_2, axle_loads[i]);
		for (size_t j = 0; j < size(x_V); j++)
			V_E[j] += two_spans_V(x_V[j], axle_dists[i], L_1, L_2, axle_loads[i]);
	}
}

// Calculate absolute values of a 2D vector
std::vector<std::vector<double>> abs(const std::vector<std::vector<double>>& x) {
	std::vector<std::vector<double>> res(size(x));
	for (size_t i = 0; i < size(x); i++) {
		res[i].resize(size(x[i]));
		for (size_t j = 0; j < size(x[i]); j++)
			res[i][j] = std::abs(x[i][j]);
	}
	return res;
}

// Perform reliability analysis for the given cross-section locations and proof load effects
void reliability_analysis(const std::vector<double>& x_M, const std::vector<double>& x_V,
		const std::vector<std::vector<double>> M_PL, const std::vector<std::vector<double>> V_PL) {	
	// Model uncertainty
	auto rv_θ_E = rdv::rv_lognormal::mv(1.0, 0.1);
	auto rv_θ_E_PL = rdv::rv_lognormal::mv(1.0, 0.1);  // plus correlation 0.7
	
	// Model uncertainty traffic load
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.1, 0.1);  // includes DAF bias
	
	//
	// Permanent load effect
	//
	
	// Random variables moment
	auto rv_M_G_11 = rdv::rv_normal::mv(405e3, 0.07);
	auto rv_M_G_12 = rdv::rv_normal::mv(453e3, 0.07);
	auto rv_M_G_13 = rdv::rv_normal::mv(356e3, 0.07);
	auto rv_M_G_2  = rdv::rv_normal::mv(810e3, 0.07);
	
	// Random variables shear
	auto rv_V_G_11 = rdv::rv_normal::mv(146e3, 0.07);
	auto rv_V_G_12 = rdv::rv_normal::mv(113e3, 0.07);
	auto rv_V_G_13 = rdv::rv_normal::mv( 81e3, 0.07);
	auto rv_V_G_21 = rdv::rv_normal::mv(243e3, 0.07);
	auto rv_V_G_22 = rdv::rv_normal::mv(275e3, 0.07);
	auto rv_V_G_23 = rdv::rv_normal::mv(308e3, 0.07);
	
	// Global correlation matrix
	Eigen::Matrix7d R_G_g = Eigen::Matrix7d::Constant(0.7);
	R_G_g.diagonal().setOnes();
	Eigen::LLT<Eigen::Matrix7d> llt_G_g(R_G_g);
	auto L_G_g = llt_G_g.matrixL();
	
	// Local correlation matrix
	Eigen::Matrix3d R_G_l = Eigen::Matrix3d::Constant(0.7);
	R_G_l.diagonal().setOnes();
	Eigen::LLT<Eigen::Matrix3d> llt_G_l(R_G_l);
	auto L_G_l = llt_G_l.matrixL();
	
	//
	// Variable traffic load effect
	//
	
	// Random variables moment
	auto rv_M_Q_11 = rdv::rv_gumbel_max::mv(789e3, 0.025);
	auto rv_M_Q_12 = rdv::rv_gumbel_max::mv(901e3, 0.025);
	auto rv_M_Q_13 = rdv::rv_gumbel_max::mv(819e3, 0.025);
	auto rv_M_Q_2  = rdv::rv_gumbel_max::mv(897e3, 0.02);
	
	// Random variables shear
	auto rv_V_Q_11 = rdv::rv_gumbel_max::mv(365e3, 0.03);
	auto rv_V_Q_12 = rdv::rv_gumbel_max::mv(329e3, 0.03);
	auto rv_V_Q_13 = rdv::rv_gumbel_max::mv(298e3, 0.03);
	auto rv_V_Q_21 = rdv::rv_gumbel_max::mv(387e3, 0.02);
	auto rv_V_Q_22 = rdv::rv_gumbel_max::mv(425e3, 0.02);
	auto rv_V_Q_23 = rdv::rv_gumbel_max::mv(463e3, 0.02);
	
	// Global correlation matrix
	Eigen::Matrix4d R_MV_Q;
	R_MV_Q << 1.00, 0.50, 0.35, 0.50,
	          0.50, 1.00, 0.50, 0.35,
	          0.35, 0.50, 1.00, 0.35, 
	          0.50, 0.35, 0.35, 1.00;
	Eigen::LLT<Eigen::Matrix4d> llt_MV_Q(R_MV_Q);
	auto L_MV_Q = llt_MV_Q.matrixL();
	
	// Local correlation matrix for moments
	Eigen::Matrix3d R_M_Q;
	R_M_Q << 1.00, 0.85, 0.85,
	         0.85, 1.00, 0.82,
	         0.85, 0.82, 1.00;
	Eigen::LLT<Eigen::Matrix3d> llt_M_Q(R_M_Q);
	auto L_M_Q = llt_M_Q.matrixL();
	
	// Local correlation matrix for shears
	Eigen::Matrix3d R_V_Q;
	R_V_Q << 1.00, 0.85, 0.82,
	         0.85, 1.00, 0.85,
		     0.82, 0.85, 1.00;
	Eigen::LLT<Eigen::Matrix3d> llt_V_Q(R_V_Q);
	auto L_V_Q = llt_V_Q.matrixL();
	
	//
	// Resistance definition
	//
	
	// (Statistical) model uncertainty
	auto rv_θ_M_R = rdv::rv_lognormal::mv(1.0, 0.05);
	auto rv_θ_V_R = rdv::rv_lognormal::mv(1.0, 0.1);
	
	// Random variables moment
	double m_m_M_R_spa = 1.5 * (rv_M_G_12.mean() + rv_M_Q_12.mean());
	double m_m_M_R_sup = 1.5 * (rv_M_G_2.mean() + rv_M_Q_2.mean());
	double m_m_M_R_sps = 1.5 * ( rv_M_G_11.mean() + rv_M_Q_11.mean() +
	                             rv_M_G_13.mean() + rv_M_Q_13.mean() ) / 2.0;
	auto rv_m_M_R_spa = rdv::rv_lognormal::mv(m_m_M_R_spa, 0.5);
	auto rv_m_M_R_sup = rdv::rv_lognormal::mv(m_m_M_R_sup, 0.5);
	auto rv_m_M_R_sps = rdv::rv_lognormal::mv(m_m_M_R_sps, 0.5);
	auto rv_V_M_R = rdv::rv_lognormal::mv(0.05, 0.1);
	
	// Random variables shear
	double m_m_V_R_end = 1.5 * (rv_V_G_11.mean() + rv_V_Q_11.mean());
	double m_m_V_R_mid = 1.5 * (rv_V_G_23.mean() + rv_V_Q_23.mean());
	auto rv_m_V_R_end = rdv::rv_lognormal::mv(m_m_V_R_end, 0.5);
	auto rv_m_V_R_mid = rdv::rv_lognormal::mv(m_m_V_R_mid, 0.5);
	auto rv_V_V_R = rdv::rv_lognormal::mv(0.1, 0.1);
	
	// Correlation global moment
	Eigen::Matrix3d R_MM_R;
	R_MM_R << 1.0, 0.4, 0.4,
	          0.4, 1.0, 0.4,
	          0.4, 0.4, 1.0;
	Eigen::LLT<Eigen::Matrix3d> llt_MM_R(R_MM_R);
	auto L_MM_R = llt_MM_R.matrixL();
	
	// Correlation local moment
	Eigen::Matrix3d R_M_R;
	R_M_R << 1.0, 0.7, 0.7,
	         0.7, 1.0, 0.9,
	         0.7, 0.9, 1.0;
	Eigen::LLT<Eigen::Matrix3d> llt_M_R(R_M_R);
	auto L_M_R = llt_M_R.matrixL();
	
	//
	// Simulation procedure
	//
	
	// Prepare random sampling
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(2e8);
	size_t post_samples = 0, failures = 0;
	std::vector<size_t> failures_M(size(x_M), 0);
	std::vector<size_t> failures_V(size(x_V), 0);
	rdv::rv_normal rv_U(0.0, 1.0);
	
	// Bayesian Monte Carlo simulation
	#pragma omp parallel reduction(+:post_samples, failures, failures_M, failures_V)
	{
		// Thread objects
		std::mt19937 gen(seeds[omp_get_thread_num()]);
		std::vector<double> M_R(size(x_M)), V_R(size(x_V));
		std::vector<double> M_G(size(x_M)), V_G(size(x_V));
		std::vector<double> M_Q(size(x_M)), V_Q(size(x_V));
		
		// Sampling procedure
		#pragma omp for
		for (size_t i_sample = 0; i_sample < samples; i_sample++) {
			// Fully correlated between cross-sections
			double u_θ_E = rv_U(gen);
			double u_θ_E_PL = 0.7 * u_θ_E + 0.714143 * rv_U(gen);
			double θ_E = rv_θ_E.from_std_norm(u_θ_E);
			double θ_E_PL = rv_θ_E_PL.from_std_norm(u_θ_E_PL);
			double C_0Q = rv_C_0Q(gen);
			double θ_M_R = rv_θ_M_R(gen);
			double m_M_R_spa = rv_m_M_R_spa(gen);
			double m_M_R_sup = rv_m_M_R_sup(gen);
			double m_M_R_sps = rv_m_M_R_sps(gen);
			double V_M_R_spa = rv_V_M_R(gen);
			double V_M_R_sup = rv_V_M_R(gen);
			double V_M_R_sps = rv_V_M_R(gen);
			double θ_V_R = rv_θ_V_R(gen);
			double m_V_R_end = rv_m_V_R_end(gen);
			double m_V_R_mid = rv_m_V_R_mid(gen);
			double V_V_R_end = rv_V_V_R(gen);
			double V_V_R_mid = rv_V_V_R(gen);
			
			// Per cross-section realisations of resistance
			{
				auto rv_M_R_spa = rdv::rv_lognormal::mv(m_M_R_spa, V_M_R_spa);
				auto rv_M_R_sup = rdv::rv_lognormal::mv(m_M_R_sup, V_M_R_sup);
				auto rv_M_R_sps = rdv::rv_lognormal::mv(m_M_R_sps, V_M_R_sps);
				Eigen::Vector3d u_R_g = L_MM_R * Eigen::Vector3d(rv_U(gen), rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_R_1 = L_M_R * Eigen::Vector3d(u_R_g[0], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_R_3 = L_M_R * Eigen::Vector3d(u_R_g[2], rv_U(gen), rv_U(gen));
				M_R[0] = rv_M_R_sps.from_std_norm(u_R_1[1]);  // M_R(2.5)
				M_R[1] = rv_M_R_spa.from_std_norm(u_R_1[0]);  // M_R(4.0)
				M_R[2] = rv_M_R_sps.from_std_norm(u_R_1[2]);  // M_R(5.5)
				M_R[3] = rv_M_R_sup.from_std_norm(u_R_g[1]);  // M_R(10.0)
				M_R[4] = rv_M_R_sps.from_std_norm(u_R_3[2]);  // M_R(14.5)
				M_R[5] = rv_M_R_spa.from_std_norm(u_R_3[0]);  // M_R(16.0)
				M_R[6] = rv_M_R_sps.from_std_norm(u_R_3[1]);  // M_R(17.5)
				
				auto rv_V_R_end = rdv::rv_lognormal::mv(m_V_R_end, V_V_R_end);
				auto rv_V_R_mid = rdv::rv_lognormal::mv(m_V_R_mid, V_V_R_mid);
				V_R[0]  = rv_V_R_end(gen);  // V_R(1.5)
				V_R[1]  = rv_V_R_end(gen);  // V_R(2.0)
				V_R[2]  = rv_V_R_end(gen);  // V_R(2.5)
				V_R[3]  = rv_V_R_mid(gen);  // V_R(7.5)
				V_R[4]  = rv_V_R_mid(gen);  // V_R(8.0)
				V_R[5]  = rv_V_R_mid(gen);  // V_R(8.5)
				V_R[6]  = rv_V_R_mid(gen);  // V_R(11.5)
				V_R[7]  = rv_V_R_mid(gen);  // V_R(12.0)
				V_R[8]  = rv_V_R_mid(gen);  // V_R(12.5)
				V_R[9]  = rv_V_R_end(gen);  // V_R(17.5)
				V_R[10] = rv_V_R_end(gen);  // V_R(18.0)
				V_R[11] = rv_V_R_end(gen);  // V_R(18.5)
			}
			
			// Per cross-section realisations of permanent load effects
			{
				Eigen::Vector7d u_G_g = L_G_g * Eigen::Vector7d(rv_U(gen), 
						rv_U(gen), rv_U(gen), rv_U(gen), rv_U(gen), rv_U(gen), rv_U(gen));
				
				Eigen::Vector3d u_G_1 = L_G_l * Eigen::Vector3d(u_G_g[0], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_G_2 = L_G_l * Eigen::Vector3d(u_G_g[2], rv_U(gen), rv_U(gen));
				M_G[0] = rv_M_G_11.from_std_norm(u_G_1[1]);
				M_G[1] = rv_M_G_12.from_std_norm(u_G_1[0]);
				M_G[2] = rv_M_G_13.from_std_norm(u_G_1[2]);
				M_G[3] = rv_M_G_2 .from_std_norm(u_G_g[1]);
				M_G[4] = rv_M_G_13.from_std_norm(u_G_2[2]);
				M_G[5] = rv_M_G_12.from_std_norm(u_G_2[0]);
				M_G[6] = rv_M_G_11.from_std_norm(u_G_2[1]);
				
				Eigen::Vector3d u_G_3 = L_G_l * Eigen::Vector3d(u_G_g[3], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_G_4 = L_G_l * Eigen::Vector3d(u_G_g[4], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_G_5 = L_G_l * Eigen::Vector3d(u_G_g[5], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_G_6 = L_G_l * Eigen::Vector3d(u_G_g[6], rv_U(gen), rv_U(gen));
				V_G[0]  = rv_V_G_11.from_std_norm(u_G_3[0]);
				V_G[1]  = rv_V_G_12.from_std_norm(u_G_3[1]);
				V_G[2]  = rv_V_G_13.from_std_norm(u_G_3[2]);
				V_G[3]  = rv_V_G_21.from_std_norm(u_G_4[2]);
				V_G[4]  = rv_V_G_22.from_std_norm(u_G_4[1]);
				V_G[5]  = rv_V_G_23.from_std_norm(u_G_4[0]);
				V_G[6]  = rv_V_G_23.from_std_norm(u_G_5[0]);
				V_G[7]  = rv_V_G_22.from_std_norm(u_G_5[1]);
				V_G[8]  = rv_V_G_21.from_std_norm(u_G_5[2]);
				V_G[9]  = rv_V_G_13.from_std_norm(u_G_6[2]);
				V_G[10] = rv_V_G_12.from_std_norm(u_G_6[1]);
				V_G[11] = rv_V_G_11.from_std_norm(u_G_6[0]);
			}
			
			//
			// Subject all cross-sections to a year of traffic load (in-service proven strength)
			//
			
			// Per cross-section realisations of variable traffic load
			{
				Eigen::Vector4d u_MV_Q = L_MV_Q * Eigen::Vector4d(
							rv_U(gen), rv_U(gen), rv_U(gen), rv_U(gen));
				
				Eigen::Vector3d u_M_Q_1 = L_M_Q * Eigen::Vector3d(u_MV_Q[0], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_M_Q_3 = L_M_Q * Eigen::Vector3d(u_MV_Q[1], rv_U(gen), rv_U(gen));
				M_Q[0] = rv_M_Q_11.from_std_norm(u_M_Q_1[1]);
				M_Q[1] = rv_M_Q_12.from_std_norm(u_M_Q_1[0]);
				M_Q[2] = rv_M_Q_13.from_std_norm(u_M_Q_1[2]);
				M_Q[3] = rv_M_Q_2(gen);
				M_Q[4] = rv_M_Q_13.from_std_norm(u_M_Q_3[2]);
				M_Q[5] = rv_M_Q_12.from_std_norm(u_M_Q_3[0]);
				M_Q[6] = rv_M_Q_11.from_std_norm(u_M_Q_3[1]);
				
				Eigen::Vector3d u_V_Q_1 = L_V_Q * Eigen::Vector3d(rv_U(gen), rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_2 = L_V_Q * Eigen::Vector3d(u_MV_Q[2], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_3 = L_V_Q * Eigen::Vector3d(u_MV_Q[3], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_4 = L_V_Q * Eigen::Vector3d(rv_U(gen), rv_U(gen), rv_U(gen));
				V_Q[0]  = rv_V_Q_11.from_std_norm(u_V_Q_1[0]);
				V_Q[1]  = rv_V_Q_12.from_std_norm(u_V_Q_1[1]);
				V_Q[2]  = rv_V_Q_13.from_std_norm(u_V_Q_1[2]);
				V_Q[3]  = rv_V_Q_21.from_std_norm(u_V_Q_2[2]);
				V_Q[4]  = rv_V_Q_22.from_std_norm(u_V_Q_2[1]);
				V_Q[5]  = rv_V_Q_23.from_std_norm(u_V_Q_2[0]);
				V_Q[6]  = rv_V_Q_23.from_std_norm(u_V_Q_3[0]);
				V_Q[7]  = rv_V_Q_22.from_std_norm(u_V_Q_3[1]);
				V_Q[8]  = rv_V_Q_21.from_std_norm(u_V_Q_3[2]);
				V_Q[9]  = rv_V_Q_13.from_std_norm(u_V_Q_4[2]);
				V_Q[10] = rv_V_Q_12.from_std_norm(u_V_Q_4[1]);
				V_Q[11] = rv_V_Q_11.from_std_norm(u_V_Q_4[0]);
			}
			
			// Evaluate LSF for all cross-sections and load effects
			bool failed = false;
			for (size_t i = 0; i < size(x_M) && !failed; i++)
				failed = θ_M_R * M_R[i] - θ_E * (M_G[i] + C_0Q * M_Q[i]) < 0.0;
			for (size_t i = 0; i < size(x_V) && !failed; i++)
				failed = θ_V_R * V_R[i] - θ_E * (V_G[i] + C_0Q * V_Q[i]) < 0.0;
			
			// Ignore this sample when the structure fails (system failure)
			if (failed)
				continue;
			
			//
			// Subject all cross-sections to load effect caused by proof load
			//
			
			// Evaluate LSF for all cross-sections and load effects
			failed = false;
			for (size_t j = 0; j < size(M_PL) && !failed; j++) {
				double c_PL = std::max(1.0 + 0.02 * rv_U(gen), 0.001);
				for (size_t i = 0; i < size(x_M) && !failed; i++)
					failed = θ_M_R * M_R[i] - (θ_E * M_G[i] + θ_E_PL * c_PL * M_PL[j][i]) < 0.0;
				for (size_t i = 0; i < size(x_V) && !failed; i++)
					failed = θ_V_R * V_R[i] - (θ_E * V_G[i] + θ_E_PL * c_PL * V_PL[j][i]) < 0.0;
			}
			
			// Ignore this sample when the structure fails (system failure)
			if (failed)
				continue;
			
			//
			// Subject all cross-sections to a year of traffic load (future)
			//
			
			// This structure has survived and matches our observation (assume system survival)
			post_samples++;
			
			// Per cross-section realisations of variable traffic load
			{
				Eigen::Vector4d u_MV_Q = L_MV_Q * Eigen::Vector4d(
						rv_U(gen), rv_U(gen), rv_U(gen), rv_U(gen));
				
				Eigen::Vector3d u_M_Q_1 = L_M_Q * Eigen::Vector3d(u_MV_Q[0], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_M_Q_3 = L_M_Q * Eigen::Vector3d(u_MV_Q[1], rv_U(gen), rv_U(gen));
				M_Q[0] = rv_M_Q_11.from_std_norm(u_M_Q_1[1]);
				M_Q[1] = rv_M_Q_12.from_std_norm(u_M_Q_1[0]);
				M_Q[2] = rv_M_Q_13.from_std_norm(u_M_Q_1[2]);
				M_Q[3] = rv_M_Q_2(gen);
				M_Q[4] = rv_M_Q_13.from_std_norm(u_M_Q_3[2]);
				M_Q[5] = rv_M_Q_12.from_std_norm(u_M_Q_3[0]);
				M_Q[6] = rv_M_Q_11.from_std_norm(u_M_Q_3[1]);
				
				Eigen::Vector3d u_V_Q_1 = L_V_Q * Eigen::Vector3d(rv_U(gen), rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_2 = L_V_Q * Eigen::Vector3d(u_MV_Q[2], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_3 = L_V_Q * Eigen::Vector3d(u_MV_Q[3], rv_U(gen), rv_U(gen));
				Eigen::Vector3d u_V_Q_4 = L_V_Q * Eigen::Vector3d(rv_U(gen), rv_U(gen), rv_U(gen));
				V_Q[0]  = rv_V_Q_11.from_std_norm(u_V_Q_1[0]);
				V_Q[1]  = rv_V_Q_12.from_std_norm(u_V_Q_1[1]);
				V_Q[2]  = rv_V_Q_13.from_std_norm(u_V_Q_1[2]);
				V_Q[3]  = rv_V_Q_21.from_std_norm(u_V_Q_2[2]);
				V_Q[4]  = rv_V_Q_22.from_std_norm(u_V_Q_2[1]);
				V_Q[5]  = rv_V_Q_23.from_std_norm(u_V_Q_2[0]);
				V_Q[6]  = rv_V_Q_23.from_std_norm(u_V_Q_3[0]);
				V_Q[7]  = rv_V_Q_22.from_std_norm(u_V_Q_3[1]);
				V_Q[8]  = rv_V_Q_21.from_std_norm(u_V_Q_3[2]);
				V_Q[9]  = rv_V_Q_13.from_std_norm(u_V_Q_4[2]);
				V_Q[10] = rv_V_Q_12.from_std_norm(u_V_Q_4[1]);
				V_Q[11] = rv_V_Q_11.from_std_norm(u_V_Q_4[0]);
			}
			
			// Evaluate LSF for all cross-sections and load effects INDIVIDUALLY
			failed = false;
			for (size_t i = 0; i < size(x_M); i++) {
				if (θ_M_R * M_R[i] - θ_E * (M_G[i] + C_0Q * M_Q[i]) < 0.0) {
					failures_M[i]++;
					failed = true;
				}
			}
			for (size_t i = 0; i < size(x_V); i++) {
				if (θ_V_R * V_R[i] - θ_E * (V_G[i] + C_0Q * V_Q[i]) < 0.0) {
					failures_V[i]++;
					failed = true;
				}
			}
			
			// Keep track of failures
			if (failed)
				failures++;
		}
	}
	
	// Display result of reliability calculation
	for (size_t i = 0; i < size(x_M); i++)
		std::cout << std::setw(15) << "M(" + to_string(x_M[i]) + ")" << 
				-rdv::normal_quantile(double(failures_M[i]) / post_samples) << "\n";
	for (size_t i = 0; i < size(x_V); i++)
		std::cout << std::setw(15) << "V(" + to_string(x_V[i]) + ")" << 
				-rdv::normal_quantile(double(failures_V[i]) / post_samples) << "\n";
	std::cout << std::setw(15) << "System:" << 
			-rdv::normal_quantile(double(failures) / post_samples) << "\n";
}

// Perform reliability analysis for proof load strategy 1: Test each critical section
void reliability_analysis_strategy_1(double L_1, double L_2,
		const std::vector<double>& x_M, const std::vector<double>& x_V) {
	std::cout << "Proof load strategy 1: Test each critical section\n";
	
	// Define axle loads
	std::vector<double> two_axles(2, 300e3);
	std::vector<double> four_axles(4, 300e3);
	
	// Calculate load effects for all proof load events (placing tandem just after 1.5 m and before 
	// 18.5 m to get maximum shear)
	std::vector<std::vector<double>> M_PL(7), V_PL(7);
	calc_load_effect({1.501, 2.7}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[0], V_PL[0]);
	calc_load_effect({3.0, 4.2}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[1], V_PL[1]);
	calc_load_effect({6.8, 8.0}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[2], V_PL[2]);
	calc_load_effect({5.8, 7.0, 13.0, 14.2}, 1.25 * four_axles, x_M, x_V, L_1, L_2, M_PL[3], V_PL[3]);
	calc_load_effect({12.0, 13.2}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[4], V_PL[4]);
	calc_load_effect({15.8, 17.0}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[5], V_PL[5]);
	calc_load_effect({17.3, 18.499}, 1.45 * two_axles, x_M, x_V, L_1, L_2, M_PL[6], V_PL[6]);
	
	// Only calculate the reliability for one particular test
	//M_PL = {M_PL[2]};
	//V_PL = {V_PL[2]};
	
	// Use absolute values of proof load effect
	M_PL = abs(M_PL);
	V_PL = abs(V_PL);
	
	// Perform reliability analysis
	reliability_analysis(x_M, x_V, M_PL, V_PL);
}

// Perform reliability analysis for proof load strategy 2: Optimised locations
void reliability_analysis_strategy_2(double L_1, double L_2,
		const std::vector<double>& x_M, const std::vector<double>& x_V) {
	std::cout << "Proof load strategy 2: Optimised locations\n";
	
	// Define axle loads
	std::vector<double> two_axles(2, 300e3);
	std::vector<double> four_axles(4, 300e3);
	
	// Calculate load effects for all proof load events
	std::vector<std::vector<double>> M_PL(3), V_PL(3);
	calc_load_effect({ 2.5,  3.7}, 1.55 * two_axles, x_M, x_V, L_1, L_2, M_PL[0], V_PL[0]);
	calc_load_effect({ 6.3,  7.5, 12.5, 13.7}, 1.35 * four_axles, x_M, x_V, L_1, L_2, M_PL[1], V_PL[1]);
	calc_load_effect({16.3, 17.5}, 1.55 * two_axles, x_M, x_V, L_1, L_2, M_PL[2], V_PL[2]);
	
	// Only calculate the reliability for one particular test
	//M_PL = {M_PL[2]};
	//V_PL = {V_PL[2]};
	
	// Use absolute values of proof load effect
	M_PL = abs(M_PL);
	V_PL = abs(V_PL);
	
	// Perform reliability analysis
	reliability_analysis(x_M, x_V, M_PL, V_PL);
}

// Perform reliability analysis for a moving vehicle on a two-span bridge
void reliability_analysis_vehicle(double L_1, double L_2, 
		const std::vector<double>& x_M, const std::vector<double>& x_V, 
		const std::vector<double>& axle_dists, const std::vector<double>& axle_loads) {
	// Just one load event; the envelope
	std::vector<std::vector<double>> M_PL = {std::vector<double>(size(x_M), 0.0)};
	std::vector<std::vector<double>> V_PL = {std::vector<double>(size(x_V), 0.0)};
	for (double x = 0.0; x < L_1 + L_2 + sum(axle_dists); x += 0.001) {
		// Calculate load effect of truck at this position
		std::vector<double> M(size(x_M), 0.0);
		std::vector<double> V(size(x_V), 0.0);
		double x_axle = x;
		for (size_t i = 0; i < size(axle_dists); i++) {
			x_axle -= axle_dists[i];
			if (x_axle > 0.0 && x_axle < L_1 + L_2) {
				for (size_t j = 0; j < size(x_M); j++)
					M[j] += two_spans_M(x_M[j], x_axle, L_1, L_2, axle_loads[i]);
				for (size_t j = 0; j < size(x_V); j++)
					V[j] += two_spans_V(x_V[j], x_axle, L_1, L_2, axle_loads[i]);
			}
		}
		
		// Save if maximum
		for (size_t i = 0; i < size(x_M); i++)
			M_PL[0][i] = std::max(M_PL[0][i], (i == 3 ? -1 : 1) * M[i]);
		for (size_t i = 0; i < size(x_V); i++)
			V_PL[0][i] = std::max(V_PL[0][i], abs(V[i]));
	}
	
	// Perform reliability analysis
	reliability_analysis(x_M, x_V, M_PL, V_PL);
}

// Reliability analysis for vehicle type 1 (tipper truck)
void reliability_analysis_vehicle_1(double L_1, double L_2,
		const std::vector<double>& x_M, const std::vector<double>& x_V) {
	// Moving vehicle definition
	double P = 1.32 * 300e3;
	std::cout << "Vehicle type 1: Tipper truck, axle load = " << P / 9810 << " t\n";
	std::vector<double> axle_dists = {0.0, 3.2, 1.3, 6.0, 1.3};
	std::vector<double> axle_loads = {78e3, P, P, P, P};
	
	// Perform reliability analysis
	reliability_analysis_vehicle(L_1, L_2, x_M, x_V, axle_dists, axle_loads);
}

// Reliability analysis for vehicle type 2 (semi-low trailer)
void reliability_analysis_vehicle_2(double L_1, double L_2,
		const std::vector<double>& x_M, const std::vector<double>& x_V) {
	// Moving vehicle definition
	double P = 0.81 * 300e3;
	std::cout << "Vehicle type 2: Semi-low trailer, axle load = " << P / 9810 << " t\n";
	std::vector<double> axle_dists = { 0.0, 3.2, 1.3, 3.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
	std::vector<double> axle_loads = {78e3, 1e3, 1e3,   P,   P,   P,   P,   P,   P,   P,   P};
	
	// Perform reliability analysis
	reliability_analysis_vehicle(L_1, L_2, x_M, x_V, axle_dists, axle_loads);
}

// Program entry point
int main(int argc, char *argv[]) {
	// Initialisation
	rdv::tic();
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	std::cout.setf(std::ios::left);
	std::cout.precision(3);
	
	// Bridge properties
	double L_1 = 10.0;
	double L_2 = 10.0;
	
	// Cross-section locations
	std::vector<double> x_M = {2.5, 4.0, 5.5, 10.0, 14.5, 16.0, 17.5};
	std::vector<double> x_V = {1.5, 2.0, 2.5, 7.5, 8.0, 8.5, 11.5, 12.0, 12.5, 17.5, 18.0, 18.5};
	
	// Perform reliability analyses for two proof load strategies
	reliability_analysis_strategy_1(L_1, L_2, x_M, x_V);
	reliability_analysis_strategy_2(L_1, L_2, x_M, x_V);
	
	// Perform reliability analyses for two vehicle types
	reliability_analysis_vehicle_1(L_1, L_2, x_M, x_V);
	reliability_analysis_vehicle_2(L_1, L_2, x_M, x_V);
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
