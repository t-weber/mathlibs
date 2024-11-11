/**
 * container-agnostic math algorithms
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 2017-2021
 * @license see 'LICENSE' file
 *
 * @see references for algorithms:
 * 	- (Arens15): T. Arens et al., ISBN: 978-3-642-44919-2, DOI: 10.1007/978-3-642-44919-2 (2015).
 * 	- (Arfken13): G. B. Arfken et al., ISBN: 978-0-12-384654-9, DOI: 10.1016/C2009-0-30629-7 (2013).
 * 	- (DesktopBronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its HTML version "Desktop Bronstein"].
 * 	- (Bronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its paperback version].
 * 	- (Merziger06): G. Merziger and T. Wirth, ISBN: 3923923333 (2006).
 * 	- (Scarpino11): M. Scarpino, ISBN: 978-1-6172-9017-6 (2011).
 * 	- (Shirane02): G. Shirane et al., ISBN: 978-0-5214-1126-4 (2002).
 * 	- (Kuipers02): J. B. Kuipers, ISBN: 0-691-05872-5 (2002).
 * 	- (FUH 2021): A. Schulz and J. Rollin, "Effiziente Algorithmen", Kurs 1684 (2021), Fernuni Hagen (https://vu.fernuni-hagen.de/lvuweb/lvu/app/Kurs/01684).
 * 	- (Sellers 2014): G. Sellers et al., ISBN: 978-0-321-90294-8 (2014).
 */

#ifndef __PHYS_ALGOS_H__
#define __PHYS_ALGOS_H__

#include <functional>
#include <limits>

#include "matrix_algos.h"



// math
namespace m {

/**
 * general structure factor calculation
 * e.g. type T as vector (complex number) for magnetic (nuclear) structure factor
 * Ms_or_bs:
	- nuclear scattering lengths for nuclear neutron scattering or
	- atomic form factors for x-ray scattering
	- magnetisation (* magnetic form factor) for magnetic neutron scattering
 * Rs: atomic positions
 * Q: scattering vector G for nuclear scattering or G+k for magnetic scattering with propagation vector k
 * fs: optional magnetic form factors
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50002-1
 * @see (Shirane02), p. 25, equ. 2.26 for nuclear structure factor
 * @see (Shirane02), p. 40, equ. 2.81 for magnetic structure factor
 */
template<class t_vec, class T = t_vec, template<class...> class t_cont = std::vector,
	class t_cplx = std::complex<double>>
T structure_factor(const t_cont<T>& Ms_or_bs, const t_cont<t_vec>& Rs, const t_vec& Q, const t_vec* fs=nullptr)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx cI(0,1);
	constexpr t_real twopi = pi<t_real> * t_real(2);
	constexpr t_real expsign = -1;

	T F;
	if constexpr(is_vec<T>)
		F = zero<T>(Rs.begin()->size());	// always 3 dims...
	else if constexpr(is_complex<T>)
		F = T(0);

	auto iterM_or_b = Ms_or_bs.begin();
	auto iterR = Rs.begin();
	typename t_vec::const_iterator iterf;
	if(fs) iterf = fs->begin();

	while(iterM_or_b != Ms_or_bs.end() && iterR != Rs.end())
	{
		// if form factors are given, use them, otherwise set to 1
		t_real f = t_real(1);
		if(fs)
		{
			auto fval = *iterf;
			if constexpr(is_complex<decltype(fval)>)
				f = fval.real();
			else
				f = fval;
		}

		// structure factor
		F += (*iterM_or_b) * (f * std::exp(expsign * cI * twopi * inner<t_vec>(Q, *iterR)));

		// next M or b if available (otherwise keep current)
		auto iterM_or_b_next = std::next(iterM_or_b, 1);
		if(iterM_or_b_next != Ms_or_bs.end())
			iterM_or_b = iterM_or_b_next;

		if(fs)
		{
			// next form factor if available (otherwise keep current)
			auto iterf_next = std::next(iterf, 1);
			if(iterf_next != fs->end())
				iterf = iterf_next;
		}

		// next atomic position
		std::advance(iterR, 1);
	}

	return F;
}



/**
 * Blume-Maleev equation
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - pp. 225-226
 * @see lecture notes by P. J. Brown, 2006
 */
template<class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec> blume_maleev(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_vec<t_vec>
{
	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);
	constexpr t_cplx imag(0, 1);

	// ------------------------------------------------------------------------
	// intensity
	// nuclear
	t_cplx I = N*NConj;

	// nuclear-magnetic
	I += NConj*inner<t_vec>(P_i, Mperp);
	I += N*inner<t_vec>(Mperp, P_i);

	// magnetic, non-chiral
	I += inner<t_vec>(Mperp, Mperp);

	// magnetic, chiral
	I += -imag * inner<t_vec>(P_i, cross<t_vec>({ Mperp, MperpConj }));
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// polarisation vector
	// nuclear
	t_vec P_f = P_i * N*NConj;

	// nuclear-magnetic
	P_f += NConj * Mperp;
	P_f += N * MperpConj;
	P_f += imag * N * cross<t_vec>({ P_i, MperpConj });
	P_f += -imag * NConj * cross<t_vec>({ P_i, Mperp });

	// magnetic, non-chiral
	P_f += Mperp * inner<t_vec>(Mperp, P_i);
	P_f += MperpConj * inner<t_vec>(P_i, Mperp);
	P_f += -P_i * inner<t_vec>(Mperp, Mperp);

	// magnetic, chiral
	P_f += imag * cross<t_vec>({ Mperp, MperpConj });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev equation
 * calculated indirectly with density matrix
 *
 * V   = N*1 + <Mperp|sigma>
 * I   = 0.5 * tr( V^H V rho )
 * P_f = 0.5 * tr( V^H sigma V rho ) / I
 *
 * @returns scattering intensity and final polarisation vector
 *
 * @see lecture notes by P. J. Brown, 2006
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - pp. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec> blume_maleev_indir(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	// spin-1/2
	constexpr t_cplx c = 0.5;

	// vector of pauli matrices
	const auto sigma = su2_matrices<std::vector<t_mat>>(false);

	// density matrix
	const auto density = bloch_density_op<t_vec, t_mat>(P_i, c);

	// potential
	const auto V_mag = proj_su2<t_vec, t_mat>(Mperp, true);
	const auto V_nuc = N * unit<t_mat>(2);
	const auto V = V_nuc + V_mag;
	const auto VConj = herm(V);

	// scattering intensity
	t_cplx I = c * trace(VConj*V * density/c);

	// ------------------------------------------------------------------------
	// scattered polarisation vector
	const auto m0 = (VConj * sigma[0]) * V * density/c;
	const auto m1 = (VConj * sigma[1]) * V * density/c;
	const auto m2 = (VConj * sigma[2]) * V * density/c;

	t_vec P_f = create<t_vec>({ c*trace(m0), c*trace(m1), c*trace(m2) });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}



/**
 * calculates the berry connection
 * @see equ. 7 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type>
std::vector<t_vec> berry_connection(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat> && is_vec<t_vec> && is_vec<t_vec_real> && is_complex<t_cplx>
{
	using t_size = decltype(Q.size());
	constexpr const t_cplx imag{0, 1};

	const t_mat evecs = get_evecs(Q);
	constexpr const t_size BANDS = evecs.size1();
	/*constexpr*/ const t_size DIM = Q.size();

	// to ensure correct commutators
	t_mat comm = unit<t_mat>(BANDS);

	std::vector<t_vec> connections{};
	connections.reserve(BANDS);
	for(t_size band = 0; band < BANDS; ++band)
	{
		connections.emplace_back(create<t_vec>(DIM));
		if(band >= BANDS / 2)
			comm(band, band) = -1;
	}

	for(t_size dim = 0; dim < DIM; ++dim)
	{
		t_vec_real Q1 = Q;
		Q1[dim] += delta;

		// differentiate eigenvector matrix
		t_mat evecs_diff = (get_evecs(Q1) - evecs) / delta;

		t_mat evecs_H = m::herm(evecs);
		t_mat C = comm * evecs_H * comm * evecs_diff;

		for(t_size band = 0; band < BANDS; ++band)
			connections[band][dim] = C(band, band) * imag;
	}

	return connections;
}



/**
 * calculates the 2d berry curvature
 * @see equ. 8 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type>
std::vector<t_cplx> berry_curvature_2d(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon(),
	decltype(t_vec{}.size()) dim1 = 0, decltype(t_vec{}.size()) dim2 = 1)
requires is_mat<t_mat> && is_vec<t_vec> && is_vec<t_vec_real> && is_complex<t_cplx>
{
	using t_size = decltype(Q.size());

	const t_mat evecs = get_evecs(Q);
	constexpr const t_size BANDS = evecs.size1();
	/*constexpr*/ const t_size DIM = Q.size();
	assert(DIM == 3);

	t_vec_real h = Q, k = Q;
	h[dim1] += delta;
	k[dim2] += delta;

	std::vector<t_vec> connections =
		berry_connection<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, Q, delta);
	std::vector<t_vec> connections_h =
		berry_connection<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, h, delta);
	std::vector<t_vec> connections_k =
		berry_connection<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, k, delta);

	std::vector<t_cplx> curvatures{};
	curvatures.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
	{
		// differentiate connection's y component with respect to h
		t_cplx curv1 = (connections_h[band][dim2] - connections[band][dim2]) / delta;
		// differentiate connection's x component with respect to k
		t_cplx curv2 = (connections_k[band][dim1] - connections[band][dim1]) / delta;

		curvatures.emplace_back(curv1 - curv2);
	}

	return curvatures;
}



/**
 * calculates the full berry curvature tensor
 * @see equ. 8 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat_band, class t_mat_dim, class t_vec_dim, class t_vec_dim_real,
	typename t_cplx = typename t_vec_dim::value_type,
	typename t_real = typename t_cplx::value_type>
std::vector<t_mat_dim> berry_curvature(
	const std::function<t_mat_band(const t_vec_dim_real& Q)>& get_evecs,
	const t_vec_dim_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat_band> && is_mat<t_mat_dim> &&
	is_vec<t_vec_dim> && is_vec<t_vec_dim_real> && is_complex<t_cplx>
{
	using t_size = decltype(Q.size());

	const t_mat_band evecs = get_evecs(Q);
	constexpr const t_size BANDS = evecs.size1();
	/*constexpr*/ const t_size DIM = Q.size();

	// berry connections at Q + [0, ..., 0, delta, 0, ..., 0]
	std::vector<std::vector<t_vec_dim>> connections_plus;
	connections_plus.reserve(DIM);
	for(t_size dim = 0; dim < DIM; ++dim)
	{
		t_vec_dim_real Q_plus = Q;
		Q_plus[dim] += delta;

		std::vector<t_vec_dim> connection =
			berry_connection<t_mat_band, t_vec_dim, t_vec_dim_real, t_cplx, t_real>(
				get_evecs, Q_plus, delta);

		/*for(const t_vec_dim& vec : connection)
		{
			using namespace m_ops;
			std::cout << "dim = " << dim << ": " << vec << std::endl;
		}*/

		connections_plus.emplace_back(std::move(connection));
	}

	// berry connection at Q
	std::vector<t_vec_dim> connections =
		berry_connection<t_mat_band, t_vec_dim, t_vec_dim_real, t_cplx, t_real>(
			get_evecs, Q, delta);

	std::vector<t_mat_dim> curvatures{};
	curvatures.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
	{
		t_mat_dim curvature = create<t_mat_dim>(DIM, DIM);

		for(t_size dim1 = 0; dim1 < DIM; ++dim1)
		{
			const std::vector<t_vec_dim>& connections_plus1 = connections_plus[dim1];

			for(t_size dim2 = 0; dim2 < DIM; ++dim2)
			{
				const std::vector<t_vec_dim>& connections_plus2 = connections_plus[dim2];

				// differentiate connection's y component with respect to h
				t_cplx curv1 = (connections_plus1[band][dim2] - connections[band][dim2]) / delta;

				// differentiate connection's x component with respect to k
				t_cplx curv2 = (connections_plus2[band][dim1] - connections[band][dim1]) / delta;

				// curvature tensor element
				curvature(dim1, dim2) = curv1 - curv2;
			}
		}

		curvatures.emplace_back(std::move(curvature));
	}

	return curvatures;
}


}
#endif

