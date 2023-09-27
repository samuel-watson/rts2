// Generated by rstantools.  Do not edit by hand.

/*
    rts2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rts2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rts2.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_mcml_poisson_region_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'string', line 36, column 2 to column 21)",
                                                      " (in 'string', line 39, column 9 to column 13)",
                                                      " (in 'string', line 39, column 2 to column 38)",
                                                      " (in 'string', line 40, column 9 to column 19)",
                                                      " (in 'string', line 40, column 2 to column 81)",
                                                      " (in 'string', line 41, column 2 to column 23)",
                                                      " (in 'string', line 42, column 2 to column 18)",
                                                      " (in 'string', line 23, column 2 to column 8)",
                                                      " (in 'string', line 24, column 2 to column 9)",
                                                      " (in 'string', line 25, column 2 to column 14)",
                                                      " (in 'string', line 26, column 2 to column 10)",
                                                      " (in 'string', line 27, column 9 to column 19)",
                                                      " (in 'string', line 27, column 2 to column 24)",
                                                      " (in 'string', line 28, column 9 to column 13)",
                                                      " (in 'string', line 28, column 2 to column 23)",
                                                      " (in 'string', line 29, column 9 to column 13)",
                                                      " (in 'string', line 29, column 14 to column 18)",
                                                      " (in 'string', line 29, column 2 to column 23)",
                                                      " (in 'string', line 30, column 8 to column 18)",
                                                      " (in 'string', line 30, column 2 to column 20)",
                                                      " (in 'string', line 31, column 22 to column 31)",
                                                      " (in 'string', line 31, column 2 to column 33)",
                                                      " (in 'string', line 32, column 23 to column 26)",
                                                      " (in 'string', line 32, column 2 to column 28)",
                                                      " (in 'string', line 33, column 9 to column 12)",
                                                      " (in 'string', line 33, column 2 to column 24)",
                                                      " (in 'string', line 36, column 9 to column 13)",
                                                      " (in 'string', line 5, column 11 to column 16)",
                                                      " (in 'string', line 5, column 4 to column 25)",
                                                      " (in 'string', line 6, column 11 to column 13)",
                                                      " (in 'string', line 6, column 4 to column 17)",
                                                      " (in 'string', line 8, column 6 to column 37)",
                                                      " (in 'string', line 10, column 8 to column 44)",
                                                      " (in 'string', line 11, column 16 to column 21)",
                                                      " (in 'string', line 11, column 8 to column 23)",
                                                      " (in 'string', line 12, column 15 to column 20)",
                                                      " (in 'string', line 12, column 8 to column 46)",
                                                      " (in 'string', line 14, column 10 to column 57)",
                                                      " (in 'string', line 13, column 25 to line 15, column 9)",
                                                      " (in 'string', line 13, column 8 to line 15, column 9)",
                                                      " (in 'string', line 16, column 8 to column 95)",
                                                      " (in 'string', line 9, column 20 to line 17, column 7)",
                                                      " (in 'string', line 9, column 6 to line 17, column 7)",
                                                      " (in 'string', line 7, column 18 to line 18, column 5)",
                                                      " (in 'string', line 7, column 4 to line 18, column 5)",
                                                      " (in 'string', line 19, column 4 to column 18)",
                                                      " (in 'string', line 4, column 29 to line 20, column 3)"};
template <typename T3__, typename T4__, typename T5__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T3__>, stan::value_type_t<T4__>,
stan::value_type_t<T5__>>, -1, 1>
gen_lambda(const int& nT, const int& nR, const int& nS, const T3__& xb_arg__,
           const T4__& w_arg__, const T5__& gamma_arg__,
           const std::vector<int>& cell_id, const std::vector<int>& n_cell,
           std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T3__>,
          stan::value_type_t<T4__>,
          stan::value_type_t<T5__>>;
  const auto& xb = to_ref(xb_arg__);
  const auto& w = to_ref(w_arg__);
  const auto& gamma = to_ref(gamma_arg__);
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    current_statement__ = 28;
    validate_non_negative_index("lambda", "nT * nR", (nT * nR));
    Eigen::Matrix<local_scalar_t__, -1, 1> lambda;
    lambda = Eigen::Matrix<local_scalar_t__, -1, 1>((nT * nR));
    stan::math::fill(lambda, DUMMY_VAR__);
    
    current_statement__ = 30;
    validate_non_negative_index("u", "nS", nS);
    Eigen::Matrix<local_scalar_t__, -1, 1> u;
    u = Eigen::Matrix<local_scalar_t__, -1, 1>(nS);
    stan::math::fill(u, DUMMY_VAR__);
    
    current_statement__ = 45;
    for (int t = 1; t <= nT; ++t) {
      current_statement__ = 32;
      assign(u, nil_index_list(),
        rvalue(gamma,
          cons_list(index_min_max((((t - 1) * nS) + 1), (t * nS)),
            nil_index_list()), "gamma"), "assigning variable u");
      current_statement__ = 43;
      for (int r = 1; r <= nR; ++r) {
        int lsize;
        lsize = std::numeric_limits<int>::min();
        
        current_statement__ = 33;
        lsize = (n_cell[((r + 1) - 1)] - n_cell[(r - 1)]);
        current_statement__ = 34;
        validate_non_negative_index("idx", "lsize", lsize);
        std::vector<int> idx;
        idx = std::vector<int>(lsize, std::numeric_limits<int>::min());
        
        current_statement__ = 36;
        validate_non_negative_index("f", "lsize", lsize);
        Eigen::Matrix<local_scalar_t__, -1, 1> f;
        f = Eigen::Matrix<local_scalar_t__, -1, 1>(lsize);
        stan::math::fill(f, DUMMY_VAR__);
        
        current_statement__ = 37;
        assign(f, nil_index_list(), rep_vector(0, lsize),
          "assigning variable f");
        current_statement__ = 40;
        for (int l = 1; l <= lsize; ++l) {
          current_statement__ = 38;
          assign(f, cons_list(index_uni(l), nil_index_list()),
            gamma[((cell_id[(((n_cell[(r - 1)] + l) - 1) - 1)] +
                     ((t - 1) * nS)) - 1)], "assigning variable f");}
        current_statement__ = 41;
        assign(lambda,
          cons_list(index_uni((((t - 1) * nR) + r)), nil_index_list()),
          (stan::math::exp(xb[((r + ((t - 1) * nR)) - 1)]) *
            multiply(
              transpose(
                rvalue(w,
                  cons_list(
                    index_min_max(n_cell[(r - 1)], (n_cell[((r + 1) - 1)] -
                                                     1)), nil_index_list()),
                  "w")), stan::math::exp(f))), "assigning variable lambda");}
    }
    current_statement__ = 46;
    return lambda;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct gen_lambda_functor__ {
template <typename T3__, typename T4__, typename T5__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T3__>, stan::value_type_t<T4__>,
stan::value_type_t<T5__>>, -1, 1>
operator()(const int& nT, const int& nR, const int& nS, const T3__& xb,
           const T4__& w, const T5__& gamma, const std::vector<int>& cell_id,
           const std::vector<int>& n_cell, std::ostream* pstream__)  const 
{
return gen_lambda(nT, nR, nS, xb, w, gamma, cell_id, n_cell, pstream__);
}
};
#include <stan_meta_header.hpp>
class model_mcml_poisson_region final : public model_base_crtp<model_mcml_poisson_region> {
private:
  int N;
  int nT;
  int nRegion;
  int n_Q;
  Eigen::Matrix<double, -1, 1> Xb;
  Eigen::Matrix<double, -1, 1> Xb_cell;
  Eigen::Matrix<double, -1, -1> ZL;
  std::vector<int> y;
  std::vector<int> n_cell;
  std::vector<int> cell_id;
  Eigen::Matrix<double, -1, 1> q_weights;
  int gamma_1dim__;
 
public:
  ~model_mcml_poisson_region() { }
  
  inline std::string model_name() const final { return "model_mcml_poisson_region"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_mcml_poisson_region(stan::io::var_context& context__,
                            unsigned int random_seed__ = 0,
                            std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_mcml_poisson_region_namespace::model_mcml_poisson_region";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 8;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 8;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 9;
      context__.validate_dims("data initialization","nT","int",
          context__.to_vec());
      nT = std::numeric_limits<int>::min();
      
      current_statement__ = 9;
      nT = context__.vals_i("nT")[(1 - 1)];
      current_statement__ = 10;
      context__.validate_dims("data initialization","nRegion","int",
          context__.to_vec());
      nRegion = std::numeric_limits<int>::min();
      
      current_statement__ = 10;
      nRegion = context__.vals_i("nRegion")[(1 - 1)];
      current_statement__ = 11;
      context__.validate_dims("data initialization","n_Q","int",
          context__.to_vec());
      n_Q = std::numeric_limits<int>::min();
      
      current_statement__ = 11;
      n_Q = context__.vals_i("n_Q")[(1 - 1)];
      current_statement__ = 12;
      validate_non_negative_index("Xb", "nRegion * nT", (nRegion * nT));
      current_statement__ = 13;
      context__.validate_dims("data initialization","Xb","double",
          context__.to_vec((nRegion * nT)));
      Xb = Eigen::Matrix<double, -1, 1>((nRegion * nT));
      stan::math::fill(Xb, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> Xb_flat__;
        current_statement__ = 13;
        assign(Xb_flat__, nil_index_list(), context__.vals_r("Xb"),
          "assigning variable Xb_flat__");
        current_statement__ = 13;
        pos__ = 1;
        current_statement__ = 13;
        for (int sym1__ = 1; sym1__ <= (nRegion * nT); ++sym1__) {
          current_statement__ = 13;
          assign(Xb, cons_list(index_uni(sym1__), nil_index_list()),
            Xb_flat__[(pos__ - 1)], "assigning variable Xb");
          current_statement__ = 13;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 14;
      validate_non_negative_index("Xb_cell", "N * nT", (N * nT));
      current_statement__ = 15;
      context__.validate_dims("data initialization","Xb_cell","double",
          context__.to_vec((N * nT)));
      Xb_cell = Eigen::Matrix<double, -1, 1>((N * nT));
      stan::math::fill(Xb_cell, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> Xb_cell_flat__;
        current_statement__ = 15;
        assign(Xb_cell_flat__, nil_index_list(), context__.vals_r("Xb_cell"),
          "assigning variable Xb_cell_flat__");
        current_statement__ = 15;
        pos__ = 1;
        current_statement__ = 15;
        for (int sym1__ = 1; sym1__ <= (N * nT); ++sym1__) {
          current_statement__ = 15;
          assign(Xb_cell, cons_list(index_uni(sym1__), nil_index_list()),
            Xb_cell_flat__[(pos__ - 1)], "assigning variable Xb_cell");
          current_statement__ = 15;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 16;
      validate_non_negative_index("ZL", "N * nT", (N * nT));
      current_statement__ = 17;
      validate_non_negative_index("ZL", "N * nT", (N * nT));
      current_statement__ = 18;
      context__.validate_dims("data initialization","ZL","double",
          context__.to_vec((N * nT), (N * nT)));
      ZL = Eigen::Matrix<double, -1, -1>((N * nT), (N * nT));
      stan::math::fill(ZL, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> ZL_flat__;
        current_statement__ = 18;
        assign(ZL_flat__, nil_index_list(), context__.vals_r("ZL"),
          "assigning variable ZL_flat__");
        current_statement__ = 18;
        pos__ = 1;
        current_statement__ = 18;
        for (int sym1__ = 1; sym1__ <= (N * nT); ++sym1__) {
          current_statement__ = 18;
          for (int sym2__ = 1; sym2__ <= (N * nT); ++sym2__) {
            current_statement__ = 18;
            assign(ZL,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              ZL_flat__[(pos__ - 1)], "assigning variable ZL");
            current_statement__ = 18;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 19;
      validate_non_negative_index("y", "nRegion * nT", (nRegion * nT));
      current_statement__ = 20;
      context__.validate_dims("data initialization","y","int",
          context__.to_vec((nRegion * nT)));
      y = std::vector<int>((nRegion * nT), std::numeric_limits<int>::min());
      
      current_statement__ = 20;
      assign(y, nil_index_list(), context__.vals_i("y"),
        "assigning variable y");
      current_statement__ = 21;
      validate_non_negative_index("n_cell", "nRegion + 1", (nRegion + 1));
      current_statement__ = 22;
      context__.validate_dims("data initialization","n_cell","int",
          context__.to_vec((nRegion + 1)));
      n_cell = std::vector<int>((nRegion + 1), std::numeric_limits<int>::min());
      
      current_statement__ = 22;
      assign(n_cell, nil_index_list(), context__.vals_i("n_cell"),
        "assigning variable n_cell");
      current_statement__ = 22;
      for (int sym1__ = 1; sym1__ <= (nRegion + 1); ++sym1__) {
        current_statement__ = 22;
        current_statement__ = 22;
        check_greater_or_equal(function__, "n_cell[sym1__]",
                               n_cell[(sym1__ - 1)], 1);}
      current_statement__ = 23;
      validate_non_negative_index("cell_id", "n_Q", n_Q);
      current_statement__ = 24;
      context__.validate_dims("data initialization","cell_id","int",
          context__.to_vec(n_Q));
      cell_id = std::vector<int>(n_Q, std::numeric_limits<int>::min());
      
      current_statement__ = 24;
      assign(cell_id, nil_index_list(), context__.vals_i("cell_id"),
        "assigning variable cell_id");
      current_statement__ = 24;
      for (int sym1__ = 1; sym1__ <= n_Q; ++sym1__) {
        current_statement__ = 24;
        current_statement__ = 24;
        check_greater_or_equal(function__, "cell_id[sym1__]",
                               cell_id[(sym1__ - 1)], 1);}
      current_statement__ = 25;
      validate_non_negative_index("q_weights", "n_Q", n_Q);
      current_statement__ = 26;
      context__.validate_dims("data initialization","q_weights","double",
          context__.to_vec(n_Q));
      q_weights = Eigen::Matrix<double, -1, 1>(n_Q);
      stan::math::fill(q_weights, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> q_weights_flat__;
        current_statement__ = 26;
        assign(q_weights_flat__, nil_index_list(),
          context__.vals_r("q_weights"),
          "assigning variable q_weights_flat__");
        current_statement__ = 26;
        pos__ = 1;
        current_statement__ = 26;
        for (int sym1__ = 1; sym1__ <= n_Q; ++sym1__) {
          current_statement__ = 26;
          assign(q_weights, cons_list(index_uni(sym1__), nil_index_list()),
            q_weights_flat__[(pos__ - 1)], "assigning variable q_weights");
          current_statement__ = 26;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 27;
      gamma_1dim__ = std::numeric_limits<int>::min();
      
      current_statement__ = 27;
      gamma_1dim__ = (N * nT);
      current_statement__ = 27;
      validate_non_negative_index("gamma", "N * nT", gamma_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += gamma_1dim__;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_mcml_poisson_region_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      Eigen::Matrix<local_scalar_t__, -1, 1> gamma;
      gamma = Eigen::Matrix<local_scalar_t__, -1, 1>(gamma_1dim__);
      stan::math::fill(gamma, DUMMY_VAR__);
      
      current_statement__ = 1;
      gamma = in__.vector(gamma_1dim__);
      {
        current_statement__ = 2;
        validate_non_negative_index("u", "N * nT", (N * nT));
        Eigen::Matrix<local_scalar_t__, -1, 1> u;
        u = Eigen::Matrix<local_scalar_t__, -1, 1>((N * nT));
        stan::math::fill(u, DUMMY_VAR__);
        
        current_statement__ = 3;
        assign(u, nil_index_list(), add(Xb_cell, multiply(ZL, gamma)),
          "assigning variable u");
        current_statement__ = 4;
        validate_non_negative_index("mu", "nRegion * nT", (nRegion * nT));
        Eigen::Matrix<local_scalar_t__, -1, 1> mu;
        mu = Eigen::Matrix<local_scalar_t__, -1, 1>((nRegion * nT));
        stan::math::fill(mu, DUMMY_VAR__);
        
        current_statement__ = 5;
        assign(mu, nil_index_list(),
          gen_lambda(nT, nRegion, N, Xb, q_weights, u, cell_id,
            n_cell, pstream__), "assigning variable mu");
        current_statement__ = 6;
        lp_accum__.add(std_normal_lpdf<propto__>(gamma));
        current_statement__ = 7;
        lp_accum__.add(poisson_lpmf<propto__>(y, mu));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_mcml_poisson_region_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      Eigen::Matrix<double, -1, 1> gamma;
      gamma = Eigen::Matrix<double, -1, 1>(gamma_1dim__);
      stan::math::fill(gamma, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      gamma = in__.vector(gamma_1dim__);
      for (int sym1__ = 1; sym1__ <= gamma_1dim__; ++sym1__) {
        vars__.emplace_back(gamma[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      Eigen::Matrix<double, -1, 1> gamma;
      gamma = Eigen::Matrix<double, -1, 1>(gamma_1dim__);
      stan::math::fill(gamma, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> gamma_flat__;
        current_statement__ = 1;
        assign(gamma_flat__, nil_index_list(), context__.vals_r("gamma"),
          "assigning variable gamma_flat__");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= gamma_1dim__; ++sym1__) {
          current_statement__ = 1;
          assign(gamma, cons_list(index_uni(sym1__), nil_index_list()),
            gamma_flat__[(pos__ - 1)], "assigning variable gamma");
          current_statement__ = 1;
          pos__ = (pos__ + 1);}
      }
      for (int sym1__ = 1; sym1__ <= gamma_1dim__; ++sym1__) {
        vars__.emplace_back(gamma[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("gamma");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{
                                             static_cast<size_t>(gamma_1dim__)});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= gamma_1dim__; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "gamma" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= gamma_1dim__; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "gamma" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"gamma\",\"type\":{\"name\":\"vector\",\"length\":" << gamma_1dim__ << "},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"gamma\",\"type\":{\"name\":\"vector\",\"length\":" << gamma_1dim__ << "},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_mcml_poisson_region_namespace::model_mcml_poisson_region;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_mcml_poisson_region_namespace::profiles__;
}
#endif
#endif
