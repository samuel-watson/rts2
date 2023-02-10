#ifndef EIGENEXT_H
#define EIGENEXT_H

#include <cmath> 
#include <RcppEigen.h>

// code helpfully taken from https://stackoverflow.com/questions/50027494/eigen-indices-of-dense-matrix-meeting-condition
// index code from http://eigen.tuxfamily.org/dox/TopicCustomizing_NullaryExpr.html#title1

namespace glmmr {
namespace Eigen_ext {

template<typename Func>
struct lambda_as_visitor_wrapper : Func {
  lambda_as_visitor_wrapper(const Func& f) : Func(f) {}
  template<typename S,typename I>
  void init(const S& v, I i, I j) { return Func::operator()(v,i,j); }
};

template<typename Mat, typename Func>
inline void visit_lambda(const Mat& m, const Func& f)
{
  lambda_as_visitor_wrapper<Func> visitor(f);
  m.visit(visitor);
}

inline Eigen::ArrayXi find(Eigen::ArrayXi vec, int n){
  std::vector<int> indices;
  glmmr::Eigen_ext::visit_lambda(vec,
                                 [&indices,n](int v, int i, int j) {
                                   if(v == n)
                                     indices.push_back(i);
                                 });
  return Eigen::Map<Eigen::ArrayXi>(indices.data(),indices.size());
}

template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                        RowIndexType::SizeAtCompileTime,
                        ColIndexType::SizeAtCompileTime,
                        ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                        RowIndexType::MaxSizeAtCompileTime,
                        ColIndexType::MaxSizeAtCompileTime> MatrixType;
  
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  
  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};

template <class ArgType, class RowIndexType, class ColIndexType>
Eigen::CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
mat_indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}

inline bool issympd(Eigen::MatrixXd mat){
  Eigen::LLT<Eigen::MatrixXd> lltOfA(mat);
  return lltOfA.info() == Eigen::NumericalIssue;
}

}
}



#endif