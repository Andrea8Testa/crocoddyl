///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_MULTIBODY_COST_BASE_HPP_
#define CROCODDYL_MULTIBODY_COST_BASE_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "crocoddyl/multibody/fwd.hpp"
#include "crocoddyl/multibody/states/multibody.hpp"
#include "crocoddyl/core/data-collector-base.hpp"
#include "crocoddyl/core/activation-base.hpp"
#include "crocoddyl/core/utils/to-string.hpp"
#include "crocoddyl/core/activations/quadratic.hpp"

namespace crocoddyl {

/**
 * @brief Abstract class for cost models
 *
 * In Crocoddyl, a cost model is defined by the scalar activation function \f$a(\cdot)\f$ and by the residual function
 * \f$\mathbf{r}(\cdot)\f$ as follows: \f[ cost = a(\mathbf{r}(\mathbf{x}, \mathbf{u})), \f] where
 * the residual function depends on the state point \f$\mathbf{x}\in\mathcal{X}\f$, which lies in the state manifold
 * described with a `nq`-tuple, its velocity \f$\dot{\mathbf{x}}\in T_{\mathbf{x}}\mathcal{X}\f$ that belongs to
 * the tangent space with `nv` dimension, and the control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$. The residual vector
 * is defined by \f$\mathbf{r}\in\mathbb{R}^{nr}\f$ where `nr` describes its dimension in the Euclidean space. On the
 * other hand, the activation function builds a cost value based on the definition of the residual vector. The residual
 * vector has to be specialized in a derived classes.
 *
 * The main computations are carring out in `calc` and `calcDiff` routines. `calc` computes the cost (and its residual)
 * and `calcDiff` computes the derivatives of the cost function (and its residual). Concretely speaking, `calcDiff`
 * builds a linear-quadratic approximation of the cost function with the form: \f$\mathbf{l_x}\in\mathbb{R}^{ndx}\f$,
 * \f$\mathbf{l_u}\in\mathbb{R}^{nu}\f$, \f$\mathbf{l_{xx}}\in\mathbb{R}^{ndx\times ndx}\f$,
 * \f$\mathbf{l_{xu}}\in\mathbb{R}^{ndx\times nu}\f$, \f$\mathbf{l_{uu}}\in\mathbb{R}^{nu\times nu}\f$ are the
 * Jacobians and Hessians, respectively.
 *
 * \sa `StateAbstractTpl`, `ActivationModelAbstractTpl`, `calc()`, `calcDiff()`, `createData()`
 */
template <typename _Scalar>
class CostModelAbstractTpl {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef CostDataAbstractTpl<Scalar> CostDataAbstract;
  typedef StateMultibodyTpl<Scalar> StateMultibody;
  typedef ActivationModelAbstractTpl<Scalar> ActivationModelAbstract;
  typedef ActivationModelQuadTpl<Scalar> ActivationModelQuad;
  typedef DataCollectorAbstractTpl<Scalar> DataCollectorAbstract;
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;

  /**
   * @brief Initialize the cost model
   *
   * @param[in] state       State of the multibody system
   * @param[in] activation  Activation model
   * @param[in] nu          Dimension of control vector
   */
  CostModelAbstractTpl(boost::shared_ptr<StateMultibody> state, boost::shared_ptr<ActivationModelAbstract> activation,
                       const std::size_t& nu);

  /**
   * @copybrief CostModelAbstractTpl()
   *
   * The default value for `nu` is equals to StateAbstractTpl::get_nv().
   *
   * @param[in] state       State of the multibody system
   * @param[in] activation  Activation model
   */
  CostModelAbstractTpl(boost::shared_ptr<StateMultibody> state, boost::shared_ptr<ActivationModelAbstract> activation);

  /**
   * @copybrief CostModelAbstractTpl()
   *
   * The default activation model is quadratic, i.e. ActivationModelQuadTpl::ActivationModelQuadTpl().
   *
   * @param[in] state  State of the multibody system
   * @param[in] nr     Dimension of residual vector
   * @param[in] nu     Dimension of control vector
   */
  CostModelAbstractTpl(boost::shared_ptr<StateMultibody> state, const std::size_t& nr, const std::size_t& nu);

  /**
   * @copybrief CostModelAbstractTpl()
   *
   * The default activation model is quadratic, i.e. ActivationModelQuadTpl::ActivationModelQuadTpl().
   * Furthermore, the default value for `nu` is StateAbstractTpl::get_nv().
   *
   * @param[in] state  State of the multibody system
   * @param[in] nr     Dimension of residual vector
   * @param[in] nu     Dimension of control vector
   */
  CostModelAbstractTpl(boost::shared_ptr<StateMultibody> state, const std::size_t& nr);
  virtual ~CostModelAbstractTpl();

  /**
   * @brief Compute the cost value and its residual vector
   *
   * @param[in] data  Cost data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
  virtual void calc(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u) = 0;

  /**
   * @brief Compute the Jacobian and Hessian of cost and its residual vector
   *
   * @param[in] data  Cost data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
  virtual void calcDiff(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
                        const Eigen::Ref<const VectorXs>& u) = 0;

  /**
   * @brief Create the cost data
   *
   * The default data contains objects to store the values of the cost, residual vector and their derivatives (first
   * and second order derivatives). However, it is possible to specialized this function is we need to create
   * additional data, for instance, to avoid dynamic memory allocation.
   *
   * @return the cost data
   */
  virtual boost::shared_ptr<CostDataAbstract> createData(DataCollectorAbstract* const data);

  /**
   * @copybrief calc()
   *
   * @param[in] data  Cost data
   * @param[in] x     State point
   */
  void calc(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x);

  /**
   * @copybrief calcDiff()
   *
   * @param[in] data  Cost data
   * @param[in] x     State point
   */
  void calcDiff(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x);

  /**
   * @brief Return the state
   */
  const boost::shared_ptr<StateMultibody>& get_state() const;

  /**
   * @brief Return the activation model
   */
  const boost::shared_ptr<ActivationModelAbstract>& get_activation() const;

  /**
   * @brief Return the dimension of the control input
   */
  const std::size_t& get_nu() const;

  /**
   * @brief Modify the cost reference
   */
  template <class ReferenceType>
  void set_reference(ReferenceType ref);

  /**
   * @brief Return the cost reference
   */
  template <class ReferenceType>
  ReferenceType get_reference() const;

 protected:
  /**
   * @copybrief set_reference()
   */
  virtual void set_referenceImpl(const std::type_info&, const void*);

  /**
   * @copybrief get_reference()
   */
  virtual void get_referenceImpl(const std::type_info&, void*) const;

  boost::shared_ptr<StateMultibody> state_;                //!< State description
  boost::shared_ptr<ActivationModelAbstract> activation_;  //!< Activation model
  std::size_t nu_;                                         //!< Control dimension
  VectorXs unone_;                                         //!< No control vector
};

template <typename _Scalar>
struct CostDataAbstractTpl {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ActivationDataAbstractTpl<Scalar> ActivationDataAbstract;
  typedef DataCollectorAbstractTpl<Scalar> DataCollectorAbstract;
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;

  template <template <typename Scalar> class Model>
  CostDataAbstractTpl(Model<Scalar>* const model, DataCollectorAbstract* const data)
      : shared(data),
        activation(model->get_activation()->createData()),
        cost(Scalar(0.)),
        Lx(model->get_state()->get_ndx()),
        Lu(model->get_nu()),
        Lxx(model->get_state()->get_ndx(), model->get_state()->get_ndx()),
        Lxu(model->get_state()->get_ndx(), model->get_nu()),
        Luu(model->get_nu(), model->get_nu()),
        r(model->get_activation()->get_nr()),
        Rx(model->get_activation()->get_nr(), model->get_state()->get_ndx()),
        Ru(model->get_activation()->get_nr(), model->get_nu()) {
    Lx.setZero();
    Lu.setZero();
    Lxx.setZero();
    Lxu.setZero();
    Luu.setZero();
    r.setZero();
    Rx.setZero();
    Ru.setZero();
  }
  virtual ~CostDataAbstractTpl() {}

  DataCollectorAbstract* shared;
  boost::shared_ptr<ActivationDataAbstract> activation;
  Scalar cost;
  VectorXs Lx;
  VectorXs Lu;
  MatrixXs Lxx;
  MatrixXs Lxu;
  MatrixXs Luu;
  VectorXs r;
  MatrixXs Rx;
  MatrixXs Ru;
};

}  // namespace crocoddyl

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "crocoddyl/multibody/cost-base.hxx"

#endif  // CROCODDYL_MULTIBODY_COST_BASE_HPP_
