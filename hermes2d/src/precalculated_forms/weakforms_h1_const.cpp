// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "weakforms_h1_const.h"
#include "api2d.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace ConstantWeakFormsH1
    {
      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::ConstantMatrixFormVol
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::ConstantMatrixFormVol
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormVol<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormVolTriangle.h1h1", 0, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormVolQuad.h1h1", 0, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormVolTriangle.l2l2", 0, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormVolQuad.l2l2", 0, 0);

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::~ConstantMatrixFormVol()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormVol<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantMatrixFormVol<Scalar>(*this);
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::ConstantMatrixFormDx
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::ConstantMatrixFormDx
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDx<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormDxTriangle.h1h1", 2, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormDxQuad.h1h1", 2, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormDxTriangle.l2l2", 2, 0);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormDxQuad.l2l2", 2, 0);

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::~ConstantMatrixFormDx()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDx<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDx<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDx<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantMatrixFormDx<Scalar>(*this);
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::ConstantMatrixFormDy
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::ConstantMatrixFormDy
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDy<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormDyTriangle.h1h1", 0, 2);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormDyQuad.h1h1", 0, 2);
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultMatrixFormDyTriangle.l2l2", 0, 2);
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "DefaultMatrixFormDyQuad.l2l2", 0, 2);

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::~ConstantMatrixFormDy()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDy<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDy<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDy<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantMatrixFormDy<Scalar>(*this);
      }


      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::ConstantVectorFormVol
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::ConstantVectorFormVol
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }
        
      template<typename Scalar>
      void ConstantVectorFormVol<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormVolTriangle.h1", 0, 0);
        this->set_h1_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormVolQuad.h1", 0, 0);
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormVolTriangle.l2", 0, 0);
        this->set_l2_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormVolQuad.l2", 0, 0);
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::~ConstantVectorFormVol()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormVol<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantVectorFormVol<Scalar>(*this);
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::ConstantVectorFormDx
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::ConstantVectorFormDx
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantVectorFormDx<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormDxTriangle.h1", 1, 0);
        this->set_h1_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormDxQuad.h1", 1, 0);
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormDxTriangle.l2", 1, 0);
        this->set_l2_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormDxQuad.l2", 1, 0);
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::~ConstantVectorFormDx()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormDx<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormDx<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormDx<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantVectorFormDx<Scalar>(*this);
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::ConstantVectorFormDy
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::ConstantVectorFormDy
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantVectorFormDy<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormDyTriangle.h1", 1, 0);
        this->set_h1_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormDyQuad.h1", 1, 0);
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "DefaultVectorFormDyTriangle.l2", 1, 0);
        this->set_l2_const_tables(HERMES_MODE_QUAD, "DefaultVectorFormDyQuad.l2", 1, 0);
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::~ConstantVectorFormDy()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormDy<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormDy<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormDy<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        return new ConstantVectorFormDy<Scalar>(*this);
      }

      template class HERMES_API ConstantMatrixFormVol<double>;
      template class HERMES_API ConstantMatrixFormVol<std::complex<double> >;
      template class HERMES_API ConstantMatrixFormDx<double>;
      template class HERMES_API ConstantMatrixFormDx<std::complex<double> >;
      template class HERMES_API ConstantMatrixFormDy<double>;
      template class HERMES_API ConstantMatrixFormDy<std::complex<double> >;
      template class HERMES_API ConstantVectorFormVol<double>;
      template class HERMES_API ConstantVectorFormVol<std::complex<double> >;
      template class HERMES_API ConstantVectorFormDx<double>;
      template class HERMES_API ConstantVectorFormDx<std::complex<double> >;
      template class HERMES_API ConstantVectorFormDy<double>;
      template class HERMES_API ConstantVectorFormDy<std::complex<double> >;
    };
  }
}
