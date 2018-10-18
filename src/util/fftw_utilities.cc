/* Copyright (C) 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "fftw_utilities.hh"

#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

/*!
 * Return the Discrete Fourier Transform sample frequencies.
 */
std::vector<double> fftfreq(int M, double normalization) {
  std::vector<double> result(M);

  for (int i = 0; i <= M / 2; i++) {
    result[i] = i;
  }

  for (int i = M / 2 + 1; i < M; i++) {
    result[i] = M - i;
  }

  // normalize
  for (int i = 0; i < M; i++) {
    result[i] *= normalization;
  }

  return result;
}

//! \brief Fill `input` with zeros.
void clear_fftw_array(fftw_complex *input, int Nx, int Ny) {
  VecAccessor2D<fftw_complex> fftw_in(input, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      fftw_in(i, j)[0] = 0;
      fftw_in(i, j)[1] = 0;
    }
  }
}

//! @brief Copy `source` to `destination`.
void copy_fftw_array(fftw_complex *source, fftw_complex *destination, int Nx, int Ny) {
  VecAccessor2D<fftw_complex> S(source, Nx, Ny), D(destination, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      D(i, j)[0] = S(i, j)[0];
      D(i, j)[1] = S(i, j)[1];
    }
  }
}

//! Set the real part of output to input. Input has the size of My*Mx, embedded in the
//! bigger (output) grid of size Ny*Nx. Offsets i0 and j0 specify the location of the
//! subset to set.
/*!
 * Sets the imaginary part to zero.
 */
void set_real_part(Vec input,
                   double normalization,
                   int Mx, int My,
                   int Nx, int Ny,
                   int i0, int j0,
                   fftw_complex *output) {
  petsc::VecArray2D in(input, Mx, My);
  VecAccessor2D<fftw_complex> out(output, Nx, Ny, i0, j0);

  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j)[0] = in(i, j) * normalization;
      out(i, j)[1] = 0.0;
    }
  }
}


//! \brief Get the real part of input and put it in output.
/*!
 * See set_real_part for details.
 */
void get_real_part(fftw_complex *input,
                   double normalization,
                   int Mx, int My,
                   int Nx, int Ny,
                   int i0, int j0,
                   Vec output) {
  petsc::VecArray2D out(output, Mx, My);
  VecAccessor2D<fftw_complex> in(input, Nx, Ny, i0, j0);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j) = in(i, j)[0] * normalization;
    }
  }
}

} // end of namespace pism
