// -*-c++-*-

#ifndef PINANG_CONFORMATION_H
#define PINANG_CONFORMATION_H

#include "constants.h"
#include "myvec.h"

#include <vector>
#include <cmath>
#include <Eigen/Dense>

namespace pinang {

class Conformation
{
 public:
  Conformation();
  Conformation(std::vector<Eigen::Vector3d> v);
  virtual ~Conformation() {coordinates_.clear();}

  inline void reset();

  int get_size() {return n_atom_;}

  int set_conformation(std::vector<Eigen::Vector3d> v);
  inline const Eigen::Vector3d& get_coor(int n) const;

 protected:
  std::vector<Eigen::Vector3d> coordinates_;
  int n_atom_;
};

Conformation::Conformation()
{
  n_atom_ = 0;
  coordinates_.clear();
}

Conformation::Conformation(std::vector<Eigen::Vector3d> v)
{
  coordinates_ = v;
  n_atom_ = coordinates_.size();
}

inline void Conformation::reset()
{
  n_atom_ = 0;
  coordinates_.clear();
}

int Conformation::set_conformation(std::vector<Eigen::Vector3d> v)
{
  int m = v.size();
  if (m != n_atom_ && n_atom_ > 0)
  {
    std::cerr << " ERROR: Wrong atom number when set conformation. " << std::endl;
    return 1;
  } else {
    coordinates_ = v;
    n_atom_ = m;
    return 0;
  }
}

inline const Eigen::Vector3d& Conformation::get_coor(int n) const
{
  if ( n >= n_atom_ || n < 0)
  {
    std::cerr << " ERROR: Atom index out of range in Conformation. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return coordinates_[n];
  }
}

}
#endif
