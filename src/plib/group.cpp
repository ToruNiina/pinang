/*!
  @file group.cpp
  @brief Define functions of class Group.

  Definitions of member or friend functions of class Group.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:41
  @copyright GNU Public License V3.0
*/


#include "group.hpp"

namespace pinang {

Group::Group()
{
  n_atom_ = 0;
  coordinates_.clear();
}

Group::Group(std::vector<Vec3d> v)
{
  coordinates_ = v;
  n_atom_ = coordinates_.size();
}

Group::Group(const Conformation& c, const Selection& s)
{
  if (s.get_size() > c.get_size()) {
    std::cout << " ~             PINANG :: group.hpp            ~ " << "\n";
    std::cerr << " ERROR: Sel number > conformation atom number." << "\n";
    exit(EXIT_SUCCESS);
  }

  Conformation c1 = c;
  int n = s.get_size();
  for (int i = 0; i < n; ++i) {
    int m = s.get_selection(i);
    coordinates_.push_back(c1.get_coordinate(m));
  }
  n_atom_ = n;
}

int Group::set_conformation(std::vector<Vec3d> v)
{
  int m = v.size();
  if (m != n_atom_ && n_atom_ > 0)
  {
    std::cout << " ~             PINANG :: group.hpp              ~ " << "\n";
    std::cerr << " ERROR: Wrong atom number when set conformation. " << "\n";
    return 1;
  } else {
    coordinates_ = v;
    n_atom_ = m;
    return 0;
  }
}

Vec3d Group::get_centroid() const
{
  Vec3d com;
  for (int i = 0; i < n_atom_; ++i) {
    com += coordinates_[i];
  }
  com /= n_atom_;
  return com;
}

}
