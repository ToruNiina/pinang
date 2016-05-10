// -*-c++-*-

#ifndef PINANG_MYVEC_H_
#define PINANG_MYVEC_H_

#include "constants.h"

#include <Eigen/Dense>
#include <cmath>

namespace pinang{

class MyVec
{
 public:
  static double vec_distance (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
  static double vec_angle (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
  static double vec_angle_deg (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
};

double MyVec::vec_distance (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2)
{
  // distance between two Points
  Eigen::VectorXd v3 = v1-v2;
  return v3.norm();
}

double MyVec::vec_angle (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2)
{
  double d1 = (v1.dot(v2))/(v1.norm()*v2.norm());
  if (d1 >= 0.99999)
    return 0;
  if (d1 <= -0.99999)
    return 3.1415926;
  return acos(d1);
}

double MyVec::vec_angle_deg (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2)
{
  double d1 = (v1.dot(v2))/(v1.norm()*v2.norm());
  if (d1 >= 0.99999)
    return 0;
  if (d1 <= -0.99999)
    return 180;
  double ang = acos(d1);
  return (180 * ang / k_pi);
}

}
#endif
