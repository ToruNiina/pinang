// -*-c++-*-

#ifndef PINANG_CHAIN_H_
#define PINANG_CHAIN_H_

#include <iostream>
#include "residue.h"
#include "myvec.h"

namespace pinang {

class Chain
{
 public:
  Chain();
  virtual ~Chain() {residues_.clear();}

  inline void reset();

  inline char get_chain_ID() const;
  inline void set_chain_ID(char a);

  inline ChainType get_chain_type() const;
  inline void set_chain_type(ChainType ct);

  inline Residue& get_residue(int n);
  inline int add_residue(const Residue& r);

  inline int get_chain_length() const;

  inline void pr_seq(int n) const;
  inline void output_fasta(std::ostream & f_fasta, std::string s) const;
  inline void self_check();

  void output_cg_pos(std::ostream& o, int& n);
  void output_top_mass(std::ostream& o, int& n);
  void output_top_bond(std::ostream& o, int& n);
  void output_top_angle(std::ostream& o, int& n);
  void output_top_dihedral(std::ostream& o, int& n);
  void output_top_native(std::ostream& o);
  int get_native_contact_number();

  Chain operator+(Chain& other);

 protected:
  char chain_ID_;
  ChainType chain_type_;
  int n_residue_;
  std::vector<Residue> residues_;
};


inline char Chain::get_chain_ID() const
{
  return chain_ID_;
}
inline void Chain::set_chain_ID(char a)
{
  chain_ID_ = a;
}


inline ChainType Chain::get_chain_type() const
{
  return chain_type_;
}
inline void Chain::set_chain_type(ChainType a)
{
  chain_type_ = a;
  // for (int i = 0; i < n_residue_; i++)
  //   residues_[i].set_chain_type(a);
  // if (a == DNA) {
  //   for (int i = 0; i < residues_[0].get_residue_size(); i++) {
  //     std::string s0 = residues_[0].get_atom(i).atom_name();
  //     if (s0 == "OP3" || s0 == "P  " || s0 == "OP1" || s0 == "OP2"){
  //       if (!residues_[0].delete_atom(i))
  //         i--;
  //     }
  //   }
  // }
}


inline Residue& Chain::get_residue(int n)
{
  if (residues_.empty())
  {
    std::cerr << "ERROR: No Residues found in Chain: "
              << chain_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= int(residues_.size()))
  {
    std::cerr << "ERROR: Residue index out of range in Chain: "
              << chain_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return residues_[n];
}

inline int Chain::add_residue(const Residue& r)
{
  r.self_check();
  residues_.push_back(r);
  n_residue_++;
  return 0;
}

inline int Chain::get_chain_length() const
{
  return n_residue_;
}

inline void Chain::self_check()
{
  for (Residue& r : residues_) {
    if (r.get_chain_ID() != chain_ID_ || r.get_chain_type() != chain_type_) {
      std::cerr << "ERROR: Inconsistent chain ID or chain type in Chain "
                << chain_ID_ << std::endl;
      exit(EXIT_SUCCESS);
    }
  }
  if (chain_type_ == protein) {
    residues_[0].set_term_flag(-1);
    residues_[n_residue_ - 1].set_term_flag(1);
  }
  if (chain_type_ == DNA) {
    residues_[0].set_term_flag(5);
    residues_[n_residue_ - 1].set_term_flag(3);
  }
}

inline void Chain::pr_seq(int n) const
{
  if (chain_type_ == water)
  {
    return;
  }

  int j = 0;

  std::cout << " - Chain " << chain_ID_
            << " : " << n_residue_ << " residues."
            << std::endl;

  if (n == 1)
  {
    std::cout << " ";
    for (const Residue& r : residues_) {
      std::string s_tmp =  r.get_short_name();
      std::cout << std::setw(1) << s_tmp[1];
      j++;
      if (j%10 == 5)
        std::cout << " ";
      if (j%10 == 0)
      {
        std::cout << "  " << std::setw(4) << r.get_resid_index() << std::endl;
        std::cout << " ";
      }
    }
    std::cout << std::endl;
  } else if (n == 3) {
    for (const Residue& r : residues_) {
      std::cout << std::setw(4) << r.get_resid_name();
      j++;
      if (j%10 == 0)
      {
        std::cout << "  " << std::setw(4) << r.get_resid_index() << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

inline void Chain::output_fasta(std::ostream & f_fasta, std::string s0) const
{
  if (chain_type_ == water)
    return;
  if (n_residue_ <= 3)
    return;

  f_fasta << ">" << s0 << "_chain_"
          << chain_ID_ << "_type_"
          << chain_type_
          << std::endl;

  for (const Residue& r : residues_) {
    std::string s_tmp =  r.get_short_name();
    f_fasta << std::setw(1) << s_tmp[1];
  }
  f_fasta << std::endl;
}

void Chain::output_cg_pos(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  o << " - Chain " << chain_ID_ << " : "
    << n_residue_ << std::endl;
  int i = 0;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (Residue& r : residues_) {
      Eigen::Vector3d coor_CA;
      coor_CA = r.get_C_alpha().get_coordinates();
      o << std::setw(6) << ++n
        << std::setw(5) << r.get_resid_name()
        << std::setw(5) << r.get_resid_index() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << coor_CA(0) << " "
        << coor_CA(1) << " "
        << coor_CA(2) << " "
        << std::endl;
    }
  } else {
    for (Residue& r : residues_) {
      Eigen::Vector3d coor_P;
      Eigen::Vector3d coor_S;
      Eigen::Vector3d coor_B;
      if (r.get_term_flag() != 5) {
        coor_P = r.get_P().get_coordinates();
        o << std::setw(6) << ++n
          << std::setw(5) << "P"
          << std::setw(5) << r.get_resid_index() << "   "
          << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
          << std::setw(8) << coor_P.x() << " "
          << std::setw(8) << coor_P.y() << " "
          << std::setw(8) << coor_P.z() << " "
          << std::endl;
      }
      coor_S = r.get_S().get_coordinates();
      o << std::setw(6) << ++n
        << std::setw(5) << "S"
        << std::setw(5) << r.get_resid_index() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << std::setw(8) << coor_S.x() << " "
        << std::setw(8) << coor_S.y() << " "
        << std::setw(8) << coor_S.z() << " "
        << std::endl;
      coor_B = r.get_B().get_coordinates();
      o << std::setw(6) << ++n
        << std::setw(5) << r.get_short_name()
        << std::setw(5) << r.get_resid_index() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << std::setw(8) << coor_B.x() << " "
        << std::setw(8) << coor_B.y() << " "
        << std::setw(8) << coor_B.z() << " "
        << std::endl;
    }
  }
  o << std::endl;
}

void Chain::output_top_mass(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
  {
    return;
  }

  int i = 0;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (const Residue& r : residues_) {
      o << std::setw(11) << ++n
        << std::setw(8) << r.get_resid_index()
        << std::setw(8) << r.get_resid_name()
        << std::setw(8) << "CA"
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(10) << r.get_resid_mass()
        << std::setw(8)
        << r.get_resid_charge()
        << std::endl;
    }
  } else {
    for (const Residue& r : residues_) {
      if (r.get_term_flag() != 5) {
        o << std::setw(11) << ++n
          << std::setw(8) << r.get_resid_index()
          << std::setw(8) << r.get_resid_name()
          << std::setw(8) << "P"
          << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
          << std::setw(10) << 94.93
          << std::setw(8) << -0.6
          << std::endl;
      }
      o << std::setw(11) << ++n
        << std::setw(8) << r.get_resid_index()
        << std::setw(8) << r.get_resid_name()
        << std::setw(8) << "S"
        << std::setw(10) << 99.11
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(8) << 0.0
        << std::endl;
      o << std::setw(11) << ++n
        << std::setw(8) << r.get_resid_index()
        << std::setw(8) << r.get_resid_name()
        << std::setw(8) << "B"
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(10) << r.get_resid_mass()
        << std::setw(8) << 0.0
        << std::endl;
    }
  }
}

void Chain::output_top_bond(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
  {
    return;
  }

  int i = 0;
  double d = 0;
  double d_ps = 0;
  double d_sb = 0;
  double d_sp = 0;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_ - 1; i++) {
      d = resid_ca_distance(residues_[i], residues_[i+1]);
      n++;
      o << std::setw(8) << n
        << std::setw(6) << n + 1
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4)
        << std::setw(10) << d
        << std::setprecision(1)
        << std::setw(8) << k_K_bond
        << std::endl;
    }
    n++;
  } else {
    d_sb = atom_distance(residues_[0].get_S(), residues_[0].get_B());
    o << std::setw(8) << n+1
      << std::setw(6) << n+2
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
      << std::setw(10) << d_sb
      << std::setprecision(1) << std::setw(8) << k_K_bond
      << std::endl;
    d_sp = atom_distance(residues_[1].get_P(), residues_[0].get_S());
    o << std::setw(8) << n+1
      << std::setw(6) << n+3
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(4)
      << std::setw(10) << d_sp
      << std::setprecision(1)
      << std::setw(8) << k_K_bond
      << std::endl;
    for (i = 1; i < n_residue_ - 1; i++) {
      n += 3;
      d_ps = atom_distance(residues_[i].get_P(), residues_[i].get_S());
      o << std::setw(8) << n
        << std::setw(6) << n+1 << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4) << std::setw(10) << d_ps
        << std::setprecision(1) << std::setw(8) << k_K_bond
        << std::endl;
      d_sb = atom_distance(residues_[i].get_S(), residues_[i].get_B());
      o << std::setw(8) << n+1
        << std::setw(6) << n+2 << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4) << std::setw(10) << d_sb
        << std::setprecision(1) << std::setw(8) << k_K_bond
        << std::endl;
      d_sp = atom_distance(residues_[i].get_S(), residues_[i+1].get_P());
      o << std::setw(8) << n+1
        << std::setw(6) << n+3 << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4) << std::setw(10) << d_sp
        << std::setprecision(1) << std::setw(8) << k_K_bond
        << std::endl;
    }
    i = n_residue_ - 1;
    n += 3;
    d_ps = atom_distance(residues_[i].get_P(), residues_[i].get_S());
    o << std::setw(8) << n
      << std::setw(6) << n+1 << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(4) << std::setw(10) << d_ps
      << std::setprecision(1) << std::setw(8) << k_K_bond
      << std::endl;
    d_sb = atom_distance(residues_[i].get_S(), residues_[i].get_B());
    o << std::setw(8) << n+1
      << std::setw(6) << n+2 << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(4) << std::setw(10) << d_sb
      << std::setprecision(1) << std::setw(8) << k_K_bond
      << std::endl;
    n += 2;
  }
}

void Chain::output_top_angle(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double a = 0;
  Eigen::Vector3d v1, v2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-2; i++) {
      v1 = residues_[i].get_C_alpha().get_coordinates()
           - residues_[i+1].get_C_alpha().get_coordinates();
      v2 = residues_[i+2].get_C_alpha().get_coordinates()
           - residues_[i+1].get_C_alpha().get_coordinates();
      a = MyVec::vec_angle_deg (v1, v2);
      o << std::setw(8) << i+1+n
        << std::setw(6) << i+2+n
        << std::setw(6) << i+3+n
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4)
        << std::setw(12) << a
        << std::setprecision(1)
        << std::setw(8) << k_K_angle
        << std::endl;
    }
    n += n_residue_;
  } else {
    // ---------- angle BSP ----------
    v1 = residues_[0].get_B().get_coordinates()
         - residues_[0].get_S().get_coordinates();
    v2 = residues_[1].get_P().get_coordinates()
         - residues_[0].get_S().get_coordinates();
    a = MyVec::vec_angle_deg (v1, v2);
    o << std::setw(8) << n+2
      << std::setw(6) << n+1
      << std::setw(6) << n+3
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
      << std::setw(12) << a << std::setprecision(1)
      << std::setw(8) << k_K_angle << std::endl;
    // ---------- angle SPS ----------
    v1 = residues_[0].get_S().get_coordinates()
         - residues_[1].get_P().get_coordinates();
    v2 = residues_[1].get_S().get_coordinates()
         - residues_[1].get_P().get_coordinates();
    a = MyVec::vec_angle_deg (v1, v2);
    o << std::setw(8) << n+1
      << std::setw(6) << n+3
      << std::setw(6) << n+4
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
      << std::setw(12) << a << std::setprecision(1)
      << std::setw(8) << k_K_angle << std::endl;

    // -------------------- loop --------------------
    for (i = 1; i < n_residue_-1; i++) {
      n += 3;
      // ---------- angle PSB ----------
      v1 = residues_[i].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      v2 = residues_[i].get_B().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      a = MyVec::vec_angle_deg (v1, v2);
      o << std::setw(8) << n
        << std::setw(6) << n+1
        << std::setw(6) << n+2
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
        << std::setw(12) << a << std::setprecision(1)
        << std::setw(8) << k_K_angle << std::endl;
      // ---------- angle PSP ----------
      v1 = residues_[i].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      v2 = residues_[i+1].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      a = MyVec::vec_angle_deg (v1, v2);
      o << std::setw(8) << n
        << std::setw(6) << n+1
        << std::setw(6) << n+3
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
        << std::setw(12) << a << std::setprecision(1)
        << std::setw(8) << k_K_angle << std::endl;
      // ---------- angle BSP ----------
      v1 = residues_[i].get_B().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      v2 = residues_[i+1].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      a = MyVec::vec_angle_deg (v1, v2);
      o << std::setw(8) << n+2
        << std::setw(6) << n+1
        << std::setw(6) << n+3
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
        << std::setw(12) << a << std::setprecision(1)
        << std::setw(8) << k_K_angle << std::endl;
      // ---------- angle SPS ----------
      v1 = residues_[i].get_S().get_coordinates()
           - residues_[i+1].get_P().get_coordinates();
      v2 = residues_[i+1].get_S().get_coordinates()
           - residues_[i+1].get_P().get_coordinates();
      a = MyVec::vec_angle_deg (v1, v2);
      o << std::setw(8) << n+1
        << std::setw(6) << n+3
        << std::setw(6) << n+4
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
        << std::setw(12) << a << std::setprecision(1)
        << std::setw(8) << k_K_angle << std::endl;
    }
    n += 3;
    i = n_residue_ - 1;
    // ---------- angle PSB ----------
    v1 = residues_[i].get_P().get_coordinates()
         - residues_[i].get_S().get_coordinates();
    v2 = residues_[i].get_B().get_coordinates()
         - residues_[i].get_S().get_coordinates();
    a = MyVec::vec_angle_deg (v1, v2);
    o << std::setw(8) << n
      << std::setw(6) << n+1
      << std::setw(6) << n+2
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
      << std::setw(12) << a << std::setprecision(1)
      << std::setw(8) << k_K_angle << std::endl;
    n += 2;
  }
}

void Chain::output_top_dihedral(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double d = 0;           // dihedral
  Eigen::Vector3d v1, v2, v3, n1, n2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-3; i++) {
      v1 = residues_[i].get_C_alpha().get_coordinates()
           - residues_[i+1].get_C_alpha().get_coordinates();
      v2 = residues_[i+2].get_C_alpha().get_coordinates()
           - residues_[i+1].get_C_alpha().get_coordinates();
      v3 = residues_[i+2].get_C_alpha().get_coordinates()
           - residues_[i+3].get_C_alpha().get_coordinates();
      n1 = v1.cross(v2);
      n2 = v2.cross(v3);
      d = MyVec::vec_angle_deg (n1, n2);
      o << std::setw(8) << i+1+n
        << std::setw(6) << i+2+n
        << std::setw(6) << i+3+n
        << std::setw(6) << i+4+n
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4)
        << std::setw(12) << d
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1
        << std::setw(8) << k_K_dihedral_3
        << std::endl;
    }
    n += n_residue_;
  } else {
    v1 = residues_[0].get_S().get_coordinates()
         - residues_[1].get_P().get_coordinates();
    v2 = residues_[1].get_S().get_coordinates()
         - residues_[1].get_P().get_coordinates();
    v3 = residues_[1].get_S().get_coordinates()
         - residues_[2].get_P().get_coordinates();
    n1 = v1.cross(v2);
    n2 = v2.cross(v3);
    d = MyVec::vec_angle_deg (n1, n2);
    o << std::setw(8) << n + 1
      << std::setw(6) << n + 3
      << std::setw(6) << n + 4
      << std::setw(6) << n + 6
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(4)
      << std::setw(12) << d
      << std::setprecision(1)
      << std::setw(8) << k_K_dihedral_1
      << std::setw(8) << k_K_dihedral_3
      << std::endl;

    for (i = 1; i < n_residue_-1; i++) {
      n += 3;
      // ---------- PSPS ----------
      v1 = residues_[i].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      v2 = residues_[i+1].get_P().get_coordinates()
           - residues_[i].get_S().get_coordinates();
      v3 = residues_[i+1].get_P().get_coordinates()
           - residues_[i+1].get_S().get_coordinates();
      n1 = v1.cross(v2);
      n2 = v2.cross(v3);
      d = MyVec::vec_angle_deg (n1, n2);
      o << std::setw(8) << n
        << std::setw(6) << n + 1
        << std::setw(6) << n + 3
        << std::setw(6) << 4 + n
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4)
        << std::setw(12) << d
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1
        << std::setw(8) << k_K_dihedral_3
        << std::endl;

      if (i == n_residue_ - 2)
        break;
      // ---------- SPSP ----------
      v1 = residues_[i].get_S().get_coordinates()
           - residues_[i+1].get_P().get_coordinates();
      v2 = residues_[i+1].get_S().get_coordinates()
           - residues_[i+1].get_P().get_coordinates();
      v3 = residues_[i+1].get_S().get_coordinates()
           - residues_[i+2].get_P().get_coordinates();
      n1 = v1.cross(v2);
      n2 = v2.cross(v3);
      d = MyVec::vec_angle_deg (n1, n2);
      o << std::setw(8) << n + 1
        << std::setw(6) << n + 3
        << std::setw(6) << n + 4
        << std::setw(6) << n + 6
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(4)
        << std::setw(12) << d
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1
        << std::setw(8) << k_K_dihedral_3
        << std::endl;
    }
    n += 5;
  }
}

void Chain::output_top_native(std::ostream& o)
{
  int i = 0, j = 0;
  double d = -1, f = -1;
  ChainType cti, ctj;
  for (i = 0; i < n_residue_-4; i++) {
    cti = residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; j++) {
      ctj = residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = resid_min_distance(residues_[i], residues_[j]);
      if ( d < g_cutoff)
      {
        f = resid_ca_distance(residues_[i], residues_[j]);
        o << std::setw(8) << i+1
          << std::setw(6) << j+1
          << std::setiosflags(std::ios_base::fixed)
          << std::setprecision(2)
          << std::setw(8) << k_K_native
          << std::setprecision(4)
          << std::setw(10) << f
          << std::endl;
      }
    }
  }
}

int Chain::get_native_contact_number()
{
  int i = 0, j = 0;
  double d = -1;
  int n = 0;
  ChainType cti, ctj;

  for (i = 0; i < n_residue_-4; i++) {
    cti = residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; j++) {
      ctj = residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = resid_min_distance(residues_[i], residues_[j]);
      if ( d < g_cutoff)
        n++;
    }
  }
  return n;
}

// Chain -------------------------------------------------------------------
inline Chain::Chain()
{
  chain_ID_ = -1;
  chain_type_ = none;
  residues_.clear();
  n_residue_ = 0;
}

inline void Chain::reset()
{
  chain_ID_ = -1;
  chain_type_ = none;
  residues_.clear();
  n_residue_ = 0;
}

Chain Chain::operator+(Chain& other)
{
  int i = 0;
  Chain c1;
  c1.set_chain_ID(other.chain_ID_);
  if (n_residue_ > 0)
    for (i = 0; i < n_residue_; i++) {
      c1.residues_.push_back(residues_[i]);
      c1.n_residue_++;
    }
  if (other.n_residue_ > 0)
    for (i = 0; i < other.n_residue_; i++) {
      c1.residues_.push_back(other.residues_[i]);
      c1.n_residue_++;
    }

  return c1;
}

inline std::ostream& operator<<(std::ostream& o, Chain& c)
{
  int i = 0;
  for (i = 0; i < c.get_chain_length(); i++) {
    o << c.get_residue(i) ;
  }
  o << "TER   " << std::endl;
  return o;
}

}

#endif
