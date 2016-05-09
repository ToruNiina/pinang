// -*-c++-*-

#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Core>

#include "constants.h"
#include "atom.h"

namespace pinang {

class Residue
{
 public:
  Residue();
  virtual ~Residue() {atoms_.clear();}

  inline void reset();

  inline std::string get_resid_name() const;
  inline std::string get_short_name() const;
  void set_resid_name(const std::string& s);
  void set_residue_by_name(const std::string& s);

  inline char get_chain_ID() const;
  inline void set_chain_ID(char a);

  inline ChainType get_chain_type() const;
  inline void set_ChainTypeype(ChainType ct);

  inline int get_resid_index() const;
  inline void set_resid_index(int i);

  inline double get_resid_charge() const;
  inline void set_resid_charge(double c);

  inline double get_resid_mass() const;
  inline void set_resid_mass(double c);

  inline const Atom& get_atom(int n) const;
  inline int set_atom(int n, const Atom& a);
  inline int add_atom(const Atom& a);
  inline int delete_atom(const int i);

  inline int get_residue_size() const;

  inline const Atom& get_C_alpha() const;
  inline const Atom& get_C_beta() const;
  inline const Atom& get_P() const;
  inline const Atom& get_S() const;
  inline const Atom& get_B() const;

  void set_cg_na();

 protected:
  std::string resid_name_;
  std::string short_name_;
  char chain_ID_;
  int resid_index_;
  std::vector<Atom> atoms_;
  int n_atom_;
  double charge_;
  double mass_;

  Atom C_alpha_;
  Atom C_beta_;
  Atom P_;
  Atom S_;
  Atom B_;
  ChainType chain_type_;
};


inline std::string Residue::get_resid_name() const
{
  return resid_name_;
}
inline std::string Residue::get_short_name() const
{
  return short_name_;
}
void Residue::set_resid_name(const std::string& s)
{
  resid_name_ = s;
}
void Residue::set_residue_by_name(const std::string& s)
{
  PhysicalProperty p;
  resid_name_ = s;
  short_name_ = p.get_short_name(s);
  chain_type_ = p.get_chain_type(s);
  mass_ = p.get_mass(s);
  charge_ = p.get_charge(s);
}


inline char Residue::get_chain_ID() const
{
  return chain_ID_;
}
inline void Residue::set_chain_ID(char a)
{
  chain_ID_ = a;
}


inline ChainType Residue::get_chain_type() const
{
  return ChainTypeype_;
}
inline void Residue::set_ChainTypeype(ChainType a)
{
  ChainTypeype_ = a;
  // if (a == DNA) {
  //   if      (resid_name_ == "A" || resid_name_ == "A3" || resid_name_ == "A5") {
  //     resid_name_ = "DA"; short_name_ = 'A';  mass_ = 134.12;}
  //   else if (resid_name_ == "C" || resid_name_ == "C3" || resid_name_ == "C5") {
  //     resid_name_ = "DC"; short_name_ = 'C';   mass_ = 110.09;}
  //   else if (resid_name_ == "G" || resid_name_ == "G3" || resid_name_ == "G5") {
  //     resid_name_ = "DG"; short_name_ = 'G';   mass_ = 150.12;}
  //   else if (resid_name_ == "T" || resid_name_ == "T3" || resid_name_ == "T5") {
  //     resid_name_ = "DT"; short_name_ = 'T';   mass_ = 125.091;}
  //   for (int i = 0; i < n_atom_; i++) {
  //     atoms_[i].set_resid_name(resid_name_);
  //   }
  // }
}


inline int Residue::get_resid_index() const
{
  return resid_index_;
}
inline void Residue::set_resid_index(int i)
{
  resid_index_ = i;
}


inline double Residue::get_resid_charge() const
{
  return charge_;
}
inline void Residue::set_resid_charge(double c)
{
  charge_ = c;
}

inline double Residue::get_resid_mass() const
{
  return mass_;
}
inline void Residue::set_resid_mass(double m)
{
  mass_ = m;
}


inline const Atom& Residue::get_atom(int n) cosnt
{
  if (atoms_.empty())
  {
    std::cerr << "ERROR: No Atoms in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  } 
  if (n >= int(atoms_.size()))
  {
    std::cerr << "ERROR: Atom index out of range in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return atoms_[n];
}

inline int Residue::add_atom(const Atom& a)
{
  if (a.resid_index() != resid_index_)
  {
    return 1;
  }
  for (int i = 0; i < n_atom_; i++) {
    Atom b = atoms_[i];
    if (a.atom_name() == b.atom_name())
      return 0;
  }
  // -------------------- add atom --------------------
  atoms_.push_back(a);
  // -------------------- add atom --------------------
  if (a.atom_name() == "CA ")
  {
    C_alpha_ = a;
  }
  if (a.atom_name() == "C3'" || a.atom_name() == "S  " || a.atom_name() == "DS ")
  {
    S_ = a;
  }
  if (a.atom_name() == "P  " || a.atom_name() == "DP ")
  {
    P_ = a;
  }
  if (a.atom_name() == "N1 " || a.atom_name() == "B  " || a.atom_name() == "DB ")
  {
    B_ = a;
  }
  if (a.atom_flag() == "HETATM" && a.element() != "H")
  {
    C_alpha_ = a;
  }
  n_atom_++;
  return 0;
}
inline int Residue::delete_atom(const int i)
{
  if (i >= n_atom_){
    return 1;
  }
  atoms_.erase(atoms_.begin()+i);
  n_atom_--;
  return 0;
}



inline const Atom& Residue::get_C_alpha() const
{
  return C_alpha_;
}

inline const Atom& Residue::get_C_beta() const
{
  return C_beta_;
}

inline const Atom& Residue::get_P() const
{
  return P_;
}

inline const Atom& Residue::get_S() const
{
  return S_;
}

inline const Atom& Residue::get_B() const
{
  return B_;
}

void Residue::set_cg_na()
{
  int i = 0;
  Vec3d coor_P(0,0,0);
  Vec3d coor_C5p(0,0,0),
      coor_C4p(0,0,0),
      coor_O4p(0,0,0),
      coor_C1p(0,0,0),
      coor_C2p(0,0,0),
      coor_C3p(0,0,0),
      coor_O2p(0,0,0),
      coor_O3p(0,0,0),
      coor_O5p(0,0,0);
  Vec3d coor_N1(0,0,0),
      coor_C2(0,0,0),
      coor_N3(0,0,0),
      coor_C4(0,0,0),
      coor_C5(0,0,0),
      coor_C6(0,0,0),
      coor_O2(0,0,0),
      coor_N4(0,0,0),
      coor_O4(0,0,0),
      coor_N6(0,0,0),
      coor_O6(0,0,0),
      coor_N2(0,0,0),
      coor_N7(0,0,0),
      coor_C8(0,0,0),
      coor_N9(0,0,0);
  Vec3d com_P(0,0,0);
  Vec3d com_S(0,0,0);
  Vec3d com_B(0,0,0);
  double mass_C = 12.011;
  double mass_O = 15.999;
  double mass_N = 14.001;
  int n_cs=0;
  int n_cb=0;
  int n_os=0;
  int n_ob=0;
  int n_nb=0;

  if (ChainTypeype_ != DNA && ChainTypeype_ != RNA && ChainTypeype_ != na)
  {
    return;
  }
  for (i = 0; i < n_atom_; i++) {
    std::string aname = atoms_[i].atom_name();
    char c = aname[0];
    switch (c) {
      case 'C':
        if (aname == "C5'") {coor_C5p = atoms_[i].coordinates(); n_cs++;}
        else if (aname == "C1'") {coor_C1p = atoms_[i].coordinates(); n_cs++;}
        else if (aname == "C2'") {coor_C2p = atoms_[i].coordinates(); n_cs++;}
        else if (aname == "C3'") {coor_C3p = atoms_[i].coordinates(); n_cs++;}
        else if (aname == "C4'") {coor_C4p = atoms_[i].coordinates(); n_cs++;}
        else if (aname == "C2 ") { coor_C2 = atoms_[i].coordinates(); n_cb++;}
        else if (aname == "C4 ") { coor_C4 = atoms_[i].coordinates(); n_cb++;}
        else if (aname == "C5 ") { coor_C5 = atoms_[i].coordinates(); n_cb++;}
        else if (aname == "C6 ") { coor_C6 = atoms_[i].coordinates(); n_cb++;}
        else if (aname == "C8 ") { coor_C8 = atoms_[i].coordinates(); n_cb++;}
        break;
      case 'O':
        if (aname == "O4'") {coor_O4p = atoms_[i].coordinates(); n_os++;}
        else if (aname == "O2'") {coor_O2p = atoms_[i].coordinates(); n_os++;}
        else if (aname == "O3'") {coor_O3p = atoms_[i].coordinates();}
        else if (aname == "O5'") {coor_O5p = atoms_[i].coordinates();}
        else if (aname == "O2 ") {coor_O2 = atoms_[i].coordinates(); n_ob++;}
        else if (aname == "O4 ") {coor_O4 = atoms_[i].coordinates(); n_ob++;}
        else if (aname == "O6 ") {coor_O6 = atoms_[i].coordinates(); n_ob++;}
        break;
      case 'N':
        if (aname == "N1 ") {coor_N1 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N2 ") {coor_N2 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N3 ") {coor_N3 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N4 ") {coor_N4 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N6 ") {coor_N6 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N7 ") {coor_N7 = atoms_[i].coordinates(); n_nb++;}
        else if (aname == "N9 ") {coor_N9 = atoms_[i].coordinates(); n_nb++;}
        break;
      default:
        if (aname == "P  ") coor_P = atoms_[i].coordinates();
    }
  }
  com_P = coor_P;
  com_S = ( (coor_C1p + coor_C2p + coor_C3p + coor_C4p + coor_C5p)
            * mass_C
            + ( coor_O2p + coor_O4p ) * mass_O ) * (1/( n_os * mass_O + n_cs * mass_C ));
  com_B = ( (coor_N1 + coor_N2 + coor_N3 + coor_N4 + coor_N6 + coor_N7 + coor_N9)
            * mass_N
            + (coor_O2 + coor_O4 + coor_O6) * mass_O
            + (coor_C2 + coor_C4 + coor_C5 + coor_C6 + coor_C8)
            * mass_C )
          * (1/(n_nb*mass_N + n_ob*mass_O + n_cb*mass_C));
  P_.set_coords(com_P);
  if (S_.atom_name() != "S  ")
    S_.set_coords(com_S);
  if (B_.atom_name() != "B  ")
    B_.set_coords(com_B);
}


inline int Residue::get_residue_size() const
{
  return n_atom_;
}

// Residue------------------------------------------------------------------
inline Residue::Residue()
{
  resid_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  resid_index_ = -1;
  atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;

  C_alpha_.reset();
  P_.reset();
  S_.reset();
  B_.reset();
  ChainTypeype_ = none;
}

inline void Residue::reset()
{
  resid_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  resid_index_ = -1;
  atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;

  ChainTypeype_ = none;
  P_.reset();
  S_.reset();
  B_.reset();
  C_alpha_.reset();
}


// Other functions -----------------------------------------------------------
inline std::ostream& operator<<(std::ostream& o, Residue& r)
{
  // o << "Residue "
  //   << std::setw(4) << r.resid_index() << ":  "
  //   << std::setw(3) << r.resid_name() << std::endl;
  int i = 0;
  for (i = 0; i < r.m_residue_size(); i++) {
    o << r.m_atom(i) << std::endl;
  }
  return o;
}

inline double resid_min_distance (Residue& r1, Residue& r2)
{
  int i, j;
  double d = -1;           // distance;
  double f = 0;           // tmp distance;
  for (i = 0; i < r1.m_residue_size(); i++) {
    if (r1.m_atom(i).element() == "H")
      continue;
    for (j = 0; j < r2.m_residue_size(); j++) {
      if (r2.m_atom(j).element() == "H")
        continue;
      f = atom_distance(r1.m_atom(i), r2.m_atom(j));
      if (d < 0 || d > f)
      {
        d = f;
      }
    }
  }
  return d;
}

inline double resid_ca_distance (Residue& r1, Residue& r2)
{
  double d = -1;           // distance;
  d = atom_distance(r1.m_C_alpha(), r2.m_C_alpha());
  return d;
}

}

#endif
