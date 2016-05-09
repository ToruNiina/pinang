// -*-c++-*-

#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Core>

#include "atom.h"

namespace pinang {

class Residue
{
 public:
  Residue();
  virtual ~Residue() {atoms_.clear();}

  inline void reset();

  inline std::string get_resid_name() const;
  inline char get_short_name() const;
  void set_resid_name(const std::string& s);

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

  inline Atom& get_atom(int n);
  inline int add_atom(const Atom& a);
  inline int delete_atom(const int i);

  inline int get_residue_size() const;

  inline Atom& get_C_alpha();
  inline Atom& get_C_beta();
  inline Atom& get_P();
  inline Atom& get_S();
  inline Atom& get_B();

  void set_cg_na();

 protected:
  std::string resid_name_;
  char short_name_;
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

/*                _     _
//  _ __ ___  ___(_) __| |  _ __   __ _ _ __ ___   ___
// | '__/ _ \/ __| |/ _` | | '_ \ / _` | '_ ` _ \ / _ \
// | | |  __/\__ \ | (_| | | | | | (_| | | | | | |  __/
// |_|  \___||___/_|\__,_| |_| |_|\__,_|_| |_| |_|\___|
*/
inline std::string Residue::resid_name() const
{
  return resid_name_;
}
inline char Residue::short_name() const
{
  return short_name_;
}
void Residue::set_resid_name(const std::string& s)
{
  resid_name_ = s;
  char c = s[0];

  switch (c) {
    case 'A':
      if (resid_name_ == "ARG") {short_name_ = 'R'; charge_ = 1.0; mass_ = 156.19; chain_type_ = protein;}
      else if (resid_name_ == "ASP") {short_name_ = 'D'; charge_ = -1.0;
        mass_ = 115.09; chain_type_ = protein;}
      else if (resid_name_ == "ASN") {short_name_ = 'N'; mass_ = 114.11; chain_type_ = protein;}
      else if (resid_name_ == "ALA") {short_name_ = 'A'; mass_ = 71.09; chain_type_ = protein;}
      else if (resid_name_ == "A") {short_name_ = 'A'; chain_type_ = na;mass_ = 134.12;}
      break;
    case 'C':
      if (resid_name_ == "CYS") {short_name_ = 'C';  mass_ = 103.15;  chain_type_ = protein;}
      else if (resid_name_ == "C") {short_name_ = 'C'; chain_type_ = na;  mass_ = 110.09;}
      else if (resid_name_ == "CA") {short_name_ = 'c'; charge_ = 2.0; mass_ = 40.08; chain_type_ = ion;}
      break;
    case 'D':
      if      (resid_name_ == "DA" || resid_name_ == "DA3" || resid_name_ == "DA5") {
        short_name_ = 'A'; chain_type_ = DNA; mass_ = 134.12;}
      else if (resid_name_ == "DC" || resid_name_ == "DC3" || resid_name_ == "DC5") {
        short_name_ = 'C'; chain_type_ = DNA;  mass_ = 110.09;}
      else if (resid_name_ == "DG" || resid_name_ == "DG3" || resid_name_ == "DG5") {
        short_name_ = 'G'; chain_type_ = DNA;  mass_ = 150.12;}
      else if (resid_name_ == "DT" || resid_name_ == "DT3" || resid_name_ == "DT5") {
        short_name_ = 'T'; chain_type_ = DNA;  mass_ = 125.091;}
      break;
    case 'G':
      if (resid_name_ == "GLU") {short_name_ = 'E'; charge_ = -1.0; mass_ = 129.12; chain_type_ = protein;}
      else if (resid_name_ == "GLN") {short_name_ = 'Q';  mass_ = 128.14;  chain_type_ = protein;}
      else if (resid_name_ == "GLY") {short_name_ = 'G';  mass_ = 57.05; chain_type_ = protein;}
      else if (resid_name_ == "G") {short_name_ = 'G'; chain_type_ = na;  mass_ = 150.12;}
      break;
    case 'H':
      if (resid_name_ == "HIS") {short_name_ = 'H';  mass_ = 137.14; charge_ = 1.0; chain_type_ = protein;}
      else if (resid_name_ == "HOH") {short_name_ = 'w';   chain_type_ = water;}
      break;
    case 'I':
      if (resid_name_ == "ILE") {short_name_ = 'I';  mass_ = 113.16; chain_type_ = protein;}
      break;
    case 'L':
      if (resid_name_ == "LYS") {short_name_ = 'K'; charge_ = 1.0; mass_ = 128.17; chain_type_ = protein;}
      else if (resid_name_ == "LEU") {short_name_ = 'L';  mass_ = 113.16; chain_type_ = protein;}
      break;
    case 'M':
      if (resid_name_ == "MET") {short_name_ = 'M';  mass_ = 131.19; chain_type_ = protein;}
      else if (resid_name_ == "MG") {short_name_ = 'm';  mass_ = 24.305; charge_ = 2.0; chain_type_ = ion;}
      break;
    case 'P':
      if (resid_name_ == "PRO") {short_name_ = 'P';  mass_ = 97.12; chain_type_ = protein;}
      else if (resid_name_ == "PHE") {short_name_ = 'F';  mass_ = 147.18; chain_type_ = protein;}
      break;
    case 'R':
      if (resid_name_ == "RA") {short_name_ = 'A'; chain_type_ = RNA; mass_ = 134.12;}
      else if (resid_name_ == "RU") {short_name_ = 'U'; chain_type_ = RNA; mass_ = 111.08;}
      else if (resid_name_ == "RC") {short_name_ = 'C'; chain_type_ = RNA;  mass_ = 110.09;}
      else if (resid_name_ == "RG") {short_name_ = 'G';  chain_type_ = RNA; mass_ = 150.12;}
      break;
    case 'S':
      if (resid_name_ == "SER") {short_name_ = 'S';  mass_ = 87.08; chain_type_ = protein;}
      else if (resid_name_ == "SEC") {short_name_ = 'U'; chain_type_ = protein;}
      break;
    case 'T':
      if (resid_name_ == "TYR") {short_name_ = 'Y';  mass_ = 163.18; chain_type_ = protein;}
      else if (resid_name_ == "TRP") {short_name_ = 'W'; mass_ = 186.21; chain_type_ = protein;}
      else if (resid_name_ == "THR") {short_name_ = 'T'; mass_ = 101.11; chain_type_ = protein;}
      else if (resid_name_ == "T") {short_name_ = 'T';  chain_type_ = DNA; mass_ = 125.091;}
      break;
    case 'U':
      if (resid_name_ == "U") {short_name_ = 'U'; chain_type_ = RNA; mass_ = 111.08;}
      break;
    case 'V':
      if (resid_name_ == "VAL") {short_name_ = 'V'; mass_ = 99.14; chain_type_ = protein;}
      break;
    default:
      if (resid_name_ == "ZN") {short_name_ = 'z'; charge_ = 2.0; mass_ = 65.409; chain_type_ = ion;}
  }
}

/*      _           _         ___ ____
//  ___| |__   __ _(_)_ __   |_ _|  _ \
// / __| '_ \ / _` | | '_ \   | || | | |
//| (__| | | | (_| | | | | |  | || |_| |
// \___|_| |_|\__,_|_|_| |_| |___|____/
*/
inline char Residue::chain_ID() const
{
  return chain_ID_;
}
inline void Residue::set_chain_ID(char a)
{
  chain_ID_ = a;
}

/*       _           _         _
//   ___| |__   __ _(_)_ __   | |_ _   _ _ __   ___
//  / __| '_ \ / _` | | '_ \  | __| | | | '_ \ / _ \
// | (__| | | | (_| | | | | | | |_| |_| | |_) |  __/
//  \___|_| |_|\__,_|_|_| |_|  \__|\__, | .__/ \___|
//                                 |___/|_|
*/
inline ChainType Residue::chain_type() const
{
  return ChainTypeype_;
}
inline void Residue::set_ChainTypeype(ChainType a)
{
  ChainTypeype_ = a;
  if (a == DNA) {
    // std::string s = resid_name_;
    if      (resid_name_ == "A" || resid_name_ == "A3" || resid_name_ == "A5") {
      resid_name_ = "DA"; short_name_ = 'A';  mass_ = 134.12;}
    else if (resid_name_ == "C" || resid_name_ == "C3" || resid_name_ == "C5") {
      resid_name_ = "DC"; short_name_ = 'C';   mass_ = 110.09;}
    else if (resid_name_ == "G" || resid_name_ == "G3" || resid_name_ == "G5") {
      resid_name_ = "DG"; short_name_ = 'G';   mass_ = 150.12;}
    else if (resid_name_ == "T" || resid_name_ == "T3" || resid_name_ == "T5") {
      resid_name_ = "DT"; short_name_ = 'T';   mass_ = 125.091;}
    // std::cout << resid_name_ << std::endl;
    for (int i = 0; i < n_atom_; i++) {
      atoms_[i].set_resid_name(resid_name_);
    }
  }
}

/*                _     _   _           _
//  _ __ ___  ___(_) __| | (_)_ __   __| | _____  __
// | '__/ _ \/ __| |/ _` | | | '_ \ / _` |/ _ \ \/ /
// | | |  __/\__ \ | (_| | | | | | | (_| |  __/>  <
// |_|  \___||___/_|\__,_| |_|_| |_|\__,_|\___/_/\_\
*/
inline int Residue::resid_index() const
{
  return resid_index_;
}
inline void Residue::set_resid_index(int i)
{
  resid_index_ = i;
}

/*                _     _        _
//  _ __ ___  ___(_) __| |   ___| |__   __ _ _ __ __ _  ___
// | '__/ _ \/ __| |/ _` |  / __| '_ \ / _` | '__/ _` |/ _ \
// | | |  __/\__ \ | (_| | | (__| | | | (_| | | | (_| |  __/
// |_|  \___||___/_|\__,_|  \___|_| |_|\__,_|_|  \__, |\___|
//                                               |___/
*/
inline double Residue::resid_charge() const
{
  return charge_;
}
inline void Residue::set_resid_charge(double c)
{
  charge_ = c;
}

inline double Residue::resid_mass() const
{
  return mass_;
}
inline void Residue::set_resid_mass(double m)
{
  mass_ = m;
}

/*                     _
//  _ __ ___      __ _| |_ ___  _ __ ___
// | '_ ` _ \    / _` | __/ _ \| '_ ` _ \
// | | | | | |  | (_| | || (_) | | | | | |
// |_| |_| |_|___\__,_|\__\___/|_| |_| |_|
//          |_____|
*/
inline Atom& Residue::m_atom(int n)
{
  if (atoms_.empty())
  {
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~             PINANG :: Residue              ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cerr << "ERROR: No Atoms found in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    if (n >= int(atoms_.size()))
    {
      std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
      std::cout << " ~             PINANG :: Residue              ~ " << std::endl;
      std::cout << " ============================================== " << std::endl;
      std::cerr << "ERROR: Atom index out of range in Residue: "
                << resid_index_ << std::endl;
      exit(EXIT_SUCCESS);
    } else {
      return atoms_[n];
    }
  }
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


/*   ____         _       _
//  / ___|   __ _| |_ __ | |__   __ _
// | |      / _` | | '_ \| '_ \ / _` |
// | |___  | (_| | | |_) | | | | (_| |
//  \____|  \__,_|_| .__/|_| |_|\__,_|
//                 |_|
*/
inline Atom& Residue::m_C_alpha()
{
  return C_alpha_;
}

inline Atom& Residue::m_P()
{
  return P_;
}

inline Atom& Residue::m_S()
{
  return S_;
}

inline Atom& Residue::m_B()
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

/*                _     _       _
//  _ __ ___  ___(_) __| |  ___(_)_______
// | '__/ _ \/ __| |/ _` | / __| |_  / _ \
// | | |  __/\__ \ | (_| | \__ \ |/ /  __/
// |_|  \___||___/_|\__,_| |___/_/___\___|
*/
inline int Residue::m_residue_size() const
{
  return n_atom_;
}

// Residue------------------------------------------------------------------
inline Residue::Residue()
{
  resid_name_ = "";
  short_name_ = '0';
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
  short_name_ = '0';
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

/*              _               __
//   ___  _   _| |_ ___ _ __   / _|_   _ _ __   ___
//  / _ \| | | | __/ _ \ '__| | |_| | | | '_ \ / __|
// | (_) | |_| | ||  __/ |    |  _| |_| | | | | (__
//  \___/ \__,_|\__\___|_|    |_|  \__,_|_| |_|\___|
*/
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
