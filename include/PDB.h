// -*-c++-*-

#ifndef PINANG_PDB_H
#define PINANG_PDB_H

#include "model.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

namespace pinang{

class PDB
{
 public:
  PDB(const std::string& s);
  virtual ~PDB() {models_.clear();}

  inline std::string get_pdb_name() const;

  inline Model& get_model(unsigned int n);
  inline int get_n_models() const;

  inline void print_sequence(int n) const;
  inline void output_fasta(std::ostream & f_fasta) const;

  // void contact_map(double c);

 protected:
  std::string PDB_file_name_;
  std::vector<Model> models_;
  int n_model_;
};

PDB::PDB(const std::string& s)
{
  Atom atom_tmp;
  Residue resid_tmp;
  Chain chain_tmp;
  Model model_tmp;

  PDB_file_name_ = s;
  n_model_ = 0;
  models_.clear();

  std::ifstream ifile(PDB_file_name_.c_str());
  if (!ifile.is_open())
  {
    std::cerr << " ERROR: Cannot read file: " << s << std::endl;
    exit(EXIT_FAILURE);
  }

  while (ifile.good()) {
    ifile >> atom_tmp;
    if (ifile.fail())
    {
      break;
    }

    if (atom_tmp.get_atom_flag() == "MODEL ")
    {
      model_tmp.reset();
      model_tmp.set_model_ID(atom_tmp.get_serial());

      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "TER   ")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      model_tmp.add_chain(chain_tmp);

      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "ENDMDL")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      if (chain_tmp.get_chain_length() != 0)
      {
        model_tmp.add_chain(chain_tmp);
      }
      models_.push_back(model_tmp);
      n_model_++;

      model_tmp.reset();
      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "END   ")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      if (chain_tmp.get_chain_length() != 0)
      {
        model_tmp.add_chain(chain_tmp);
      }
      if (model_tmp.get_model_size() != 0)
      {
        models_.push_back(model_tmp);
        n_model_++;
      }

      model_tmp.reset();
      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "ATOM  " )
    {
      if (resid_tmp.add_atom(atom_tmp))
      {
        if (resid_tmp.get_residue_size() != 0)
        {
          chain_tmp.add_residue(resid_tmp);
          if (resid_tmp.get_atom(0).get_atom_flag() == "HETATM")
          {
            chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
            chain_tmp.set_chain_type(resid_tmp.get_chain_type());
            model_tmp.add_chain(chain_tmp);
            chain_tmp.reset();
          }

          resid_tmp.reset();
        }
        resid_tmp.set_residue_by_name(atom_tmp.get_resid_name());
        resid_tmp.set_chain_ID(atom_tmp.get_chain_ID());
        resid_tmp.set_resid_index(atom_tmp.get_resid_index());

        resid_tmp.add_atom(atom_tmp);
      }
    }
    if (atom_tmp.get_atom_flag() == "HETATM")
    {
      if (resid_tmp.add_atom(atom_tmp))
      {
        if (resid_tmp.get_residue_size() != 0)
        {
          chain_tmp.add_residue(resid_tmp);
          chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
          chain_tmp.set_chain_type(resid_tmp.get_chain_type());
          model_tmp.add_chain(chain_tmp);

          chain_tmp.reset();
          resid_tmp.reset();
        }
        resid_tmp.set_resid_index(atom_tmp.get_resid_index());
        resid_tmp.set_resid_name(atom_tmp.get_resid_name());
        resid_tmp.set_chain_ID(atom_tmp.get_chain_ID());
        resid_tmp.add_atom(atom_tmp);
      }
    }
  }

  ifile.close();
}

inline Model& PDB::get_model(unsigned int n)
{
  if (models_.empty()) {
    std::cerr << "ERROR: No Model found in this PDB: "
              << PDB_file_name_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= models_.size()) {
    std::cerr << "ERROR: Model number out of range in PDB: "
              << PDB_file_name_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return models_[n];
}

inline std::string PDB::get_pdb_name() const
{
  return PDB_file_name_;
}
inline int PDB::get_n_models() const
{
  return n_model_;
}


inline void PDB::print_sequence(int n) const
{
  if (n != 1 && n != 3)
  {
    std::cerr << " Usage: PINANG::PDB.print_sequence(): \n"
              << "       n = 1: 1-char residue name;\n"
              << "       n = 3: 3-char residue name.\n"
              << std::endl;
    exit(EXIT_SUCCESS);
  }
  models_[0].print_sequence(n);
}

inline void PDB::output_fasta(std::ostream & f_fasta) const
{
  std::string s = PDB_file_name_;
  for (int i = 0; i < 4; i++) {
    s.pop_back();
  }

  models_[0].output_fasta(f_fasta, s);
}


inline std::ostream& operator<<(std::ostream& o, PDB& p)
{
  int i = 0;
  for (i = 0; i < p.get_n_models(); i++) {
    o << p.get_model(i) << std::endl;
  }
  return o;
}

}

#endif
