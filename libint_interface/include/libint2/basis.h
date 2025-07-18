/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_basis_h_
#define _libint2_src_lib_libint_basis_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "libint2/basis.h requires C++11 support"
#endif
#include <libint2/cxxapi.h>

#include <cerrno>
#include <iostream>
#include <fstream>
#include <locale>
#include <vector>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <libint2.h>
#include <libint2/shell.h>
#include <libint2/atom.h>

namespace libint2 {

  /// Computes the number of basis functions in a range of shells
  /// @tparam ShellRange a range type
  /// @param[in] shells a sequence of shells
  /// @return the number of basis functions
  template <typename ShellRange> size_t nbf(const ShellRange& shells) {
    size_t n = 0;
    for (auto&& shell: shells)
      n += shell.size();
    return n;
  }

  /// Computes the maximum number of primitives in any Shell among a range of shells
  /// @tparam ShellRange a range type
  /// @param[in] shells a sequence of shells
  /// @return the maximum number of primitives
  template <typename ShellRange> size_t max_nprim(const ShellRange& shells) {
    size_t n = 0;
    for (auto&& shell: shells)
      n = std::max(shell.nprim(), n);
    return n;
  }

  /// Computes the maximum angular momentum quantum number @c l in any Shell among a range of shells
  /// @tparam ShellRange a range type
  /// @param[in] shells a sequence of shells
  /// @return the maximum angular momentum
  template <typename ShellRange> int max_l(const ShellRange& shells) {
    int l = 0;
    for (auto&& shell: shells)
      for (auto&& c: shell.contr)
        l = std::max(c.l, l);
    return l;
  }

  /// BasisSet is a slightly decorated \c std::vector of \c libint2::Shell objects.
  class BasisSet : private std::vector<libint2::Shell> {
    public:
      using base_type = std::vector<libint2::Shell>;

      BasisSet() : name_(""), nbf_(-1), max_nprim_(0), max_l_(-1) {}
      BasisSet(const BasisSet&) = default;
      BasisSet(BasisSet&& other) :
        std::vector<libint2::Shell>(std::move(other)),
        name_(std::move(other.name_)),
        nbf_(other.nbf_),
        max_nprim_(other.max_nprim_),
        max_l_(other.max_l_),
        shell2bf_(std::move(other.shell2bf_))
      {
      }
      BasisSet(const base_type& other) : base_type(other) { init(); }
      BasisSet(base_type&& other) : base_type(std::move(other)) { init(); }
      ~BasisSet() = default;
      BasisSet& operator=(const BasisSet&) = default;

      /// @brief Construct from the basis set name and a vector of atoms.

      /**
       * @param[in] name the basis set name
       * @param[in] atoms \c std::vector of Atom objects
       * @param[in] throw_if_no_match If true, and the basis is not found for this atomic number, throw a std::logic_error.
       *       Otherwise omit the basis quietly.
       * @throw std::logic_error if throw_if_no_match is true and for at least one atom no matching basis is found.
       * @throw std::ios_base::failure if throw_if_no_match is true and the basis file could not be read
       * \note All instances of the same chemical element receive the same basis set.
       * \note \c name will be "canonicalized" using BasisSet::canonicalize(name) to
       *       produce the file name where the basis will be sought. This file needs to contain
       *       the basis set definition in Gaussian94 format (see \c lib/basis directory for examples).
       *       The expected location of the file is determined by BasisSet::data_path as follows:
       *       <ol>
       *         <li> specified by LIBINT_DATA_PATH environmental variable, if defined </li>
       *         <li> specified by DATADIR macro variable, if defined </li>
       *         <li> specified by SRCDATADIR macro variable, if defined </li>
       *         <li> hardwired to directory \c NONE/share/libint/2.8.0/basis </li>
       *       </ol>
       */
      BasisSet(std::string name,
               const std::vector<Atom>& atoms,
               const bool throw_if_no_match = false) : name_(std::move(name)) {

        // read in the library file contents
        std::string basis_lib_path = data_path();
        auto canonical_name = canonicalize_name(name_);
        // some basis sets use cartesian d shells by convention, the convention is taken from Gaussian09
        auto force_cartesian_d = gaussian_cartesian_d_convention(canonical_name);

        // parse the name into components
        std::vector<std::string> basis_component_names = decompose_name_into_components(canonical_name);

        // ref_shells[component_idx][Z] => vector of Shells
        std::vector<std::vector<std::vector<libint2::Shell>>> component_basis_sets;
        component_basis_sets.reserve(basis_component_names.size());

        // read in ALL basis set components
        for(const auto& basis_component_name: basis_component_names) {
          auto file_dot_g94 = basis_lib_path + "/" + basis_component_name + ".g94";

          // use same cartesian_d convention for all components!
          component_basis_sets.emplace_back(read_g94_basis_library(file_dot_g94, force_cartesian_d, throw_if_no_match));
        }

        // for each atom find the corresponding basis components
        for(auto a=0ul; a<atoms.size(); ++a) {

          const std::size_t Z = atoms[a].atomic_number;

          // add each component in order
          for(auto comp_idx=0ul; comp_idx!=component_basis_sets.size(); ++comp_idx) {
            const auto& component_basis_set = component_basis_sets[comp_idx];
            if (!component_basis_set.at(Z).empty()) {  // found? add shells in order
              for(auto s: component_basis_set.at(Z)) {
                this->push_back(std::move(s));
                this->back().move({{atoms[a].x, atoms[a].y, atoms[a].z}});
              } // shell loop
            }
            else if (throw_if_no_match) {  // not found? throw, if needed
              std::string errmsg(std::string("did not find the basis for this Z in ") +
                      basis_lib_path + "/" + basis_component_names[comp_idx] + ".g94");
              throw std::logic_error(errmsg);
            }
          } // basis component loop
        } // atom loop

        init();
      }

      /// @brief Construct from a vector of atoms and the per-element basis set specification

      /**
       * @param[in] atoms \c std::vector of Atom objects
       * @param[in] element_bases vector of shell sequences for each element; if @c throw_if_no_match is false,
       *            ok for @c element_bases[Z] to be nonempty or
       *            not exist (e.g. if @c element_bases.size()<=Z )
       * @param[in] name the basis set name
       * @param[in] throw_if_no_match If true, and the basis is not found for this atomic number, throw a std::logic_error.
       *       Otherwise omit the basis quietly.
       * @throw std::logic_error if throw_if_no_match is true and for at least one atom with no matching (or empty) basis is found.
       * \note All instances of the same chemical element receive the same basis set.
       */
      BasisSet(const std::vector<Atom>& atoms,
               const std::vector<std::vector<Shell>>& element_bases,
               std::string name = "",
               const bool throw_if_no_match = false) : name_(std::move(name)) {
        // for each atom find the corresponding basis components
        for(auto a=0ul; a<atoms.size(); ++a) {

          auto Z = atoms[a].atomic_number;

          if (decltype(Z)(element_bases.size()) > Z && !element_bases.at(Z).empty()) {  // found? add shells in order
            for(auto s: element_bases.at(Z)) {
              this->push_back(std::move(s));
              this->back().move({{atoms[a].x, atoms[a].y, atoms[a].z}});
            } // shell loop
          }
          else if (throw_if_no_match) {  // not found? throw, if needed
            throw std::logic_error(std::string("did not find the basis for Z=") + std::to_string(Z) + " in the element_bases");
          }
        } // atom loop

        init();
      }

      const base_type& shells() const {
        return static_cast<const base_type&>(*this);
      }

      /// @return the number of shells in the basis
      std::size_t size() const {
        return shells().size();
      }

      /// @return iterator pointing to the first shell
      base_type::const_iterator begin() const {
        return shells().begin();
      }

      /// @return iterator pointing past the last shell
      base_type::const_iterator end() const {
        return shells().end();
      }

      using base_type::cbegin;
      using base_type::cend;
      using base_type::empty;

      /// @param[in] i index of the element to access
      /// @return const reference to `i`th element
      const Shell& operator[](std::size_t i) const {
        return shells()[i];
      }

      /// @param[in] i index of the element to access
      /// @return const reference to `i`th element
      const Shell& at(std::size_t i) const {
        return shells().at(i);
      }

      /// forces solid harmonics/Cartesian Gaussians
      /// @param solid if true, force all shells with L>1 to be solid harmonics, otherwise force all shells to Cartesian
      void set_pure(bool solid) {
        for(size_t s=0; s!=size(); ++s) {
          static_cast<base_type&>(*this)[s].contr[0].pure = solid;
        }
        init();
      }

      /// @return the number of basis functions in the basis; -1 if uninitialized
      long nbf() const {
        return nbf_;
      }
      /// @return the maximum number of primitives in a contracted Shell, i.e. maximum contraction length; 0 if uninitialized
      size_t max_nprim() const {
        return max_nprim_;
      }
      /// @return the maximum angular momentum of a contraction; -1 if uninitialized
      long max_l() const {
        return max_l_;
      }
      /// @return the map from shell index to index of the first basis function from this shell
      /// \note basis functions are ordered as shells, i.e. shell2bf[i] >= shell2bf[j] iff i >= j
      const std::vector<size_t>& shell2bf() const {
        return shell2bf_;
      }
      /// Computes the map from this object's shells to the corresponding atoms in \c atoms. If no atom matches the origin of a shell, it is mapped to -1.
      /// @note coordinates must match \em exactly , i.e. shell2atom[k] == l iff atoms[l].x == *this[k].O[0] && atoms[l].y == *this[k].O[1] &&  atoms[l].z == *this[k].O[2]
      /// @return the map from shell index to the atom in the list \c atoms that coincides with its origin;
      std::vector<long> shell2atom(const std::vector<Atom>& atoms) const {
        return shell2atom(*this, atoms, false);
      }
      /// Computes the map from \c atoms to the corresponding shells in this object. Coordinates are compared bit-wise (@sa BasisSet::shell2atom() )
      /// @return the map from atom index to the vector of shell indices whose origins conincide with the atom;
      /// @note this does not assume that \c shells are ordered in the order of atoms, as does BasisSet
      std::vector<std::vector<long>> atom2shell(const std::vector<Atom>& atoms) const {
        return atom2shell(atoms, *this);
      }

      /// Computes the map from \c shells to the corresponding atoms in \c atoms. Coordinates are compared bit-wise, i.e.
      /// shell2atom[k] == l iff atoms[l].x == *this[k].O[0] && atoms[l].y == *this[k].O[1] &&  atoms[l].z == *this[k].O[2]
      /// @param throw_if_no_match If true, and no atom matches the origin of a shell, throw a std::logic_error.
      ///        Otherwise such shells will be mapped to -1.
      /// @return the map from shell index to the atom in the list \c atoms that coincides with its origin;
      /// @throw std::logic_error if throw_if_no_match is true and for at least one shell no matching atom is found.
      static std::vector<long> shell2atom(const std::vector<Shell>& shells, const std::vector<Atom>& atoms, bool throw_if_no_match = false) {
        std::vector<long> result;
        result.reserve(shells.size());
        for(const auto& s: shells) {
          auto a = std::find_if(atoms.begin(), atoms.end(), [&s](const Atom& a){ return s.O[0] == a.x && s.O[1] == a.y && s.O[2] == a.z; } );
          const auto found_match = (a != atoms.end());
          if (throw_if_no_match && !found_match)
            throw std::logic_error("shell2atom: no matching atom found");
          result.push_back( found_match ? a - atoms.begin() : -1);
        }
        return result;
      }
      /// Computes the map from \c atoms to the corresponding shells in \c shells. Coordinates are compared bit-wise (@sa BasisSet::shell2atom() )
      /// @return the map from atom index to the vector of shell indices whose origins conincide with the atom;
      /// @note this does not assume that \c shells are ordered in the order of atoms, as does BasisSet
      static std::vector<std::vector<long>> atom2shell(const std::vector<Atom>& atoms, const std::vector<Shell>& shells) {
        std::vector<std::vector<long>> result;
        result.resize(atoms.size());
        size_t iatom = 0;
        for(const auto& a: atoms) {
          auto s = shells.begin();
          while (s != shells.end()) {
            s = std::find_if(s, shells.end(), [&a](const Shell& s){ return s.O[0] == a.x && s.O[1] == a.y && s.O[2] == a.z; } );
            if (s != shells.end()) {
              result[iatom].push_back( s - shells.begin());
              ++s;
            }
          }
          ++iatom;
        }
        return result;
      }

    private:
      std::string name_;
      long nbf_;
      size_t max_nprim_;
      int max_l_;
      std::vector<size_t> shell2bf_;

      void init() {
        nbf_ = libint2::nbf(*this);
        max_nprim_ = libint2::max_nprim(*this);
        max_l_ = libint2::max_l(*this);
        shell2bf_ = compute_shell2bf(*this);
      }

      struct canonicalizer {
          char operator()(char c) {
            char cc = ::tolower(c);
            switch (cc) {
              case '/': cc = 'I'; break;
            }
            return cc;
          }
      };

      static std::string canonicalize_name(const std::string& name) {
        auto result = name;
        std::transform(name.begin(), name.end(),
                       result.begin(), BasisSet::canonicalizer());
        return result;
      }

      // see http://gaussian.com/basissets/
      bool gaussian_cartesian_d_convention(const std::string& canonical_name) {
        // 3-21??g??, 4-31g??
        if (canonical_name.find("3-21")    == 0 ||
            canonical_name.find("4-31g")   == 0)
          return true;
        // 6-31??g?? but not 6-311 OR 6-31g()
        if (canonical_name.find("6-31") == 0 && canonical_name[4] != '1') {
          // to exclude 6-31??g() find the g, then check the next character
          auto g_pos = canonical_name.find('g');
          if (g_pos == std::string::npos) // wtf, I don't even know what this is, assume spherical d is OK
            return false;
          if (g_pos+1 == canonical_name.size()) // 6-31??g uses cartesian d
            return true;
          if (canonical_name[g_pos+1] == '*') // 6-31??g*? uses cartesian d
            return true;
        }
        return false;
      }

      /// decompose basis set name into components
      std::vector<std::string> decompose_name_into_components(std::string name) {
        std::vector<std::string> component_names;
        // aug-cc-pvxz* = cc-pvxz* + augmentation-... , except aug-cc-pvxz-cabs
        if ( (name.find("aug-cc-pv") == 0) && (name.find("cabs")==std::string::npos)  ) {
          std::string base_name = name.substr(4);
          component_names.push_back(base_name);
          component_names.push_back(std::string("augmentation-") + base_name);
        }
        else
          component_names.push_back(name);

        return component_names;
      }

      /** determines the path to the data directory, as follows:
       *       <ol>
       *         <li> specified by LIBINT_DATA_PATH environmental variable, if defined </li>
       *         <li> specified by DATADIR macro variable, if defined </li>
       *         <li> specified by SRCDATADIR macro variable, if defined </li>
       *         <li> hardwired to directory \c NONE/share/libint/2.8.0/basis </li>
       *       </ol>
       *  @throw std::system_error if the path is not valid, or cannot be determined
       *  @return valid path to the data directory
       */
      static std::string data_path() {
        std::string path;
        const char* data_path_env = getenv("LIBINT_DATA_PATH");
        if (data_path_env) {
          path = data_path_env;
        }
        else {
#if defined(DATADIR)
          path = std::string{DATADIR};
#elif defined(SRCDATADIR)
          path = std::string{SRCDATADIR};
#else
          path = std::string("NONE/share/libint/2.8.0");
#endif
        }
        // validate basis_path = path + "/basis"
        std::string basis_path = path + std::string("/basis");
        bool error = true;
        std::error_code ec;
        auto validate_basis_path = [&basis_path, &error, &ec]() -> void {
          if (not basis_path.empty()) {
            struct stat sb;
            error = (::stat(basis_path.c_str(), &sb) == -1);
            error = error || not S_ISDIR(sb.st_mode);
            if (error)
              ec = std::error_code(errno, std::generic_category());
          }
        };
        validate_basis_path();
        if (error) { // try without "/basis"
          basis_path = path;
          validate_basis_path();
        }
        if (error) {
          std::ostringstream oss; oss << "BasisSet::data_path(): path \"" << path << "{/basis}\" is not valid";
          throw std::system_error(ec, oss.str());
        }
        return basis_path;
      }

      /// converts fortran scientific-notation floats that use d/D instead of e/E in \c str
      /// @param[in,out] str string in which chars 'd' and 'D' are replaced with 'e' and 'E',
      ///                respectively
      static void fortran_dfloats_to_efloats(std::string& str) {
        for(auto& ch: str) {
          if (ch == 'd') ch = 'e';
          if (ch == 'D') ch = 'E';
        }
      }

    public:

      /** reads in all basis sets from a Gaussian94-formatted basis set file (see https://bse.pnl.gov/bse/portal)
       *  @param[in] file_dot_g94 file name
       *  @param[in] force_cartesian_d force use of Cartesian d shells, if true
       *  @param[in] locale_name specifies the locale to use
       *  @throw std::ios_base::failure if the path is not valid, or cannot be determined
       *  @throw std::logic_error if the contents of the file cannot be interpreted
       *  @return vector of basis sets for each element
       *  @warning the included library basis sets should be parsed using POSIX locale
       */
      static std::vector<std::vector<libint2::Shell>> read_g94_basis_library(std::string file_dot_g94,
                                                                             bool force_cartesian_d = false,
                                                                             bool throw_if_missing = true,
                                                                             std::string locale_name = std::string("POSIX")) {

        std::locale locale(locale_name.c_str());  // TODO omit c_str() with up-to-date stdlib
        std::vector<std::vector<libint2::Shell>> ref_shells(118); // 118 = number of chemical elements
        std::ifstream is(file_dot_g94);
        is.imbue(locale);

        if (is.good()) {
          if (libint2::verbose())
            libint2::verbose_stream() << "Will read basis set from " << file_dot_g94 << std::endl;

          std::string line, rest;

          auto LIBINT2_LINE_TO_STRINGSTREAM = [&](const std::string& line) {
            std::istringstream iss{line};
            iss.imbue(locale);
            return iss;
          };

          size_t Z;
          auto nextbasis = true, nextshell = false;
          bool first_element = true;
          // read lines till end
          do {
            // skipping empties and starting with '!' (the comment delimiter)
            if (line.empty() || line[0] == '!') continue;
            if (line == "****") {
              // old (EMSL) basis set exchange g94 format marks the beginning of data by ****
              // new (MolSSI) basis set exchange g94 format does not start with ****
              // so if found **** and still waiting for the first element, this is new g94 format, skip to next line
              if (first_element)
                continue;
              nextbasis = true;
              nextshell = false;
              continue;
            }
            if (nextbasis) {
              nextbasis = false;
              first_element = false;
              auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
              std::string elemsymbol;
              iss >> elemsymbol >> rest;

              bool found = false;
              for (const auto &e: libint2::chemistry::get_element_info()) {
                if (strcaseequal(e.symbol, elemsymbol)) {
                  Z = e.Z;
                  found = true;
                  break;
                }
              }
              if (not found) {
                std::ostringstream oss;
                oss << "in file " << file_dot_g94
                    << " found G94 basis set for element symbol \""
                    << elemsymbol << "\", not found in Periodic Table.";
                throw std::logic_error(oss.str());
              }

              nextshell = true;
              continue;
            }
            if (nextshell) {
              auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
              std::string amlabel;
              std::size_t nprim;
              iss >> amlabel >> nprim >> rest;
              if (amlabel != "SP" && amlabel != "sp") {
                assert(amlabel.size() == 1);
                auto l = Shell::am_symbol_to_l(amlabel[0]);
                svector<double> exps;
                svector<double> coeffs;
                for (decltype(nprim) p = 0; p != nprim; ++p) {
                  while (std::getline(is, line) && (line.empty() || line[0] == '!')) {}
                  fortran_dfloats_to_efloats(line);
                  auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                  double e, c;
                  iss >> e >> c;
                  exps.emplace_back(e);
                  coeffs.emplace_back(c);
                }
                auto pure = force_cartesian_d ? (l > 2) : (l > 1);
                ref_shells.at(Z).push_back(
                    libint2::Shell{
                        std::move(exps),
                        {
                            {l, pure, std::move(coeffs)}
                        },
                        {{0, 0, 0}}
                    }
                );
              } else { // split the SP shells
                svector<double> exps;
                svector<double> coeffs_s, coeffs_p;
                for (decltype(nprim) p = 0; p != nprim; ++p) {
                  while (std::getline(is, line) && (line.empty() || line[0] == '!')) {}
                  fortran_dfloats_to_efloats(line);
                  auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                  double e, c1, c2;
                  iss >> e >> c1 >> c2;
                  exps.emplace_back(e);
                  coeffs_s.emplace_back(c1);
                  coeffs_p.emplace_back(c2);
                }
                ref_shells.at(Z).push_back(
                    libint2::Shell{exps,
                                   {
                                       {0, false, coeffs_s}
                                   },
                                   {{0, 0, 0}}
                    }
                );
                ref_shells.at(Z).push_back(
                    libint2::Shell{std::move(exps),
                                   {
                                       {1, false, std::move(coeffs_p)}
                                   },
                                   {{0, 0, 0}}
                    }
                );
              }
            }
          } while (std::getline(is, line));

        }
        else {  // !is.good()
          if (throw_if_missing) {
            std::ostringstream oss;
            oss << "BasisSet::read_g94_basis_library(): could not open \"" << file_dot_g94 << "\"" << std::endl;
            throw std::ios_base::failure(oss.str());
          }
        }

        return ref_shells;
      }

      DEPRECATED static size_t nbf(const std::vector<libint2::Shell>& shells) {
        return libint2::nbf(shells);
      }

      DEPRECATED static size_t max_nprim(const std::vector<libint2::Shell>& shells) {
        return libint2::max_nprim(shells);
      }

      DEPRECATED static int max_l(const std::vector<libint2::Shell>& shells) {
        return libint2::max_l(shells);
      }

      static std::vector<size_t> compute_shell2bf(const std::vector<libint2::Shell>& shells) {
        std::vector<size_t> result;
        result.reserve(shells.size());

        size_t n = 0;
        for (auto shell: shells) {
          result.push_back(n);
          n += shell.size();
        }

        return result;
      }

     private:

      friend inline bool operator==(const BasisSet&, const BasisSet&);
      friend inline bool operator==(const BasisSet&, const base_type&);
      friend inline bool operator==(const base_type&, const BasisSet&);

  }; // BasisSet

  inline bool operator==(const BasisSet& bs1, const BasisSet& bs2) {
    return bs1.shells() == bs2.shells();
  }
  inline bool operator==(const BasisSet& bs1, const BasisSet::base_type& bs2) {
    return bs1.shells() == bs2;
  }
  inline bool operator==(const BasisSet::base_type& bs1, const BasisSet& bs2) {
    return bs1 == bs2.shells();
  }

} // namespace libint2

#endif /* _libint2_src_lib_libint_basis_h_ */
