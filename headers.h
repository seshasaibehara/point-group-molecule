#ifndef HEADERS_H
#define HEADERS_H
#include "eigen-git-mirror/Eigen/Dense"
#include "eigen-git-mirror/Eigen/Eigenvalues"
#include "eigen-git-mirror/Eigen/src/Core/Matrix.h"
#include <unordered_map>
#include <vector>

#define PREC 1e-3

typedef Eigen::Matrix2d Matrix2d;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;

typedef Eigen::SelfAdjointEigenSolver<Matrix3d> selfadjoint_eigensolver;
typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> PermutationMatrix;
typedef Eigen::EigenSolver<Matrix3d> eigensolver;
Eigen::IOFormat CleanFmt(7, 0, ", ", "\n", "[", "]");

class sym_op_class
{
public:
    Matrix3d coord_matrix;
    std::string label;
};

struct atom
{
    std::string atomic_symbol;
    Vector3d coordinate_vector;
};

class molecule
{
public:
    molecule(const std::vector<atom>& molecule_atoms) : m_all_atoms(molecule_atoms)
    {
        for (const atom& current_atom : molecule_atoms)
        {
            std::string current_type = current_atom.atomic_symbol;
            if (m_binned_atoms.count(current_type) == 0)
            {
                m_binned_atoms[current_type] = std::vector<atom>();
            }

            m_binned_atoms[current_type].push_back(current_atom);
        }
            
    shift_to_center_of_mass(&m_all_atoms);
    return;
    }

    static molecule read_xyz(const std::string filename);

    std::vector<atom> atoms_of_type(const std::string& atom_type) const
    {
        return m_binned_atoms.at(atom_type);
    }

    std::vector<atom> all_atoms() const
    {
        return m_all_atoms;
    }

    std::vector<std::string> atom_types() const
    {
        std::vector<std::string> all_keys;
        for(const auto& pair : m_binned_atoms)
        {
            all_keys.push_back(pair.first);
        }
        return all_keys;
    }

private:
    std::vector<atom> m_all_atoms;
    std::unordered_map<std::string, std::vector<atom>> m_binned_atoms;
    Vector3d compute_center_of_mass(const std::vector<atom>& coordinates_of_atoms){
    // computing centre of mass and returning the adjusted coordinates-----
        Vector3d center_of_mass(0, 0, 0);
        Vector3d test(0, 0, 0);

        for (auto& x : coordinates_of_atoms){
        for (int j = 0; j < 3; ++j){
            center_of_mass(j) = center_of_mass(j) + x.coordinate_vector(j);
        }
    }
    center_of_mass = center_of_mass / coordinates_of_atoms.size();
    return center_of_mass;
    }
    void shift_to_center_of_mass(std::vector<atom>* coordinates_of_atoms_ptr){
        auto& coordinates_of_atoms = *coordinates_of_atoms_ptr;
        auto center_of_mass = compute_center_of_mass(coordinates_of_atoms);

        for (auto& y : coordinates_of_atoms){
            for (int j = 0; j < 3; j++){
                y.coordinate_vector(j) = y.coordinate_vector(j) - center_of_mass(j);
            }
        }
    return;
    }
};






#endif
