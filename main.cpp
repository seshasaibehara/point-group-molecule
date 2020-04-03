#include "eigen-git-mirror/Eigen/Dense"
#include "headers.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
// Function to compare doubles
bool compare_floats(double a, double b, double absEpsilon, double relEpsilon)
{
    // Check if the numbers are really close -- needed when comparing numbers near zero.
    double diff{fabs(a - b)};
    if (diff <= absEpsilon)
        return true;
    // Otherwise fall back to Knuth's algorithm
    return diff <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * relEpsilon);
    //
}

molecule molecule::read_xyz(std::string file_name)
{
    std::ifstream input;
    input.open(file_name, std::ios::in);
    std::string line;

    std::vector<atom> coordinates_of_atoms;
    atom coordinates_of_atom;
    Vector3d coordinates;

    while (getline(input, line))
    {
        std::istringstream ss(line);
        std::string s;
        float x, y, z;

        ss >> s >> x >> y >> z;

        coordinates(0) = x;
        coordinates(1) = y;
        coordinates(2) = z;
        coordinates_of_atom.atomic_symbol = s;
        coordinates_of_atom.coordinate_vector = coordinates;
        coordinates_of_atoms.push_back(coordinates_of_atom);
    }
    input.close();
    molecule all_atoms(coordinates_of_atoms);
    return all_atoms;
}

std::tuple<int,int,Matrix3d> _eigen_stuff(const Matrix3d orthonormal_matrix){
    eigensolver ev(orthonormal_matrix,true);
    Vector3d evs,eigen_vector;
    Matrix3d eigen_vectors;

    evs = ev.eigenvalues().real();
    eigen_vectors = ev.eigenvectors().real();
    int one_counter = 0;
    int zero_counter = 0;
    for (int i = 0; i < 3; i++){
        if (compare_floats(evs(i,0),1,PREC,PREC)==true){
            one_counter = one_counter + 1;
        }
        if (compare_floats(evs(i,0),0,PREC,PREC)==true){
            zero_counter = zero_counter + 1;
        }
    }
    std::tuple<int,int,Matrix3d> tuple_to_return = std::make_tuple(one_counter,zero_counter,eigen_vectors);
    return tuple_to_return; 
}

Matrix3d _get_projection_operator(Vector3d position_vector){
    return position_vector*position_vector.transpose();
}

bool _is_planar(const molecule everything){ 
    std::vector<atom> coordinates_of_atoms = everything.all_atoms();
    std::vector<Vector3d> coords;
    Matrix3d check = Matrix3d::Zero();
    for (auto &x: coordinates_of_atoms){
        check = check + _get_projection_operator(x.coordinate_vector);
    }
    int zero_eigen_values = std::get<1>(_eigen_stuff(check));
    if (zero_eigen_values == 1){
        return true;
    }
    return false;
}

bool check_if_molecule_is_linear(const molecule everything){
    std::vector<atom> coordinates_of_atoms = everything.all_atoms();
    std::vector<Vector3d> coords;
    Matrix3d check = Matrix3d::Zero();
    for (auto &x: coordinates_of_atoms){
        check = check + _get_projection_operator(x.coordinate_vector);
    }
    int zero_eigen_values = std::get<1>(_eigen_stuff(check));
    if (zero_eigen_values == 2){
        return true;
    }
    return false;
}

Matrix3d compute_rotation_matrix(const std::vector<atom>& coordinates_of_atoms){
    Matrix3d check = Eigen::Matrix3d::Zero();
    for (auto &x: coordinates_of_atoms){
        check = check + _get_projection_operator(x.coordinate_vector);
    }
    Matrix3d eigen_vector_space = std::get<2>(_eigen_stuff(check));
   
    Matrix3d rotation_1 = eigen_vector_space.transpose();
    Matrix3d rotation_matrix, rotation_2, rotation_3, rotation_4;
    rotation_2 << 0,1,0,
               0,0,1,
               1,0,0;
    rotation_3 << 1,0,0,
               0,0,1,
               0,1,0;
    rotation_4 << 1,0,0,
               0,1,0,
               0,0,1;
        
    int counter_x = 0;
    int counter_y = 0;
    int counter_z = 0;
    for (std::size_t i=0; i<coordinates_of_atoms.size(); ++i){
        Vector3d check_0 = rotation_1*coordinates_of_atoms[i].coordinate_vector;
        if (compare_floats(check_0(0),0,PREC,PREC)==true){
            ++counter_x;
        }
        if (compare_floats(check_0(1),0,PREC,PREC)==true){
            ++counter_y;
        }
        if (compare_floats(check_0(2),0,PREC,PREC)==true){
            ++counter_z;
        }
    }
    if (counter_x == coordinates_of_atoms.size()){
        rotation_matrix = rotation_2*rotation_1;
    }
    if (counter_y == coordinates_of_atoms.size()){
        rotation_matrix = rotation_3*rotation_1;
    }
    if (counter_z == coordinates_of_atoms.size()){
        rotation_matrix = rotation_4*rotation_1;
    }
    return rotation_matrix;
}

std::vector<atom> compute_transformed_coordinates(const std::vector<atom>& coordinates_of_atoms, Matrix3d rotation_matrix){
    std::vector<atom> transformed_coordinates;
    atom transformed_coordinate;
    for (auto& x : coordinates_of_atoms){
        transformed_coordinate.coordinate_vector = rotation_matrix * x.coordinate_vector;
        transformed_coordinate.atomic_symbol = x.atomic_symbol;
        transformed_coordinates.push_back(transformed_coordinate);
    }
    return transformed_coordinates;
}
  
std::string label_sym_ops(const Matrix3d orthonormal_matrix){
    double trace = orthonormal_matrix.trace();
    double det = orthonormal_matrix.determinant();
    int eigen_values = std::get<0>(_eigen_stuff(orthonormal_matrix));
    std::string label;
  
    if (compare_floats(det, 1, PREC, PREC) == true && compare_floats(trace, 3, PREC, PREC) == true){
        label = "Identity";
    }

    if (compare_floats(det, -1, PREC, PREC) == true && compare_floats(trace, -3, PREC, PREC) == true){
        label = "Inverse";
    }

    if (compare_floats(det, 1, PREC, PREC) == true && eigen_values == 1){
        label = "Rotation";
    }

    if (compare_floats(det, -1, PREC, PREC) == true && compare_floats(trace, -3, PREC, PREC) == false && eigen_values == 0){
        label = "Improper Rotation";
    }

    if (compare_floats(det, -1, PREC, PREC) == true && eigen_values == 2){
        label = "Mirror";
    }
    return label;
}

//use STL
bool check_if(std::vector<Matrix3d>& sym_ops, const Matrix3d possible_sym_op){
    auto it_k = std::find_if(sym_ops.begin(), sym_ops.end(), [&](const Eigen::Matrix3d& sym){return sym.isApprox(possible_sym_op,PREC);});
    
    if (it_k != sym_ops.end()){
        return true;
    }
    return false;
}

bool compare_vectors(Vector3d v1, Vector3d v2){
    for (int i = 0; i < 3; i++){
        if (compare_floats(v1(i), v2(i), PREC, PREC) == false){
            return false;
        }
    }
    return true;
}

bool checking_structure_equivalency(const std::vector<atom>& coordinates_before_symmetry, const std::vector<atom>& coordinates_after_symmetry){
    for (auto& x : coordinates_after_symmetry){
        int counter = 0;
        for (auto& y : coordinates_before_symmetry)
        {
            if (x.atomic_symbol == y.atomic_symbol && compare_vectors(x.coordinate_vector, y.coordinate_vector) == true)
            {
                counter = counter + 1;
            }
        }
        if (counter == 0)
        {
            return false;
        }
    }
    return true;
}

bool is_sym_op(Matrix3d& possible_sym_op, std::vector<atom>& coordinates_of_atoms)
{
    if (compare_floats(possible_sym_op.determinant(), 1, PREC, PREC) == true ||
        compare_floats(possible_sym_op.determinant(), -1, PREC, PREC) == true)
    {
        Matrix3d identity_matrix = possible_sym_op.transpose() * possible_sym_op;

        if (identity_matrix.isIdentity(PREC) == true)
        {
            std::vector<atom> transformed_coordinates =
                compute_transformed_coordinates(coordinates_of_atoms, possible_sym_op);

            if (checking_structure_equivalency(coordinates_of_atoms, transformed_coordinates) == true)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::vector<int>> all_permutations(std::vector<int>& starting_permutation)
{

    std::vector<std::vector<int>> all_the_things;
    std::sort(starting_permutation.begin(), starting_permutation.end());
    do
    {
        all_the_things.push_back(starting_permutation);
    } while (next_permutation(starting_permutation.begin(), starting_permutation.end()));

    return all_the_things;
}

std::vector<int> first_combination(int size)
{

    assert(size > 0);
    std::vector<int> combination;
    for (int i = 0; i < size; ++i)
    {
        combination.push_back(i);
    }
    return combination;
}

bool next_combination(std::vector<int>& v, int m, int n)
{

    for (int i = m - 1; i >= 0; i--)
    {
        if (v[i] + m - i < n)
        {
            v[i]++;
            for (int j = i + 1; j < m; j++)
            {
                v[j] = v[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

Matrix2d make_2dmatrix(std::vector<int>& curr_comb, std::vector<atom>& vector_of_similar_type_atoms)
{

    Matrix2d matrix_after_sym_op;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            matrix_after_sym_op(i, j) = vector_of_similar_type_atoms[curr_comb[j]].coordinate_vector(i);
        }
    }
    return matrix_after_sym_op;
}

Matrix3d make_3dmatrix(std::vector<int>& curr_comb, std::vector<atom>& vector_of_similar_type_atoms)
{

    Matrix3d matrix_after_sym_op;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            matrix_after_sym_op(i, j) = vector_of_similar_type_atoms[curr_comb[j]].coordinate_vector(i);
        }
    }
    return matrix_after_sym_op;
}

MatrixXd make_xdmatrix(std::vector<int>& curr_comb, std::vector<atom>& vector_of_similar_type_atoms)
{

    MatrixXd matrix_after_sym_op(3, 2);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            matrix_after_sym_op(i, j) = vector_of_similar_type_atoms[curr_comb[j]].coordinate_vector(i);
        }
    }
    return matrix_after_sym_op;
}

void print_vector(const std::vector<int>& v)
{

    for (auto& x : v)
    {
        std::cout << x << "  ";
    }
    std::cout << std::endl;
}

std::vector<sym_op_class> calculate_sym_ops(molecule everything){
    std::vector<std::string> atom_types = everything.atom_types();        
    std::vector<sym_op_class> point_group;
    sym_op_class temp_sym_op;
    std::vector<Matrix3d> sym_ops;
    std::vector<atom> coordinates_of_atoms = everything.all_atoms();
    for (auto& x: atom_types)
    {
        std::vector<atom> vector_of_similar_type_atoms = everything.atoms_of_type(x);
        if (vector_of_similar_type_atoms.size() == 1)
        {
            continue;
        }

        if (vector_of_similar_type_atoms.size() == 2)
        {
            MatrixXd random_matrix(3, 2);
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    random_matrix(i, j) = vector_of_similar_type_atoms[vector_of_similar_type_atoms.size() - j - 1].coordinate_vector(i);
                }
            }
            int max_index = vector_of_similar_type_atoms.size();
            int size = 2;
            std::vector<int> curr_combo = first_combination(size);
            do
            {
                std::vector<std::vector<int>> all_the_things = all_permutations(curr_combo);
                for (auto& x : all_the_things)
                {
                    MatrixXd matrix_after_sym_op = make_xdmatrix(x, vector_of_similar_type_atoms);
                    Matrix3d inverse_matrix = random_matrix * random_matrix.transpose();
                    Matrix3d matrix = matrix_after_sym_op * random_matrix.transpose();
                    Matrix3d possible_sym_op = matrix * inverse_matrix.inverse();
                    if (is_sym_op(possible_sym_op, coordinates_of_atoms) == true &&
                      check_if(sym_ops,possible_sym_op) == false)
                    {
                        sym_ops.push_back(possible_sym_op);
                    }
                }
            } while (next_combination(curr_combo, size, max_index));
        }

        Matrix3d random_matrix;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                random_matrix(i, j) = vector_of_similar_type_atoms[vector_of_similar_type_atoms.size() - j - 1].coordinate_vector(i);
            }
        }
        int max_index = vector_of_similar_type_atoms.size();
        int size = 3;
        std::vector<int> curr_combo = first_combination(size);
        do
        {
            std::vector<std::vector<int>> all_the_things = all_permutations(curr_combo);
            for (auto& x : all_the_things)
            {
                Matrix3d matrix_after_sym_op = make_3dmatrix(x, vector_of_similar_type_atoms);
                Matrix3d possible_sym_op = matrix_after_sym_op * random_matrix.inverse();
                if (is_sym_op(possible_sym_op, coordinates_of_atoms) == true && check_if(sym_ops, possible_sym_op) == false)
                {
                    sym_ops.push_back(possible_sym_op);
                }
            }
        } while (next_combination(curr_combo, size, max_index));
    }

    for (auto& x : sym_ops)
    {
        temp_sym_op.coord_matrix = x;
        temp_sym_op.label = label_sym_ops(x);
        point_group.push_back(temp_sym_op);
    }

    if (point_group.size() == 0)
    {
        temp_sym_op.coord_matrix = Matrix3d::Identity();
        temp_sym_op.label = label_sym_ops(Matrix3d::Identity());
        point_group.push_back(temp_sym_op);
    }
    return point_group;
}

std::vector<sym_op_class> calculate_sym_ops_planar(molecule everything)
{

    std::vector<sym_op_class> point_group_planar;
    sym_op_class temp_sym_op;
    std::vector<Matrix3d> sym_ops;
    std::vector<Matrix3d> possible_sym_ops;
    Matrix3d temp_matrix;
    std::vector<atom> coordinates_of_atoms = everything.all_atoms();

    Matrix3d rotation_matrix = compute_rotation_matrix(coordinates_of_atoms);
    std::vector<atom> transformed_coordinates = compute_transformed_coordinates(coordinates_of_atoms, rotation_matrix);
   
    molecule all_atoms_planar(transformed_coordinates);
    std::vector<std::string> atom_types_planar = all_atoms_planar.atom_types();
    
    for (auto& x: atom_types_planar)
    {
        std::vector<atom> vector_of_similar_type_atoms = all_atoms_planar.atoms_of_type(x);
        if (vector_of_similar_type_atoms.size() == 1)
        {
            continue;
        }

        Matrix2d random_matrix;
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                random_matrix(i, j) = vector_of_similar_type_atoms[vector_of_similar_type_atoms.size() - j - 1].coordinate_vector(i);
            }
        }
        int max_index = vector_of_similar_type_atoms.size();
        int size = 2;
        std::vector<int> curr_combo = first_combination(size);

        do
        {
            std::vector<std::vector<int>> all_the_things = all_permutations(curr_combo);
            for (auto& x : all_the_things)
            {
                Matrix2d matrix_after_sym_op = make_2dmatrix(x, vector_of_similar_type_atoms);
                Matrix2d possible_sym_op = matrix_after_sym_op * random_matrix.inverse();
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        temp_matrix(i, j) = possible_sym_op(i, j);
                    }
                    temp_matrix(i, 2) = 0;
                    temp_matrix(2, i) = 0;
                }
                temp_matrix(2, 2) = 1;
                possible_sym_ops.push_back(temp_matrix);
                temp_matrix(2, 2) = -1;
                possible_sym_ops.push_back(temp_matrix);
            }
        } while (next_combination(curr_combo, size, max_index));
    }

    for (auto& possible_sym_op : possible_sym_ops)
    {
        if (is_sym_op(possible_sym_op, transformed_coordinates) == true && check_if(sym_ops, possible_sym_op) == false)
        {
            sym_ops.push_back(possible_sym_op);
        }
    }

    for (auto& x : sym_ops)
    {
        temp_sym_op.coord_matrix = rotation_matrix.inverse() * x * rotation_matrix;
        temp_sym_op.label = label_sym_ops(x);
        point_group_planar.push_back(temp_sym_op);
    }
    sym_op_class temp_sym;
    if (point_group_planar.size() == 0)
    {
        temp_sym.coord_matrix = Matrix3d::Identity();
        temp_sym.label = label_sym_ops(Matrix3d::Identity());
        point_group_planar.push_back(temp_sym_op);
    }

    return point_group_planar;
}

void print_molecule_point_group(std::vector<sym_op_class> point_group)
{

    for (auto& x : point_group)
    {
        std::cout << "---------" << std::endl;
        std::cout << x.label << std::endl;
        std::cout << x.coord_matrix.format(CleanFmt) << std::endl;
    }

    std::cout << "---------" << std::endl;
    std::cout << "Number of point group operations: " << point_group.size() << std::endl;
}

bool check_for_inversion(const molecule everything){

    std::vector<atom> coordinates_of_atoms = everything.all_atoms();
    std::vector<atom> after_inversion;
    Vector3d zero_vector;
    atom temp_atom;
    zero_vector << 0, 0, 0;

    for (auto& x : coordinates_of_atoms)
    {
        temp_atom.coordinate_vector = zero_vector - x.coordinate_vector;
        temp_atom.atomic_symbol = x.atomic_symbol;
        after_inversion.push_back(temp_atom);
    }

    if (checking_structure_equivalency(coordinates_of_atoms, after_inversion) == true)
    {
        return true;
    }

    return false;
}

int main(int argc, char* argv[])
{


    if (argc == 1)
    {
        std::cerr << "No molecule specified" << std::endl;
        return 1;
    }
    
    molecule everything = molecule::read_xyz(argv[1]); 
    if (check_if_molecule_is_linear(everything) == true)
    {
        if (check_for_inversion(everything) == true)
        {
            std::cout << "The molecule belongs to D\u221Eh point group" << std::endl;
            return 0;
        }
        std::cout << "The molecule belongs to C\u221Ev point group" << std::endl;
        return 0;
    }

    if (_is_planar(everything) == true &&
        check_if_molecule_is_linear(everything) == false)
    {
        std::vector<sym_op_class> point_group = calculate_sym_ops_planar(everything);
        print_molecule_point_group(point_group);
        return 0;
    }

    std::vector<sym_op_class> point_group = calculate_sym_ops(everything);
    print_molecule_point_group(point_group);
    return 0;
}
