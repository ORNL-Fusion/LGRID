//
//  vector3d.cpp
//  lgrid
//
//  Created by Cianciosa, Mark R. on 2/7/16.
//  Copyright © 2016 Cianciosa, Mark R. All rights reserved.
//

#include "vector3d.hpp"
#include <cmath>

//------------------------------------------------------------------------------
///  @brief @ref vector3d constructor.
///
///  Constructs a vector at a x, y, z position.
///
///  @param[in] x_ X position.
///  @param[in] y_ Y position.
///  @param[in] z_ Z position.
//------------------------------------------------------------------------------
vector3d::vector3d(const double x_, const double y_, const double z_) :
x(x_), y(y_), z(z_) {
}

//------------------------------------------------------------------------------
///  @brief Length of a @ref vector3d instance.
///
///  Length is defined as Sqrt(V.V).
///
///  @return Length of vector.
//------------------------------------------------------------------------------
const double vector3d::length() const {
    return sqrt(dot_product(*this, *this));
}

//------------------------------------------------------------------------------
///  @brief Normalize @ref vector3d instance.
///
///  The normalized vector is defined as V/|V|.
///
///  @return Normalized vector.
//------------------------------------------------------------------------------
const vector3d vector3d::normalize() const {
    return *this/this->length();
}

//------------------------------------------------------------------------------
///  @brief Takes the cross product of two vectors.
///
///  V1 X V2.
///
///  @param[in] v1 First vector.
///  @param[in] v2 Second vector.
///  @return Cross product of v1 and v2.
//------------------------------------------------------------------------------
const vector3d vector3d::cross_product(const vector3d &v1, const vector3d &v2) {
    return vector3d(v1.y*v2.z - v1.z*v2.y,
                    v1.z*v2.x - v1.x*v2.z,
                    v1.x*v2.y - v1.y*v2.x);
}

//------------------------------------------------------------------------------
///  @brief Takes the doc product of two vectors.
///
///  V1 . V2.
///
///  @param[in] v1 First vector.
///  @param[in] v2 Second vector.
///  @return Dot product of v1 and v2.
//------------------------------------------------------------------------------
const double vector3d::dot_product(const vector3d &v1, const vector3d &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//------------------------------------------------------------------------------
///  @brief Adds two vectors.
///
///  V1 + V2.
///
///  @param[in] v1 First vector.
///  @param[in] v2 Second vector.
///  @return Addition of v1 and v2.
//------------------------------------------------------------------------------
const vector3d operator+(const vector3d &v1, const vector3d &v2) {
    return vector3d(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

//------------------------------------------------------------------------------
///  @brief Subtracts two vectors.
///
///  V1 - V2.
///
///  @param[in] v1 First vector.
///  @param[in] v2 Second vector.
///  @return Subtraction of v1 and v2.
//------------------------------------------------------------------------------
const vector3d operator-(const vector3d &v1, const vector3d &v2) {
    return vector3d(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

//------------------------------------------------------------------------------
///  @brief Divides a vector by a scalar.
///
///  V1 / d.
///
///  @param[in] v1 Vector.
///  @param[in] d Scalar.
///  @return Divids all elements of v1 by d.
//------------------------------------------------------------------------------
const vector3d operator/(const vector3d &v1, const double d) {
    return vector3d(v1.x/d, v1.y/d, v1.z/d);
}
